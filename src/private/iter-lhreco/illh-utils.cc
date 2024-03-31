/**
 * multi-llh
 *
 * @version $Id: $
 *
 * @date: $Date: $
 *
 * @author Juan Carlos Diaz-Velez <juan.diazvelez@alumnos.udg.mx>
 *
*/

#include <iter-lhreco-proj/illh-utils.h>

namespace fs = boost::filesystem; 
using namespace std;
using namespace boost::numeric::ublas;
using boost::format;



/*
 * Sum over all pixels
 *
 */
double 
illh::Sum(const SkyMap& map,unsigned int npix) 
{
    double sumval = 0.;
    for (unsigned int i=0; i < npix; i++) 
    {
         sumval += map[i];
    }
    return sumval;
}


#if __cplusplus > 199711L
void illh::fluctuate(std::vector<SkyMapPtr>& Nmap, boost::mt19937 rng)
{
    log_info("adding random fluctuations... " );

    // create a generator 
    //Mersenne Twister generator 

    typedef boost::variate_generator< 
                    boost::mt19937, boost::poisson_distribution<> 
                      > rnd_poisson_t; 
    
    for(std::vector<SkyMapPtr>::iterator it = Nmap.begin(); it != Nmap.end(); ++it) { 
            SkyMapPtr lmap = *it;
            unsigned nside = lmap->Nside();
            unsigned npix = 12*nside*nside; 
            for (unsigned i=0;i<npix;++i) 
            {
                    double lambda = (*lmap)[i];
                    if (lambda > 0)
                    {
                       rnd_poisson_t rnd_poisson(rng, boost::poisson_distribution<>( lambda )); 
                       (*lmap)[i] = rnd_poisson();
                    }
            }
    }
}



/*
 * Generate isotropic maps with Poisson noise
 *
 */
void illh::isotropic(std::vector<SkyMapPtr>& Nmap, boost::mt19937 rng)
{
    log_info("generating isotropic maps... " );

    // create a generator 
    //Mersenne Twister generator 

    typedef boost::variate_generator< 
                    boost::mt19937, boost::poisson_distribution<> 
                      > rnd_poisson_t; 
    
    unsigned nside = Nmap[0]->Nside();
    unsigned npix = 12*nside*nside; 

    for (unsigned i=0;i<npix;++i) 
    {
            double lambda= 0.;
            unsigned count = 0;

            for(std::vector<SkyMapPtr>::iterator it = Nmap.begin(); it != Nmap.end(); ++it) 
            {
                    SkyMapPtr lmap = *it;
                    lambda += (*lmap)[i];
                    ++count;
            }
            if (!(lambda > 0)) 
                    continue;

            rnd_poisson_t rnd_poisson(rng, boost::poisson_distribution<>( lambda/count )); 
            for(std::vector<SkyMapPtr>::iterator it = Nmap.begin(); it != Nmap.end(); ++it) 
            { 
                    SkyMapPtr lmap = *it; 
                    (*lmap)[i] = rnd_poisson(); 
            }

    }
}
#endif



/*
 * loadMap - read local maps and generate vector of maps bined in 
 * siderial time steps
 */
void illh::loadMap(
    std::vector<SkyMapPtr>& Nmap,
    unsigned timeidxMin,
    unsigned timeidxMax, 
    unsigned nTimesteps, 
    unsigned nsideIn, 
    unsigned nsideOut, 
    std::string prefix,
    std::string suffix )
{
    fitshandle handle;
    std::string cuts;
    std::string coords;
    std::string timesys;
    double startMJD = -1.;
    double stopMJD = -1.;
    double totDur = 0.;
    double nEventsFio = 0.;
    int nTimeBins = -1;



    // Iterate over time bins and read local maps
    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
    {
        //double nEvents = 0.0; 

        // Read header info and then map
        std::stringstream input;
        input << prefix; 
        input << setfill('0')<<setw(3)<<timeidx;
        input << suffix;
        log_info("reading file : " << input.str() );
        handle.open(input.str().c_str());

        std::string cutsTemp = "";
        if (handle.key_present("CUTS")){
          handle.get_key("CUTS", cutsTemp);
        }
        string coordsTemp = "";  
        if (handle.key_present("COORDS")){
          handle.get_key("COORDS", coordsTemp);
        }
        double nEventsTemp = 0;
        if (handle.key_present("NEVENTS")){
          handle.get_key("NEVENTS",  nEventsTemp);
        }
        double startMJDTemp = 0;
        if (handle.key_present("STARTMJD")){
          handle.get_key("STARTMJD", startMJDTemp);
        }
        double stopMJDTemp = 0;
        if (handle.key_present("STOPMJD")){
          handle.get_key("STOPMJD", stopMJDTemp);
        }
        double totDurTemp = 0;
        if (handle.key_present("TOTDUR")){
          handle.get_key("TOTDUR",   totDurTemp);
        }
        std::string timesysTemp = "";
        if (handle.key_present("TIMESYS")){
          handle.get_key("TIMESYS", timesysTemp);
        }
        int bin = -1;
        if (handle.key_present("TIMEBIN")){
            handle.get_key("TIMEBIN",   bin);
        }
        int nTimeBinsTemp = -1;
        if (handle.key_present("NBINS")) {
            handle.get_key("NBINS",  nTimeBinsTemp);
        }

        handle.goto_hdu(2);
        SkyMapPtr locMap(new SkyMap);
        read_Healpix_map_from_fits(handle, *locMap, 1);
        handle.close();

        if (timeidx == timeidxMin) {
            nsideIn = locMap->Nside();
            log_info("Old Nside: " << nsideIn << ", New Nside: " << nsideOut);
            cuts = cutsTemp;
            coords = coordsTemp;
            timesys = timesysTemp;
            startMJD = startMJDTemp;
            stopMJD = stopMJDTemp;
            nTimeBins = nTimeBinsTemp;
        }

        totDur += totDurTemp;
        nEventsFio += nEventsTemp;

        if (cuts != cutsTemp) {
            log_fatal("Cuts are different");
        }
        if (coords != coordsTemp && !coordsTemp.empty()) {
            log_fatal("Coordinate systems are different");
        }
        if (timesys != timesysTemp && !timesysTemp.empty()) {
            log_fatal("Time systems are different");
        }
        if (nTimeBins != nTimeBinsTemp) {
            log_fatal("Number of time bins are different");
        }
        if (bin != timeidx && bin > 0) {
            log_fatal("Time Bin ID does not match Time Index of Loop. Maps may be out of order.");
        }

        if ( nsideIn != nsideOut) {
            SkyMapPtr locMapDegr(new SkyMap); 
            locMapDegr->SetNside(nsideOut, RING);     
            locMapDegr->fill(0.);
            for (int i = 0; i < 12*nsideIn*nsideIn; i++){
                pointing pt        = locMap->pix2ang(i);
                int j              = locMapDegr->ang2pix(pt);
                (*locMapDegr)[j] += (*locMap)[i];
            }
            Nmap.push_back(locMapDegr);
        } else {
            Nmap.push_back(locMap);
        }
    }

}

/*
 * save_iter - save results of iteration to FITS file
 */
void illh::save_iter(
                std::string foldername,
                std::vector<double> norm,
                SkyMap& Emap,
                std::string detector,
                unsigned nsideOut, 
                unsigned nTimesteps, 
                unsigned iteration)
{
    //  write N_tau^(i)
    stringstream normName;
    normName << foldername << boost::format("/norm_%s_%d_%d_iteration%02d.dat") % detector % nsideOut % nTimesteps % iteration;
    std::ofstream fileout(normName.str().c_str());
    for(std::vector<double>::iterator it = norm.begin(); it != norm.end(); ++it) {
      fileout << *it << "\n";
    }
    fileout.close();

    // write A_i^(n)
    fitshandle fitsOut; 
    stringstream expmapname;
    expmapname << foldername << boost::format("/exposure_%s_%d_%d_iteration%02d.fits.gz") % detector % nsideOut % nTimesteps % iteration;
    if (fs::exists(expmapname.str()) ) { 
            fs::remove( expmapname.str() ); 
    } 
    fitsOut.create(expmapname.str().c_str()); 
    write_Healpix_map_to_fits(fitsOut, Emap, MyDTYPE);
    fitsOut.close();
}

/*
 * Rotate from local detector coordinates to J2000 Equatorial reference frame
 * @param: i - pixel in Healpix map
 * @param: timeidx - time bin
 * @param: lat - detector latitude
 * @param: lon - detector longitude
 * @param: nTimesteps - number of time bins
 * @param: CRmap - healpix map (class needed to convert between directions and pixels
 * 
 * @return rotated pixel id
 */
unsigned 
illh::eq2log_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, SkyMap& CRmap)
{ 
        vec3 v =  CRmap.pix2vec(i); 
        double clat = cos(lat); 
        double slat = sin(lat);
        double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon; 
        double cb = cos(beta); 
        double sb = sin(beta);
        
        // rotation matrix
        matrix<double> rot (3, 3);

        rot(0,0) = cb*slat;
        rot(0,1) = sb*slat;
        rot(0,2) = -clat;

        rot(1,0) = -sb;
        rot(1,1) = cb;
        rot(1,2) = 0;

        rot(2,0) = cb*clat;
        rot(2,1) = clat*sb;
        rot(2,2) = slat;


        // rotation from local frame to Equatorial (ra,dec)
        vec3 vp;
        vp.x = rot(0,0)*v.x+rot(0,1)*v.y+rot(0,2)*v.z;
        vp.y = rot(1,0)*v.x+rot(1,1)*v.y+rot(1,2)*v.z;
        vp.z = rot(2,0)*v.x+rot(2,1)*v.y+rot(2,2)*v.z;
        
        return CRmap.vec2pix(vp); // local pixel
}

/*
 * Rotate from J2000 Equatorial reference frame to local detector coordinates 
 * @param: i - pixel in Healpix map
 * @param: timeidx - time bin
 * @param: lat - detector latitude
 * @param: lon - detector longitude
 * @param: nTimesteps - number of time bins
 * @param: CRmap - healpix map (class needed to convert between directions and pixels
 * 
 * @return rotated pixel id
 */
unsigned 
illh::loc2eq_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, SkyMap& CRmap)
{ 

        vec3 v =  CRmap.pix2vec(i); 
        double clat = cos(lat); 
        double slat = sin(lat);
        double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon; 
        double cb = cos(beta); 
        double sb = sin(beta);
        
        // rotation matrix
        matrix<double> rot (3, 3);

        rot(0,0) = cb*slat;
        rot(0,1) = sb*slat;
        rot(0,2) = -clat;

        rot(1,0) = -sb;
        rot(1,1) = cb;
        rot(1,2) = 0;

        rot(2,0) = cb*clat;
        rot(2,1) = clat*sb;
        rot(2,2) = slat;

        // rotation from local frame to Equatorial (ra,dec)
        vec3 vp;
        vp.x = rot(0,0)*v.x+rot(1,0)*v.y+rot(2,0)*v.z;
        vp.y = rot(0,1)*v.x+rot(1,1)*v.y+rot(2,1)*v.z;
        vp.z = rot(0,2)*v.x+rot(1,2)*v.y+rot(2,2)*v.z;
        
        return CRmap.vec2pix(vp); // local pixel
}


