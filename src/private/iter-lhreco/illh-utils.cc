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
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <iostream>
#include <vector>
#include <memory>

// CFITSIO
#include <fitsio.h>



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


SkyMapPtr read_map_column(
    fitsfile *fptr,
    const std::string &colname)
{
    int status = 0;
    int colnum = 0;

    // Look up column number by name (case-insensitive, CFITSIO convention)
    fits_get_colnum(fptr, CASEINSEN, (char*)colname.c_str(), &colnum, &status);
    if (status) {
        fits_report_error(stderr, status);
        throw std::runtime_error("Column " + colname + " not found.");
    }

    long nrows = 0;
    fits_get_num_rows(fptr, &nrows, &status);

    if (status) {
        fits_report_error(stderr, status);
        throw std::runtime_error("Error reading table size.");
    }

    // Construct Healpix_Map
    SkyMapPtr map(new SkyMap);
    long npix = nrows;
    int nside = map->npix2nside(npix);
    map->SetNside(nside, RING);  // or NESTED

    std::vector<double> coldata(npix);

    fits_read_col(fptr, TDOUBLE, colnum,
                  1, 1, npix,
                  nullptr,
                  coldata.data(),
                  nullptr,
                  &status);

    if (status) {
        fits_report_error(stderr, status);
        throw std::runtime_error("Error reading column " + colname);
    }

    for (int p = 0; p < npix; ++p)
        (*map)[p] = coldata[p];

    return map;
}

/*
 * loadMap - read local maps and generate vector of maps bined in 
 * siderial time steps
 */
void illh::loadMap(
    std::vector<SkyMapPtr>& Nmap,
    unsigned nTimesteps, 
    unsigned nsideIn, 
    unsigned nsideOut, 
    std::string filename,
    bool weights)
{
    int  hdu(2);
    if (weights) { 
        log_info("reading weighted maps");
        hdu=3;
    }
    fitsfile *fptr;
    int status = 0;

    fits_open_file(&fptr, filename.c_str(), READONLY, &status);
    if (status) { 
        fits_report_error(stderr, status); 
        throw std::runtime_error("Problem reading " + filename+ ".");
    }

    // Move to FITS table HDU
    log_info("reading file "<<filename << " status:"<<status);
    std::vector<SkyMapPtr> inmaps = read_group_hdu(fptr, hdu);

    // Iterate over time bins and read local maps


    for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
    { 
        std::string colname = "TBIN" + std::to_string(timeidx+1);

        //std::cout << "Reading column: " << colname << std::endl;
        SkyMapPtr locMapIn = inmaps[timeidx];
        //SkyMapPtr locMapIn = read_map_column(fptr, colname); 


        if ( nsideIn != nsideOut) {
            SkyMapPtr locMap(new SkyMap);
            locMap->SetNside(nsideOut, RING);     
            locMap->fill(0.);
            for (int i = 0; i < 12*nsideIn*nsideIn; i++){
                 pointing pt  = locMapIn->pix2ang(i);
                 int j        = locMap->ang2pix(pt);
                (*locMap)[j] += (*locMapIn)[i];
            }
            Nmap.push_back(locMap);
        } else {
            Nmap.push_back(locMapIn);
        }
    }
    fits_close_file(fptr, &status);

}

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
    std::string suffix,
    bool weights)
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



    if (weights) { 
        log_info("reading weighted maps");
    }
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
        read_Healpix_map_from_fits(handle, *locMap, weights ? 2 : 1);
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
    fitsOut.add_comment(std::string("LLH Exposure Map"));
    fitsOut.set_key("COORDS" , std::string("L"), "Local coordinate system");
    fitsOut.set_key("TTYPE1", std::string("acceptance"), "Normalized detector acceptance");
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
illh::eq2loc_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, const SkyMap& CRmap)
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
illh::loc2eq_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, const SkyMap& CRmap)
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

void illh::mlh_iteration(
        unsigned int pixmin, 
        unsigned int pixmax, 
        unsigned int nTimesteps,
        const std::vector<DetectorPtr> detectors,
        SkyMap &diffCRmap,
        SkyMap &CRmap,
        SkyMap &bkgMap,
		std::mutex &pmutex)
{ 

    for (unsigned int i=pixmin; i<pixmax;i++) 
    {
        double nEvents = 0.0; 
        double nBkg = 0.0; 
        for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
        { 
            std::vector<DetectorPtr>::const_iterator det_it;
            for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
            {
                    DetectorPtr det = *det_it;
                    int j = illh::eq2loc_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
                    if (det->FOV[j] && ((*det->Emap)[j] > 0.0)) { 
                        nEvents += (*det->Nmap[timeidx])[j]; 
                        nBkg += det->norm[timeidx]*(*det->Emap)[j]; 
                    } 
            }
        }
        if (nBkg > 0.0) { 
	        pmutex.lock();
            diffCRmap[i] = nEvents/nBkg-CRmap[i]; 
            bkgMap[i] = nBkg;
	        pmutex.unlock();
        } 
    }
}

void illh::mlh_normalization(
        unsigned int timeidxMin, 
        unsigned int timeidxMax, 
        unsigned int nTimesteps, 
        std::vector<DetectorPtr> detectors,
        SkyMap &diffCRmap,
        SkyMap &CRmap,
		std::mutex &pmutex)
{ 
    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
    {
        std::vector<DetectorPtr>::iterator det_it;
        for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
        { 
            DetectorPtr det = *det_it;
            double nEvents = 0.0; 
            double nBkg = 0.0; 
            unsigned int npix = CRmap.Npix();
        
            // Integrate over all (rotated) pixels 
            for (unsigned int i=0; i<npix; i++) 
            { 
                int j = illh::loc2eq_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
                if (det->FOV[i] && ((*det->Emap0)[i] > 0.0)) { 
                        nEvents += (*det->Nmap[timeidx])[i];
                        nBkg += (*det->Emap)[i]*(CRmap[j]+diffCRmap[j]);
                } 
            } 
            if (nBkg > 0.0) { 
	            pmutex.lock();
                det->norm[timeidx] = nEvents/nBkg; 
	            pmutex.unlock();
            } 
        } 
    }
}

void illh::mlh_significance(
        unsigned int pixmin, 
        unsigned int pixmax, 
        unsigned int nTimesteps,
        const std::vector<DetectorPtr> detectors,
        SkyMap &diffCRmap,
        SkyMap &CRmap,
        SkyMap &errormap,
        SkyMap &expectationMap,
        SkyMap &muon,
        SkyMap &muoff,
        SkyMap &Na,
        double &llhtemp,
		std::mutex &pmutex)
{
    for (unsigned int i=pixmin; i<pixmax;i++) 
    { 
        for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
        { 
            double mu(0.);
            double mu0(0.);
            double Njt(0.);
            double error = 0;
            std::vector<DetectorPtr>::const_iterator det_it;
            for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
            { 
                DetectorPtr det = *det_it;
                //rotation from Equatorial (ra,dec) to local (theta,phi)
                int j = illh::eq2loc_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);

                // global significance 
                double ej = (*det->Emap)[j];
                double ej0 = (*det->Emap0)[j];
                double nt = det->norm[timeidx];
                double nt0 = det->norm0[timeidx];
                double Njt_temp = (*det->Nmap[timeidx])[j];

                if (ej0*nt0 > 0.0 && ej*nt > 0.0) { 
                    mu += nt*ej;
                    mu0 += nt0*ej0;
                    Njt += Njt_temp;
                }
            }
            if ((mu>0) && (mu0>0))
            {
	            pmutex.lock();
                mu *= (diffCRmap[i]+CRmap[i]);
                mu0 *= CRmap[i];
                muon[i] += mu;
                muoff[i] += mu0;
                Na[i] += Njt;
                llhtemp += Njt*(mu*log(mu)-mu0*log(mu0));
                errormap[i]+=pow(mu-Njt,2);
	            pmutex.unlock();
            }
        }
        if (muon[i]> 0)
        {
	      pmutex.lock();
          errormap[i] = sqrt(errormap[i])/muon[i];
          expectationMap[i]=muon[i];
	      pmutex.unlock();
        }
    }

}


void illh::mlh_smooth(
        unsigned int pixmin, 
        unsigned int pixmax, 
        unsigned int nTimesteps,
        double smoothing_radius,
        SkyMap &significancemap,
        SkyMap &muon,
        SkyMap &muoff,
        SkyMap &Na,
		std::mutex &pmutex)
{
    for (unsigned int i=pixmin; i<pixmax;i++) 
    { 
        std::vector<int> listpix;
        pointing dir = significancemap.pix2ang(i);
        significancemap.query_disc(dir, smoothing_radius, listpix);

        double muontmp(0);
        double muofftmp(0);
        double Natmp(0);
        for(const int& k : listpix) 
        {

            if (muoff[k]) {
                muontmp += muon[k];
                muofftmp += muoff[k];
                Natmp += Na[k];
            }
        } 
        if ((muontmp <= 0) || (muofftmp <= 0)) {
            continue;
        }
        double smoothed_ri = muontmp/muofftmp - 1;
        double vtemp = -muontmp + muofftmp + Natmp*log(muontmp/muofftmp);
        pmutex.lock();
        significancemap[i] = sqrt(2.0*abs(vtemp));
        significancemap[i] *= (smoothed_ri < 0 ? -1 : 1);
        pmutex.unlock();
    }
}


void illh::mlh_acceptance(
        unsigned int pixmin, 
        unsigned int pixmax, 
        unsigned int nTimesteps,
        std::vector<DetectorPtr> detectors,
        const SkyMap &diffCRmap,
        const SkyMap &CRmap,
		std::mutex &pmutex)
{
    for (unsigned int i=pixmin; i<pixmax;i++) 
    {
        std::vector<DetectorPtr>::iterator det_it;
        for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
        {
            DetectorPtr det = *det_it;
            double nEvents = 0.0; 
            double nBkg = 0.0; 

            if (!((*det->Emap0)[i] > 0.0)) 
                continue;
    

            for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
            { 
                int j = illh::loc2eq_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
                nEvents += (*det->Nmap[timeidx])[i];
                nBkg += det->norm[timeidx]*(CRmap[j]+diffCRmap[j]); 
            }
            if (((*det->Emap0)[i] > 0.0) && (nBkg > 0.0)) 
            { 
                (*det->Emap)[i] = nEvents/nBkg; 
            } else { 
                (*det->Emap)[i] = 0.0; 
            }
        }
    }
}
