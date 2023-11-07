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
#ifdef ICECUBE
#include <healpix-cxx/healpix/alm.h>
#include <healpix-cxx/cxxsupport/xcomplex.h>
#include <healpix-cxx/healpix/alm_healpix_tools.h>
#include <healpix-cxx/cxxsupport/fitshandle.h>
#include <healpix-cxx/healpix/healpix_map.h>
#include <healpix-cxx/healpix/healpix_map_fitsio.h>
#define MyDTYPE FITSUTIL<double>::DTYPE 

#else
#include <healpix_cxx/alm.h>
#include <healpix_cxx/xcomplex.h>
#include <healpix_cxx/alm_healpix_tools.h>
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#define MyDTYPE PLANCK_FLOAT64
#endif //ICECUBE

#include <math.h>
#include <numeric>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/program_options.hpp> 
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include <boost/random.hpp> 
#include <boost/random/poisson_distribution.hpp> 
#include <boost/random/variate_generator.hpp> 
#include <boost/random/mersenne_twister.hpp> 
#include <boost/random/uniform_int.hpp> 
#include <boost/random/uniform_smallint.hpp> 


#include <memory>
#include <sstream>
#include <fstream>  
#include <iomanip>  
#include <string>
#include <stdexcept>
#include <iterator>
#include <vector>

#include <lhreco-proj/pickle.h>


#ifdef HAWCNEST
#include <hawcnest/HAWCNest.h>
#include <hawcnest/Logging.h>
#include <hawcnest/CommandLineConfigurator.h>
#include <detector-service/ConfigDirDetectorService.h>



#else
namespace po = boost::program_options; 
#define log_info(args...) cout << args << endl
#define log_fatal(args...) throw args
#endif


namespace fs = boost::filesystem; 
using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::random; 
using boost::format;

typedef Healpix_Map<double> SkyMap; 
typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer


/*
 * Sum over all pixels
 *
 */
double 
Sum(const SkyMap& map,unsigned int npix) 
{
    double sumval = 0.;
    for (unsigned int i=0; i < npix; i++) 
    {
         sumval += map[i];
    }
    return sumval;
}


#if __cplusplus > 199711L
void fluctuate(std::vector<SkyMapPtr>& Nmap, boost::mt19937 rng)
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
void isotropic(std::vector<SkyMapPtr>& Nmap, boost::mt19937 rng)
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
void loadMap(
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
void save_iter(
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
    expmapname << foldername << boost::format("/exposure_%s_%d_%d_iteration%02d.fits") % detector % nsideOut % nTimesteps % iteration;
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
loc2eq_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, SkyMap& CRmap)
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
eq2loc_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, SkyMap& CRmap)
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




int main(int argc, char* argv[])
{
    // declare all the variables
    std::string foldername;
    std::string detector1;
    std::string detector2;
    unsigned int nTimesteps;
    unsigned int timeidxMin;
    unsigned int timeidxMax;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideIn;
    int nsideOut;
    double lon1;
    double lat1;
    double lon2;
    double lat2;

    double thetamax1;
    double thetamax2;
    bool randfluct;
    bool iso;
    unsigned int seed;

    std::string prefix1;
    std::string prefix2;
    std::string suffix1;
    std::string suffix2;

    po::options_description desc("Options"); 
    po::variables_map vm; 
    
    try { 
       // Define and parse the program options 
       desc.add_options() 
             ("help,h", "produce help message") 
             ("prefix1", po::value<std::string>(&prefix1), "det1:prefix for input files") 
             ("prefix2", po::value<std::string>(&prefix2), "det2:prefix for input files") 
             ("suffix1", po::value<std::string>(&suffix1)->default_value(".fits.gz"), "det1:suffix for input files") 
             ("suffix2", po::value<std::string>(&suffix2)->default_value(".fits.gz"), "det2:suffix for input files") 
             ("outdir,o", po::value<std::string>(&foldername)->default_value("./sample/"), "Directory of output") 
             ("nsideout", po::value<int>(&nsideOut)->default_value(64), "Healpix Nside for output map") 
             ("timesteps", po::value<unsigned int>(&nTimesteps)->default_value(360), "Number of time steps") 
             ("timestepmin", po::value<unsigned int>(&timeidxMin)->default_value(0), "First time step to use") 
             ("timestepmax", po::value<unsigned int>(&timeidxMax)->default_value(0), "Last time step to use") 
             ("iterations", po::value<unsigned int>(&nIterations)->default_value(20), "Number of iterations") 
             ("sectors", po::value<unsigned int>(&nSectors)->default_value(1), "Number sectors") 
             ("lon1", po::value<double>(&lon1)->default_value(-97.0), "Longitude of detector1") 
             ("lat1", po::value<double>(&lat1)->default_value(19.0), "Latitude of detector1")
             ("lon2", po::value<double>(&lon2)->default_value(270.0), "Longitude of detector2") 
             ("lat2", po::value<double>(&lat2)->default_value(-90.0), "Latitude of detector2")
             ("thetamax1", po::value<double>(&thetamax1)->default_value(57.0), "max theta detector1")
             ("thetamax2", po::value<double>(&thetamax2)->default_value(65.0), "max theta detector2")
#if __cplusplus > 199711L
             ("fluctuate,f", po::bool_switch(&randfluct)->default_value(false), "add random fluctuations")
             ("seed", po::value<unsigned int>(&seed)->default_value(123), "RNG seed")
             ("iso", po::bool_switch(&iso)->default_value(false), "make isotropic map")
#endif
			 ("detector1", po::value<string>(&detector1)->default_value("HAWC"), "Name of detector1")
			 ("detector2", po::value<string>(&detector2)->default_value("IceCube"), "Name of detector2"); 
     
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
         
        /// --help option 
        if ( vm.count("help")  ) { 
            std::cout << "Basic Command Line Parameter App" 
                      << std::endl << desc << std::endl; 
            return 0; 
        } 
        po::notify(vm); // throws on error, so do after help in case there are any problems
        
        cout << prefix1 << " input files for " << detector1 << ":\n";
        cout << prefix2 << " input files for " << detector2 << ":\n";
     
        if (timeidxMax==0){
            timeidxMax=nTimesteps;
        }
        
    } catch(po::error& e) { 
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
            std::cerr << desc << std::endl; 
            return 1; 
        
    } catch(std::exception& e) { 
            std::cerr << "Unhandled Exception reached the top of main: " 
                      << e.what() << ", application will now exit" << std::endl; 
            return 2; 
    } 
 
    //*****************************************************************************
    ///// Initialize ///////////////////////////////////////////////////////////// 
    //*****************************************************************************
    unsigned int npix = 12*nsideOut*nsideOut; 

    // Import data : n_tau_i
    std::vector<SkyMapPtr> Nmap1;
    std::vector<SkyMapPtr> Nmap2;


    //*****************************************************************************
    ////// Read input maps //////////////////////////////////////////////////////// 
    //*****************************************************************************
    loadMap( Nmap1, timeidxMin, timeidxMax, nTimesteps, nsideIn, nsideOut, prefix1 , suffix1);
    loadMap( Nmap2, timeidxMin, timeidxMax, nTimesteps, nsideIn, nsideOut, prefix2 , suffix2);

    if (randfluct && iso)
            log_fatal("ranfluct and iso are mutually exclusive!!!");

    if (randfluct) 
    { 
#if __cplusplus > 199711L
            boost::mt19937 rng(seed);
            fluctuate(Nmap1,rng); 
            fluctuate(Nmap2,rng);
#else
            log_info("isotropic function disabled");
#endif
    }
    else if (iso) 
    { 
#if __cplusplus > 199711L
            boost::mt19937 rng(seed);
            isotropic(Nmap1,rng);
            isotropic(Nmap2,rng);
#else 
            log_info("isotropic function disabled");
#endif
    }

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }
    
    // detector position
    lon1 = lon1/180.*M_PI;
    lat1 = lat1/180.*M_PI;
    log_info("Sector 1 ("<<detector1<<"): Latitude=" << lat1 << ",  Longitude=" << lon1);

    lon2 = lon2/180.*M_PI;
    lat2 = lat2/180.*M_PI;
    log_info("Sector 2 ("<<detector2<<"): Latitude=" << lat2 << ",  Longitude=" << lon2);

    thetamax1 = thetamax1/180.*M_PI;
    thetamax2 = thetamax2/180.*M_PI;
    
    // normalization of isotropic flux
    double isovalue = 1.0;

    // initial CR anisotropy : I_a^(0) = 1.
    SkyMap CRmap;
    CRmap.SetNside(nsideOut, RING);
    CRmap.fill( isovalue );
                
    // initial exposure : A_i^(0)
    SkyMap Emap01;
    Emap01.SetNside(nsideOut, RING);
    Emap01.fill(0.);

    SkyMap Emap02;
    Emap02.SetNside(nsideOut, RING);
    Emap02.fill(0.);

    // n^th exposure : A_i^(n)
    SkyMap Emap1;
    Emap1.SetNside(nsideOut, RING);
    Emap1.fill(0.);

    SkyMap Emap2;
    Emap2.SetNside(nsideOut, RING);
    Emap2.fill(0.);

    //// window function of FOV map : F_a
    std::vector<bool> FOV1(npix);
    std::vector<bool> FOV2(npix);

    for (unsigned i=0;i<npix;++i) 
    {
		FOV1[i] = 1;
		FOV2[i] = 1;
        pointing pt = Emap2.pix2ang(i);
        if (pt.theta > thetamax1)
		   FOV1[i] = 0;
        if (pt.theta > thetamax2)
		   FOV2[i] = 0;
    }



    // Calculate initial exposure :	
    double totsector1 = 0.0;
    double totsector2 = 0.0;

    for (unsigned int i=0; i < npix;i++ ) 
    { 
            double nEvents1 = 0.0; 
            double nEvents2 = 0.0; 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    nEvents1    += (*Nmap1[timeidx])[i]; 
                    nEvents2    += (*Nmap2[timeidx])[i]; 
            } 
            if (FOV1[i]) {
                    Emap01[i] = nEvents1; 
                    totsector1 += nEvents1;
            }
            if (FOV2[i]) {
                    Emap02[i] = nEvents2; 
                    totsector2 += nEvents2;
            }
    }
            
    for (unsigned int i=0; i < npix;i++ ) 
    { 
            if (FOV1[i]) { 
                    Emap01[i] = Emap01[i]/totsector1; 
                    Emap1[i] = Emap01[i]; 
            }
            if (FOV2[i]) { 
                    Emap02[i] = Emap02[i]/totsector2; 
                    Emap2[i] = Emap02[i]; 
            }
    }
    log_info("Total events (det1): " << totsector1);
    log_info("Total events (det2): " << totsector2);
    log_info("Total events: " << totsector1+totsector2);

    // initial normalization : N_tau^(0)
    std::vector<double> norm01(nTimesteps);
    std::vector<double> norm02(nTimesteps);

    // n^th normalization : N_tau^(n)
    std::vector<double> norm1(nTimesteps);
    std::vector<double> norm2(nTimesteps);
    
    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
    { 
             double nEvents1 = 0.0;
             double nEvents2 = 0.0;
             for (unsigned int i=0; i < npix;i++ ) 
             { 
                     if (FOV1[i])
                        nEvents1 += (*Nmap1[timeidx])[i];
                     if (FOV2[i])
                        nEvents2 += (*Nmap2[timeidx])[i];
             } 
             norm1[timeidx] = nEvents1/isovalue; 
             norm01[timeidx] = norm1[timeidx]; 

             norm2[timeidx] = nEvents2/isovalue; 
             norm02[timeidx] = norm2[timeidx]; 
    }
        
    // data in equatorial coords : n_a
    SkyMap dataMap;
    dataMap.SetNside(nsideOut, RING); 
    dataMap.fill(0.);
    SkyMap bkgMap;
    bkgMap.SetNside(nsideOut, RING); 
    bkgMap.fill(0.);


    // rotate and combine maps
    for (unsigned int i=0; i < npix;i++ ) 
    { 
            double nBkg = 0.; 
            double nEvents = 0.; 
            
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    int j; 
                    j = loc2eq_idx(i, timeidx, lat1, lon1, nTimesteps, CRmap);
                    if (Emap01[j] > 0.0 ){ 
                            nEvents += (*Nmap1[timeidx])[j];
                            nBkg    += norm1[timeidx]*Emap01[j]; 
                    }
                    j = loc2eq_idx(i, timeidx, lat2, lon2, nTimesteps, CRmap);
                    if (Emap02[j] > 0.0 ){ 
                            nEvents += (*Nmap2[timeidx])[j];
                            nBkg    += norm2[timeidx]*Emap02[j]; 
                    }
                
            } 
            dataMap[i] = nEvents;
    } 
    
    //  write N_tau^(0), write A_i^(0)
    log_info("Writting initial exposure");
    save_iter(foldername, norm1, Emap01, detector1, nsideOut, nTimesteps, 0);
    save_iter(foldername, norm2, Emap02, detector2, nsideOut, nTimesteps, 0);


    // write n_a
    fitshandle fitsOut; 
    stringstream datamapname;
    datamapname <<  foldername << boost::format("/data_%s-%s_%d_%d.fits.gz") % detector1 % detector2 % nsideOut % nTimesteps;
    if (fs::exists(datamapname.str()) ) {
            fs::remove( datamapname.str() );
    }
    fitsOut.create(datamapname.str().c_str());
    write_Healpix_map_to_fits(fitsOut, dataMap, MyDTYPE);
    fitsOut.close();

    //*****************************************************************************
    ////// Iterate //////////////////////////////////////////////////////////////// 
    //*****************************************************************************

    // n^th differential CR anisotropy : I_a^(n) - isovalue
    SkyMap diffCRmap;
    diffCRmap.SetNside(nsideOut, RING);
    diffCRmap.fill(0.);

    for (unsigned int iteration = 1; iteration <= nIterations; iteration++)
    { 
            log_info("Iter " << iteration);

            // calculate new CR anisotropy
            for (unsigned int i=0; i<npix;i++) 
            { 
                    double nEvents = 0.0;
                    double nBkg = 0.0;

                    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                    { 
                            int j; 
                            j = loc2eq_idx(i, timeidx, lat1, lon1, nTimesteps, CRmap);
                            if (FOV1[j] && (Emap01[j] > 0.0)) { 
                                    nEvents += (*Nmap1[timeidx])[j];
                                    nBkg += norm1[timeidx]*Emap1[j]; 
                            } 
                            j = loc2eq_idx(i, timeidx, lat2, lon2, nTimesteps, CRmap);
                            if (FOV2[j] && (Emap02[j] > 0.0)) { 
                                    nEvents += (*Nmap2[timeidx])[j];
                                    nBkg += norm2[timeidx]*Emap2[j]; 
                            } 
                    }
                        
                    if (nBkg > 0.0) { 
                            diffCRmap[i] = nEvents/nBkg-CRmap[i]; 
                            bkgMap[i] = nBkg;
                    }
            } 

            // remove m=0 multipole moments :
            const int LMAX=180;

            // Initialize the spherical harmonic coefficients of the map
            Alm<xcomplex<double> > alm( LMAX, LMAX);
            arr<double> weight(2*nsideOut);
            weight.fill(1.);
            
	        // generate alm coefficients with terative map2alm method
            const int numIter = 3;
            map2alm_iter( diffCRmap, alm, numIter, weight);
            for (unsigned int l=0; l < LMAX+1; l++)
            {    
                alm(l,0) = 0.;
            }

            // generate diffCRmap from alm coefficients
            alm2map( alm, diffCRmap);

            // calculate new normalization 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                double nEvents1 = 0.0;
                double nBkg1 = 0.0; 
                double nEvents2 = 0.0;
                double nBkg2 = 0.0; 

                // Integrate over all (rotated) pixels 
                for (unsigned int i=0; i<npix; i++) 
                { 
                        int j;
                        j = eq2loc_idx(i, timeidx, lat1, lon1, nTimesteps, CRmap);
                        if (FOV1[i] && (Emap01[i] > 0.0)) { 
                                nEvents1 += (*Nmap1[timeidx])[i];
                                nBkg1 += Emap1[i]*(CRmap[j]+diffCRmap[j]);
                        } 

                        j = eq2loc_idx(i, timeidx, lat2, lon2, nTimesteps, CRmap);
                        if (FOV2[i] && (Emap02[i] > 0.0)) { 
                                nEvents2 += (*Nmap2[timeidx])[i];
                                nBkg2 += Emap2[i]*(CRmap[j]+diffCRmap[j]);
                        } 
                } 
                if (nBkg1 > 0.0) { 
                        norm1[timeidx] = nEvents1/nBkg1; 
                } 
                if (nBkg2 > 0.0) { 
                        norm2[timeidx] = nEvents2/nBkg2; 
                } 

            } 

            // caluculate new acceptance
            for (unsigned int i=0; i<npix; i++) 
            {
                double nEvents1 = 0.0;
                double nBkg1 = 0.0; 
                double nEvents2 = 0.0;
                double nBkg2 = 0.0; 
            
                if ( !(Emap01[i] > 0.0) && !( Emap02[i] > 0.0) ) 
                        continue;

                for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                {
                    int j;
                    j = eq2loc_idx(i, timeidx, lat1, lon1, nTimesteps, CRmap);
                    nEvents1 += (*Nmap1[timeidx])[i];
                    nBkg1 += norm1[timeidx]*(CRmap[j]+diffCRmap[j]);

                    j = eq2loc_idx(i, timeidx, lat2, lon2, nTimesteps, CRmap);
                    nEvents2 += (*Nmap2[timeidx])[i];
                    nBkg2 += norm2[timeidx]*(CRmap[j]+diffCRmap[j]);
                }
            
                if ((Emap01[i] > 0.0) && (nBkg1 > 0.0)) 
                {
                    Emap1[i] = nEvents1/nBkg1;
                } else { 
                    Emap1[i] = 0.0;
                }
                if ((Emap02[i] > 0.0) && (nBkg2 > 0.0)) 
                {
                    Emap2[i] = nEvents2/nBkg2;
                } else { 
                    Emap2[i] = 0.0;
                }

            }
            
            // Renormalize
            double sumEmap1 = Sum(Emap1, npix);
            double sumEmap2 = Sum(Emap2, npix);
            
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                norm1[timeidx] = norm1[timeidx]*sumEmap1; 
                norm2[timeidx] = norm2[timeidx]*sumEmap2; 
            } 
            for (unsigned int i=0; i < npix;i++ ) 
            { 
                    Emap1[i] = Emap1[i]/sumEmap1;
                    Emap2[i] = Emap2[i]/sumEmap2;
            } 
            
            SkyMap diffCRmapNormed; 
            diffCRmapNormed.SetNside(nsideOut, RING); 
            diffCRmapNormed.fill(0.);

            for (unsigned int i=0; i < npix;i++ ) 
            { 
                    diffCRmapNormed[i] = diffCRmap[i]/isovalue; 
            } 
            
            // write relative intensity map I_a^(n)
            stringstream namefits;
            namefits << foldername << boost::format("/CR_%s-%s_%d_%d_iteration%02d.fits.gz") % detector1 % detector2 % nsideOut % nTimesteps % iteration;
            if (fs::exists(namefits.str()) ) { 
                    fs::remove( namefits.str() ); 
            }
            fitshandle fitsOut; 
            fitsOut.create(namefits.str().c_str()); 
            write_Healpix_map_to_fits(fitsOut, diffCRmapNormed, MyDTYPE);
            fitsOut.close(); 

       
            // write N_tau^(n) normalizaton to file
            save_iter(foldername, norm1, Emap1, detector1, nsideOut, nTimesteps, iteration);
            save_iter(foldername, norm2, Emap2, detector2, nsideOut, nTimesteps, iteration);

            // calculate statistical significance : S_a^(n)
            SkyMap significancemap; 
            significancemap.SetNside(nsideOut, RING);
            significancemap.fill(0.);
            
            log_info("Significance " << iteration << " ...");
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    for (unsigned int i=0; i < npix;i++ ) 
                    { 
                            //rotation from local to Equatorial (ra,dec) 
                            int j;
                            j = loc2eq_idx(i, timeidx, lat1, lon1, nTimesteps, CRmap);

                            // global significance 
                            if (Emap01[j]*norm01[timeidx] > 0.0 && Emap1[j]*norm1[timeidx] > 0.0) { 
                                    significancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*Emap1[j]*norm1[timeidx]; 
                                    significancemap[i] += +2.0*CRmap[i]*Emap01[j]*norm01[timeidx]; 
                                    double temp1 = Emap1[j]/Emap01[j]*norm1[timeidx]/norm01[timeidx];
                                    significancemap[i] += 2.0*(*Nmap1[timeidx])[j]*log(temp1*(1.0+diffCRmap[i]/CRmap[i])); 
                            }

                            j = loc2eq_idx(i, timeidx, lat2, lon2, nTimesteps, CRmap);
                            if (Emap02[j]*norm02[timeidx] > 0.0 && Emap2[j]*norm2[timeidx] > 0.0) { 
                                    significancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*Emap2[j]*norm2[timeidx]; 
                                    significancemap[i] += +2.0*CRmap[i]*Emap02[j]*norm02[timeidx]; 
                                    double temp2 = Emap2[j]/Emap02[j]*norm2[timeidx]/norm02[timeidx];
                                    significancemap[i] += 2.0*(*Nmap2[timeidx])[j]*log(temp2*(1.0+diffCRmap[i]/CRmap[i])); 
                            }

                    } 
            }

            //write S_a^(n)
            stringstream nameSIGfits;
            nameSIGfits << foldername <<  boost::format("/significance_%s-%s_%d_%d_iteration%02d.fits.gz") % detector1 % detector2 % nsideOut % nTimesteps % iteration; 
            if (fs::exists(nameSIGfits.str()) ) {
                    fs::remove(nameSIGfits.str() );
            }
            fitsOut.create(nameSIGfits.str().c_str()); 
            write_Healpix_map_to_fits(fitsOut, significancemap, MyDTYPE);
            fitsOut.close(); 

            log_info("Finished iteration " << iteration << " of " << nIterations << "...");

    }

    fitshandle bkgfitsOut; 
    stringstream bkgmapname;
    bkgmapname <<  foldername << boost::format("/background_%s-%s_%d_%d.fits.gz") % detector1 % detector2 % nsideOut % nTimesteps;
    if (fs::exists(bkgmapname.str()) ) {
            fs::remove( bkgmapname.str() );
    }
    bkgfitsOut.create(bkgmapname.str().c_str());
    write_Healpix_map_to_fits(bkgfitsOut, bkgMap, MyDTYPE);
    bkgfitsOut.close();

    return 0;
}
