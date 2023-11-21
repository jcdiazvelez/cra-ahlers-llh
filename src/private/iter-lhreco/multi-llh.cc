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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


#include <memory>
#include <sstream>
#include <fstream>  
#include <iomanip>  
#include <string>
#include <stdexcept>
#include <iterator>
#include <vector>

#include <iter-lhreco-proj/pickle.h>


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
namespace pt = boost::property_tree;

using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::random; 
using boost::format;

typedef Healpix_Map<double> SkyMap; 
typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer



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



class Detector {
  public:         
    std::string name; 
    unsigned int nside; 
    double longitude; 
    double latitude; 
    double thetamax; 

    double totsector;

    std::string prefix;
    std::string suffix;

    std::vector<SkyMapPtr> Nmap;

    // initial exposure : A_i^(0)
    SkyMapPtr Emap0;

    // n^th exposure : A_i^(n)
    SkyMapPtr Emap;


    //// window function of FOV map : F_a
    std::vector<bool> FOV;

    // initial normalization : N_tau^(0)
    std::vector<double> norm0;

    // n^th normalization : N_tau^(n)
    std::vector<double> norm;
};




typedef boost::shared_ptr<Detector> DetectorPtr; // Detector shared pointer


int main(int argc, char* argv[])
{
    // declare all the variables
    // Create a root
    pt::ptree root;
    std::vector<DetectorPtr> detectors;
    // Load the json file in this ptree

    std::string foldername;
    std::string config_path;
    unsigned int nTimesteps;
    unsigned int timeidxMin;
    unsigned int timeidxMax;
    unsigned int nIterations;
    unsigned int nSectors;
    unsigned int nsideOut;
    unsigned int npix;
    bool randfluct;
    bool iso;
    unsigned int seed;

    po::options_description desc("Options"); 
    po::variables_map vm; 
    
    try { 
       // Define and parse the program options 
       desc.add_options() 
             ("help,h", "produce help message") 
             ("outdir,o", po::value<std::string>(&foldername)->default_value("./sample/"), "Directory of output") 
             ("nsideout", po::value<unsigned int>(&nsideOut)->default_value(64), "Healpix Nside for output map") 
             ("timesteps", po::value<unsigned int>(&nTimesteps)->default_value(360), "Number of time steps") 
             ("timestepmin", po::value<unsigned int>(&timeidxMin)->default_value(0), "First time step to use") 
             ("timestepmax", po::value<unsigned int>(&timeidxMax)->default_value(0), "Last time step to use") 
             ("iterations", po::value<unsigned int>(&nIterations)->default_value(20), "Number of iterations") 
#if __cplusplus > 199711L
             ("fluctuate,f", po::bool_switch(&randfluct)->default_value(false), "add random fluctuations")
             ("seed", po::value<unsigned int>(&seed)->default_value(123), "RNG seed")
             ("iso", po::bool_switch(&iso)->default_value(false), "make isotropic map")
#endif
			 ("config", po::value<string>(&config_path)->default_value("config.json"), "JSON config");
     
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
         
        /// --help option 
        if ( vm.count("help")  ) { 
            std::cout << "Basic Command Line Parameter App" 
                      << std::endl << desc << std::endl; 
            return 0; 
        } 
        po::notify(vm); // throws on error, so do after help in case there are any problems
        
        pt::read_json(config_path.c_str(), root);
        
        // Iterator over all detectors

        npix = 12*nsideOut*nsideOut; 

        for (pt::ptree::value_type &detector: root.get_child("detectors"))
        {
            // Detector is a std::pair of a string and a child
            // Get the label of the node 
            // Get the content of the node
            DetectorPtr det = DetectorPtr(new Detector); 
            det->name = detector.first;
            det->nside = detector.second.get<unsigned int>("nside");

            // detector position
            det->longitude = detector.second.get<double>("longitude")*M_PI/180.;
            det->latitude = detector.second.get<double>("latitude")*M_PI/180.;

            log_info("Sector  ("<<det->name<<"): Latitude=" << det->latitude*180./M_PI<< 
                    ",  Longitude=" << det->longitude*190./M_PI);

            det->thetamax = detector.second.get<double>("thetamax")*M_PI/180.;

            det->prefix = detector.second.get<std::string>("prefix");
            det->suffix = detector.second.get<std::string>("suffix");

            // initial exposure : A_i^(0)
            det->Emap0 = SkyMapPtr(new SkyMap);
            det->Emap0->SetNside(nsideOut, RING);
            det->Emap0->fill(0.);


            // n^th exposure : A_i^(n)
            det->Emap = SkyMapPtr(new SkyMap);
            det->Emap->SetNside(nsideOut, RING);
            det->Emap->fill(0.);

            det->FOV.resize(npix, false);
            det->norm0.resize(nTimesteps);
            det->norm.resize(nTimesteps);

            detectors.push_back(det);
        }

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

    // Import data : n_tau_i


    stringstream detector_names_str;
    std::vector<DetectorPtr>::iterator det_it;
    double totsectors(0);
    for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
    {
        DetectorPtr det = *det_it;
        cout << det->prefix << " input files for " << det->name << ":\n";
     
        detector_names_str  << boost::format("%s_") % det->name;

        //*****************************************************************************
        ////// Read input maps //////////////////////////////////////////////////////// 
        //*****************************************************************************
        illh::loadMap( det->Nmap, timeidxMin, timeidxMax, nTimesteps, det->nside, nsideOut, det->prefix , det->suffix);

        if (randfluct && iso)
            log_fatal("ranfluct and iso are mutually exclusive!!!");

        if (randfluct) 
        { 
#if __cplusplus > 199711L
            boost::mt19937 rng(seed);
            fluctuate(det->Nmap,rng); 
#else
            log_info("isotropic function disabled");
#endif
        }
        else if (iso) 
        { 
#if __cplusplus > 199711L
            boost::mt19937 rng(seed);
            isotropic(det->Nmap,rng);
#else 
            log_info("isotropic function disabled");
#endif
        }

        for (unsigned i=0;i<npix;++i) 
        { 
               det->FOV[i] = 1; 
               pointing pt = det->Emap->pix2ang(i); 
               if (pt.theta > det->thetamax) 
                   det->FOV[i] = 0; 
        } 
        // Calculate initial exposure :	
        det->totsector = 0; 
        for (unsigned int i=0; i < npix;i++ ) 
        { 
            double nEvents = 0.0; 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                nEvents    += (*det->Nmap[timeidx])[i]; 
            } 
            if (det->FOV[i]) { 
                (*det->Emap0)[i] = nEvents; 
                det->totsector += nEvents; 
            } 
        } 
        for (unsigned int i=0; i < npix;i++ ) 
        { 
            if (det->FOV[i]) { 
                (*det->Emap0)[i] = (*det->Emap0)[i]/det->totsector; 
                (*det->Emap)[i] = (*det->Emap0)[i]; 
            } 
        } 
        log_info("Total events ("<<det->name<<"): " << det->totsector); 
        totsectors += det->totsector;
    }
    log_info("Total events: " << totsectors);

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }
    

    // normalization of isotropic flux
    double isovalue = 1.0;

    // initial CR anisotropy : I_a^(0) = 1.
    SkyMap CRmap;
    CRmap.SetNside(nsideOut, RING);
    CRmap.fill( isovalue );

    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
    { 
        for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
        {
             DetectorPtr det = *det_it;
             double nEvents = 0.0;
             for (unsigned int i=0; i < npix;i++ ) 
             { 
                     if (det->FOV[i])
                        nEvents += (*det->Nmap[timeidx])[i];
             } 
             det->norm[timeidx] = nEvents/isovalue; 
             det->norm0[timeidx] = det->norm[timeidx]; 
        }
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
                for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
                {
                    DetectorPtr det = *det_it;
                    int j; 
                    j = illh::loc2eq_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
                    if ((*det->Emap0)[j] > 0.0 ){ 
                        nEvents += (*det->Nmap[timeidx])[j];
                        nBkg    += det->norm[timeidx]*(*det->Emap0)[j]; 
                    }
                }
            } 
            dataMap[i] = nEvents;
    } 
    
    //  write N_tau^(0), write A_i^(0)
    log_info("Writting initial exposure");
    for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
    {
        DetectorPtr det = *det_it;
        illh::save_iter(foldername, det->norm, *det->Emap0, det->name, nsideOut, nTimesteps, 0);
    }


    // write n_a
    fitshandle fitsOut; 
    stringstream datamapname;
    datamapname << foldername << boost::format("/data_%s_%d_%d.fits.gz") % detector_names_str.str() % nsideOut % nTimesteps;
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
                    for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
                    {
                            DetectorPtr det = *det_it;
                            int j = illh::loc2eq_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
                            if (det->FOV[j] && ((*det->Emap0)[j] > 0.0)) { 
                                nEvents += (*det->Nmap[timeidx])[j]; 
                                nBkg += det->norm[timeidx]*(*det->Emap)[j]; 
                            } 
                    }
                }
                if (nBkg > 0.0) { 
                    diffCRmap[i] = nEvents/nBkg-CRmap[i]; 
                    bkgMap[i] = nBkg;
                } } 

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
                for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
                { 
                    DetectorPtr det = *det_it;
                    double nEvents = 0.0; 
                    double nBkg = 0.0; 
                
                    // Integrate over all (rotated) pixels 
                    for (unsigned int i=0; i<npix; i++) 
                    { 
                        int j = illh::eq2loc_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
                        if (det->FOV[i] && ((*det->Emap0)[i] > 0.0)) { 
                                nEvents += (*det->Nmap[timeidx])[i];
                                nBkg += (*det->Emap)[i]*(CRmap[j]+diffCRmap[j]);
                        } 
                    } 
                    if (nBkg > 0.0) { 
                        det->norm[timeidx] = nEvents/nBkg; 
                    } 
                } 
            }

            // caluculate new acceptance
            for (unsigned int i=0; i<npix; i++) 
            {
                for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
                {
                    DetectorPtr det = *det_it;
                    double nEvents = 0.0; 
                    double nBkg = 0.0; 

                    if (!((*det->Emap0)[i] > 0.0)) 
                        continue;
            

                    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                    { 
                        int j = illh::eq2loc_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
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
            
            // Renormalize
            for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
            { 
                DetectorPtr det = *det_it;
                double sumEmap = illh::Sum(*det->Emap, npix); 
                for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                { 
                    det->norm[timeidx] = det->norm[timeidx]*sumEmap; 
                } 
                for (unsigned int i=0; i < npix;i++ ) 
                { 
                    (*det->Emap)[i] = (*det->Emap)[i]/sumEmap;
                } 
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
            namefits << foldername << boost::format("/CR_%s_%d_%d_iteration%02d.fits.gz") % detector_names_str.str() % nsideOut % nTimesteps % iteration;

            if (fs::exists(namefits.str()) ) { 
                    fs::remove( namefits.str() ); 
            }
            fitshandle fitsOut; 
            fitsOut.create(namefits.str().c_str()); 
            write_Healpix_map_to_fits(fitsOut, diffCRmapNormed, dataMap, bkgMap, MyDTYPE);
            fitsOut.close(); 

       
            // write N_tau^(n) normalizaton to file
            for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
            { 
                DetectorPtr det = *det_it;
                illh::save_iter(foldername, det->norm, *det->Emap, det->name, nsideOut, nTimesteps, iteration);
            }

            // calculate statistical significance : S_a^(n)
            SkyMap significancemap; 
            SkyMap sigmu_on; 
            SkyMap sigmu_off; 
            SkyMap signtot; 
            significancemap.SetNside(nsideOut, RING);
            significancemap.fill(0.);
            sigmu_on.SetNside(nsideOut, RING);
            sigmu_on.fill(0.);
            sigmu_off.SetNside(nsideOut, RING);
            sigmu_off.fill(0.);
            signtot.SetNside(nsideOut, RING);
            signtot.fill(0.);
            
            log_info("Significance " << iteration << " ...");
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    for (unsigned int i=0; i < npix;i++ ) 
                    { 
                        for (det_it = detectors.begin(); det_it != detectors.end(); det_it++)
                        { 
                            DetectorPtr det = *det_it;
                            //rotation from local to Equatorial (ra,dec) 
                            int j = illh::loc2eq_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);

                            // global significance 
                            if ((*det->Emap0)[j]*det->norm0[timeidx] > 0.0 && (*det->Emap)[j]*det->norm[timeidx] > 0.0) { 
                                significancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*(*det->Emap)[j]*det->norm[timeidx]; 
                                double temp1 = (*det->Emap)[j]/(*det->Emap0)[j]*det->norm[timeidx]/det->norm0[timeidx];
                                significancemap[i] += 2.0*(*det->Nmap[timeidx])[j]*log(temp1*(1.0+diffCRmap[i]/CRmap[i])); 
                            }
                        }
                    } 
            }



            //write S_a^(n)
            stringstream nameSIGfits;
            nameSIGfits << foldername << boost::format("/significance_%s_%d_%d_iteration%02d.fits.gz") % detector_names_str.str() % nsideOut % nTimesteps % iteration;
            if (fs::exists(nameSIGfits.str()) ) {
                    fs::remove(nameSIGfits.str() );
            }
            fitsOut.create(nameSIGfits.str().c_str()); 
            write_Healpix_map_to_fits(fitsOut, significancemap, MyDTYPE);
            fitsOut.close(); 

            log_info("Finished iteration " << iteration << " of " << nIterations << "...");

    }

    return 0;
}
