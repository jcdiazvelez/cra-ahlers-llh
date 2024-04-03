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
    bool save_iter_flag;
    bool iso;
    bool done(false);
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
             ("save-iter", po::bool_switch(&save_iter_flag)->default_value(false), "save each iteration")
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
                    j = illh::eq2loc_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
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
    /**
    fitshandle fitsOut; 
    stringstream datamapname;
    datamapname << foldername << boost::format("/data_%s_%d_%d.fits.gz") % detector_names_str.str() % nsideOut % nTimesteps;
    if (fs::exists(datamapname.str()) ) {
            fs::remove( datamapname.str() );
    }
    fitsOut.create(datamapname.str().c_str());
    write_Healpix_map_to_fits(fitsOut, dataMap, MyDTYPE);
    fitsOut.close();
    */

    //*****************************************************************************
    ////// Iterate //////////////////////////////////////////////////////////////// 
    //*****************************************************************************

    // n^th differential CR anisotropy : I_a^(n) - isovalue
    SkyMap diffCRmap;
    diffCRmap.SetNside(nsideOut, RING);
    diffCRmap.fill(0.);

    std::vector<double> llh;

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
                            int j = illh::eq2loc_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
                            if (det->FOV[j] && ((*det->Emap)[j] > 0.0)) { 
                                nEvents += (*det->Nmap[timeidx])[j]; 
                                nBkg += det->norm[timeidx]*(*det->Emap)[j]; 
                            } 
                    }
                }
                if (nBkg > 0.0) { 
                    diffCRmap[i] = nEvents/nBkg-CRmap[i]; 
                    bkgMap[i] = nBkg;
                } 
            } 

            // remove m=0 multipole moments :
            //const int LMAX=180;
            const int LMAX=2*nsideOut;

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
                        int j = illh::loc2eq_idx(i, timeidx, det->latitude, det->longitude, nTimesteps, CRmap);
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
            
                   
            // calculate statistical significance : S_a^(n)
            SkyMap variancemap; 
            variancemap.SetNside(nsideOut, RING);
            variancemap.fill(0.);

            // calculate statistical error: 
            SkyMap errormap; 
            errormap.SetNside(nsideOut, RING);
            errormap.fill(0.);

            // calculate expectation mu
            SkyMap expectationMap; 
            expectationMap.SetNside(nsideOut, RING);
            expectationMap.fill(0.);

            double llhtemp = 0;
            
            log_info("Significance " << iteration << " ...");
            for (unsigned int i=0; i < npix;i++ ) 
            { 
                    double musum(0.);
                    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                    { 
                        double mu(0.);
                        double mu0(0.);
                        double Njt(0.);
                        double error = 0;
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
                                variancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*ej*nt; 
                                variancemap[i] += 2.0*CRmap[i]*ej0*nt0; 
                                double temp1 = ej/ej0*nt/nt0;
                                variancemap[i] += 2.0*Njt_temp*log(temp1*(1.0+diffCRmap[i]/CRmap[i])); 

                                mu += nt*ej;
                                mu0 += nt0*ej0;
                                Njt += Njt_temp;
                                musum += mu;
                            }
                        }
                        if (mu>0 && mu0>0)
                        {
                            mu *= (diffCRmap[i]+CRmap[i]);
                            mu0 *= CRmap[i];
                            llhtemp+= Njt*(mu*log(mu)-mu0*log(mu0));
                            errormap[i]+=pow(mu-Njt,2);
                        }
                    }
                    if (musum> 0)
                    {
                      errormap[i] = sqrt(errormap[i])/musum;
                      expectationMap[i]=musum;
                    }
            }

            // cleaning out nonsensical values
            for (unsigned int i=0; i < npix;i++ ) 
            {
                if (variancemap[i] < 0)
                {
                    variancemap[i] = 0.;
                }
            }

            log_info("Finished iteration " << iteration << " of " << nIterations << "...");
            log_info("llh : " << llhtemp);
            double llhprev(0);

            if (llh.size()>0) 
            {
                log_info("llh ratio : " << 2*(llhtemp-llh[0]));
                llhprev = llh.back();
            }
            if (llh.size()>1) 
                log_info("llh ratio n, n-1: " << 2*(llhtemp-llh.back()));

            llh.push_back(llhtemp);

            if ((llh.size()>2) && ((llhtemp-llhprev) < 1e-4))
            {
                log_info("Converged: llh ratio n, n-1: " << 2*(llhtemp-llhprev));
                done = true;
            }

            if (iteration >= nIterations)
            {
                log_info("Reached max iterations");
                done = true;
            }


            if (save_iter_flag || done)
            { 
                // write N_tau^(n) normalizaton to file
                for (det_it = detectors.begin(); det_it != detectors.end(); det_it++) 
                { 
                    DetectorPtr det = *det_it; 
                    illh::save_iter(foldername, det->norm, *det->Emap, det->name, nsideOut, nTimesteps, iteration);
                } 

                // write relative intensity map I_a^(n)
                stringstream namefits;
                namefits << foldername << boost::format("/CR_%s_%d_%d_iteration%02d.fits.gz") % detector_names_str.str() % nsideOut % nTimesteps % iteration;

                if (fs::exists(namefits.str()) ) { 
                        fs::remove( namefits.str() ); 
                }
                fitshandle fitsOut; 
                fitsOut.create(namefits.str().c_str()); 
                fitsOut.add_comment("LLH reco Map");
                fitsOut.set_key("COORDS" , std::string("C"), "Equatorial coordinate system");
                fitsOut.set_key("TTYPE1", std::string("data"), "binned data counts");
                fitsOut.set_key("TTYPE2", std::string("background"), "llh reconstructed background");
                fitsOut.set_key("TTYPE3", std::string("relint"), "(alm-smoothed) relative intensity");
                write_Healpix_map_to_fits(fitsOut, dataMap, bkgMap, diffCRmapNormed, MyDTYPE);
                fitsOut.close(); 


                //write S_a^(n)
                stringstream nameSIGfits;
                nameSIGfits << foldername << boost::format("/variance_%s_%d_%d_iteration%02d.fits.gz") % detector_names_str.str() % nsideOut % nTimesteps % iteration;
                if (fs::exists(nameSIGfits.str()) ) {
                        fs::remove(nameSIGfits.str() );
                }
                fitsOut.create(nameSIGfits.str().c_str()); 
                fitsOut.add_comment(std::string("LLH stat Map"));
                fitsOut.set_key("COORDS" , std::string("C"), "Equatorial coordinate system");
                fitsOut.set_key("TTYPE1", std::string("sigma"), "statistical unsertainty");
                fitsOut.set_key("TTYPE2", std::string("lima"), "Li-Ma significance (squared)");
                fitsOut.set_key("TTYPE3", std::string("mu"), "statistical count expectation");
                write_Healpix_map_to_fits(fitsOut, errormap, variancemap, expectationMap, MyDTYPE);
                fitsOut.close(); 

         
                //  write llhratio^(i) 
                stringstream llhName; 
                llhName << foldername << boost::format("/llhratio_%s_%d_%d_iteration%02d.dat") % detector_names_str.str() % nsideOut % nTimesteps % iteration; 
                std::ofstream llhfileout(llhName.str().c_str()); 
                for(std::vector<double>::iterator it = llh.begin(); it != llh.end(); ++it) { 
                    llhfileout << *it << "\n"; 
                } 
                llhfileout.close();

            }
            if (done)
            {
                log_info("Finished");
                break;
            }


    }

    return 0;
}
