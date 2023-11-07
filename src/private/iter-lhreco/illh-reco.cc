/**
 * illh-reco
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

typedef Healpix_Map<double> SkyMap; 
typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer




int main(int argc, char* argv[])
{
//*****************************************************************************
////// Input ////////////////////////////////////////////////////////////////// 
//*****************************************************************************
    std::string foldername;
    std::string detector;
    std::string prefix;
    std::string suffix;
    unsigned int nTimesteps;
    unsigned int timeidxMin;
    unsigned int timeidxMax;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideIn;
    int nsideOut;
    double lon;
    double lat;
    double thetamax;
    std::vector< std::string > input;

#ifdef HAWCNEST
    CommandLineConfigurator cl;
  
    //cl.AddPositionalOption<std::vector <std::string> >("input", "HEALPix local maps (theta/phi) in timesteps of a sidereal day");
    cl.AddOption<std::string>("outdir","./sample/","Directory of output");
    cl.AddOption<std::string>("detector","HAWC","name of detector");
    cl.AddOption<int>("nsideout",64,"nSide of output maps");
    cl.AddOption<int>("iterations",20,"Number of iterations");
    cl.AddOption<int>("timesteps",360,"Timesteps in sidereal day");
    cl.AddOption<int>("timestepmin",0,"first time step to use");
    cl.AddOption<int>("timestepmax",0,"last time step to use (default = timesteps)");
    cl.AddOption<int>("sectors",1,"Number of FOV sectors to use");
    cl.AddOption<double>("lon",-97.0,"Longitude of detector");
    cl.AddOption<double>("lat",19.0,"Latitude of detector");
    if (!cl.ParseCommandLine(argc, argv)) {
      return 1;
    }
    

    //input = cl.GetArgument<std::vector <std::string> >("input");
    foldername = cl.GetArgument<std::string>("outdir");
    detector = cl.GetArgument<std::string>("detector");
    nTimesteps = cl.GetArgument<int>("timesteps");
    timeidxMin = cl.GetArgument<int>("timestepmin");
    timeidxMax = cl.GetArgument<int>("timestepmax");
    if (timeidxMax == 0){
      timeidxMax=nTimesteps;
    }
    nIterations = cl.GetArgument<int>("iterations");
    nSectors = cl.GetArgument<int>("sectors");
    nsideOut = cl.GetArgument<int>("nsideout");
    lon = cl.GetArgument<double>("lon");
    lat = cl.GetArgument<double>("lat");
#else

    po::options_description desc("Options"); 
    po::positional_options_description p;
    //p.add("input", -1);
    po::variables_map vm; 
    
    try { 
       // Define and parse the program options 
       desc.add_options() 
             ("help,h", "produce help message") 
             ("prefix", po::value<std::string>(&prefix), "prefix for input files") 
             ("suffix", po::value<std::string>(&suffix)->default_value(".fits.gz"), "suffix for input files") 
             //("input", po::value<std::vector <std::string> >(&input)->multitoken(), "Input files") 
             ("outdir,o", po::value<std::string>(&foldername)->default_value("./sample/"), "Directory of output") 
             ("nsideout", po::value<int>(&nsideOut)->default_value(64), "Healpix Nside for output map") 
             ("timesteps", po::value<unsigned int>(&nTimesteps)->default_value(360), "Number of time steps") 
             ("timestepmin", po::value<unsigned int>(&timeidxMin)->default_value(0), "First time step to use") 
             ("timestepmax", po::value<unsigned int>(&timeidxMax)->default_value(0), "Last time step to use") 
             ("iterations", po::value<unsigned int>(&nIterations)->default_value(20), "Number of iterations") 
             ("sectors", po::value<unsigned int>(&nSectors)->default_value(1), "Number sectors") 
             ("lon", po::value<double>(&lon)->default_value(-97.0), "Longitude of detector") 
             ("lat", po::value<double>(&lat)->default_value(19.0), "Latitude of detector")
             ("thetamax", po::value<double>(&thetamax)->default_value(70.0), "max theta")
			("detector", 
			 po::value<string>(&detector)->default_value("HAWC"), 
			 "Name of detector"); 
     
        //po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw 
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm); 
         
        /// --help option 
        if ( vm.count("help")  ) { 
            std::cout << "Basic Command Line Parameter App" 
                      << std::endl << desc << std::endl; 
            return 0; 
        } 
        po::notify(vm); // throws on error, so do after help in case there are any problems
        //cout << prefix << " input files for " << detector << ":\n";
        
        //if (vm.count("input")) 
        //{ 
         //   cout << "Input files: " << "\n";
        //}
     
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
#endif 
 
//*****************************************************************************
////// Initialize ///////////////////////////////////////////////////////////// 
//*****************************************************************************
    unsigned int npix = 12*nsideOut*nsideOut; 

    // Import data : n_tau_i
    std::vector<SkyMapPtr> Nmap;

    //*****************************************************************************
    ////// Read input maps //////////////////////////////////////////////////////// 
    //*****************************************************************************
    illh::loadMap( Nmap, timeidxMin, timeidxMax, nTimesteps, nsideIn, nsideOut, prefix , suffix);

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }
    
    // detector position
    lon = lon/180.*M_PI;
    lat = lat/180.*M_PI;
    thetamax = thetamax/180.*M_PI;
#ifdef HAWCNEST
    if (detector == "HAWC"){
      string configDir = "";
      if (getenv("CONFIG_HAWC")){
        configDir = getenv("CONFIG_HAWC");
      } else {
        log_fatal("config-dir is not defined by CONFIG_HAWC");
      }
      HAWCNest nest;
      nest.Service<ConfigDirDetectorService>("det")
        ("configDir",configDir);
      nest.Configure();
      const det::Detector& detector = GetService<DetectorService>("det").GetDetector(TimeStamp(1111111111));

      const LatLonAlt loc = detector.GetLatitudeLongitudeHeight();
      lon = -loc.GetLongitude();
      lat = loc.GetLatitude();
    }
#endif
    double clat = cos(lat);
    double slat = sin(lat);
    log_info("Latitude=" << lat << ",  Longitude=" << lon);
    
    // normalization of isotropic flux
    double isovalue = 1.0;
    // initial CR anisotropy : I_a^(0) = 1.
    SkyMap CRmap;
    CRmap.SetNside(nsideOut, RING);
    CRmap.fill( isovalue );
                
    // initial exposure : A_i^(0)
    SkyMap Emap0;
    Emap0.SetNside(nsideOut, RING);
    Emap0.fill(0.);
    // n^th exposure : A_i^(n)
    SkyMap Emap;
    Emap.SetNside(nsideOut, RING);
    Emap.fill(0.);
    double nEventsTot = 0.0;

    //// window function of FOV map : F_a
    std::vector<bool> FOV(npix);
    for (unsigned i=0;i<npix;++i) 
    {
		FOV[i] = 1;
        pointing pt = Emap.pix2ang(i);
        if (pt.theta > thetamax)
		   FOV[i] = 0;
    }



    for (unsigned int i=0; i < npix;i++ ) 
    { 
            double nEvents = 0.0; 
            //for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    nEvents    += (*Nmap[timeidx])[i];
            } 
            if (FOV[i]) {
                Emap0[i] = nEvents; 
                nEventsTot += nEvents;
            }
    }
            
    for (unsigned int i=0; i < npix;i++ ) 
    { 
            if (FOV[i]) { 
                Emap0[i] = Emap0[i]/nEventsTot; 
                Emap[i] = Emap0[i]; 
            }
    }
    log_info("Total events " << nEventsTot);

    // initial normalization : N_tau^(0)
    std::vector<double> norm0(nTimesteps);
    // n^th normalization : N_tau^(n)
    std::vector<double> norm(nTimesteps);
    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
    { 
             double nEvents = 0.0;
             for (unsigned int i=0; i < npix;i++ ) 
             { 
                     if (FOV[i]) 
                         nEvents += (*Nmap[timeidx])[i];
             } 
             norm[timeidx] = nEvents/isovalue; 
             norm0[timeidx] = norm[timeidx]; 
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
        
                    int j;  // local pixel
                    j = illh::loc2eq_idx(i, timeidx, lat, lon, nTimesteps, CRmap);
                    if (Emap0[j] > 0.0 ){ 
                            nEvents += (*Nmap[timeidx])[j];
                            nBkg    += norm[timeidx]*Emap0[j]; 
                    } 
            } 
            dataMap[i]      = nEvents;
    } 

    
    //  write N_tau^(0), write A_i^(0)
    log_info("Writting initial exposure");
    illh::save_iter(foldername, norm, Emap0, detector, nsideOut, nTimesteps, 0);

    // write n_a
    fitshandle fitsOut; 
    stringstream datamapname;
    datamapname <<  foldername << boost::format("/data_%s_%d_%d.fits") % detector % nsideOut % nTimesteps;
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
                            j = illh::loc2eq_idx(i, timeidx, lat, lon, nTimesteps, CRmap);

                            if (FOV[j] && (Emap0[j] > 0.0)) { 
                                    nEvents += (*Nmap[timeidx])[j];
                                    nBkg += norm[timeidx]*Emap[j]; 
                            } 
                    }
                        
                    if (nBkg > 0.0) { 
                            diffCRmap[i] = nEvents/nBkg-CRmap[i]; 
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
                double nEvents = 0.0;
                double nBkg = 0.0; 

                // Integrate over all (rotated) pixels 
                for (unsigned int i=0; i<npix; i++) 
                { 

                        int j;
                        j = illh::eq2loc_idx(i, timeidx, lat, lon, nTimesteps, CRmap);
                        if (FOV[i] && (Emap0[i] > 0.0)) { 
                                nEvents += (*Nmap[timeidx])[i];
                                nBkg += Emap[i]*(CRmap[j]+diffCRmap[j]);
                        } 
                } 
                if (nBkg > 0.0) { 
                        norm[timeidx] = nEvents/nBkg; 
                } 
            } 

            // caluculate new acceptance
            for (unsigned int i=0; i<npix; i++) 
            {
                double nEvents = 0.0;
                double nBkg = 0.0; 
            
                if ( !(Emap0[i] > 0.0) )
                        continue;

                for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                {
                    int j;
                    j = illh::eq2loc_idx(i, timeidx, lat, lon, nTimesteps, CRmap);

                    nEvents += (*Nmap[timeidx])[i];
                    nBkg += norm[timeidx]*(CRmap[j]+diffCRmap[j]);
                }
            
                if ((Emap0[i] > 0.0) && (nBkg > 0.0)) 
                {
                    Emap[i] = nEvents/nBkg;
                } else
                { 
                    Emap[i] = 0.0;
                }
            }
            
            // Renormalize
            double sumEmap = illh::Sum(Emap, npix);
            
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                norm[timeidx] = norm[timeidx]*sumEmap; 
            } 
            for (unsigned int i=0; i < npix;i++ ) 
            { 
                    Emap[i] = Emap[i]/sumEmap;
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
            namefits << foldername << boost::format("/CR_%s_%d_%d_iteration%02d.fits") % detector % nsideOut % nTimesteps % iteration;
            if (fs::exists(namefits.str()) ) { 
                    fs::remove( namefits.str() ); 
            }
            fitshandle fitsOut; 
            fitsOut.create(namefits.str().c_str()); 
            write_Healpix_map_to_fits(fitsOut, diffCRmapNormed, MyDTYPE);
            fitsOut.close(); 

       
            // write N_tau^(n) normalizaton to file
            illh::save_iter(foldername, norm, Emap, detector, nsideOut, nTimesteps, iteration);

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
                            j = illh::loc2eq_idx(i, timeidx, lat, lon, nTimesteps, CRmap);

                            // global significance 
                            if (Emap0[j]*norm0[timeidx] > 0.0 && Emap[j]*norm[timeidx] > 0.0) { 
                                    significancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*Emap[j]*norm[timeidx]; 
                                    significancemap[i] += +2.0*CRmap[i]*Emap0[j]*norm0[timeidx]; 
                                    double temp1 = Emap[j]/Emap0[j]*norm[timeidx]/norm0[timeidx];
                                    significancemap[i] += 2.0*(*Nmap[timeidx])[j]*log(temp1*(1.0+diffCRmap[i]/CRmap[i])); 
                            }
                    } 
            }

            //write S_a^(n)
            stringstream nameSIGfits;
            nameSIGfits << foldername <<  boost::format("/significance_%s_%d_%d_iteration%02d.fits") % detector % nsideOut % nTimesteps % iteration; 
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
    bkgmapname <<  foldername << boost::format("/background_%s_%d_%d.fits.gz") % detector % nsideOut % nTimesteps;
    if (fs::exists(bkgmapname.str()) ) {
            fs::remove( bkgmapname.str() );
    }
    bkgfitsOut.create(bkgmapname.str().c_str());
    write_Healpix_map_to_fits(bkgfitsOut, bkgMap, MyDTYPE);
    bkgfitsOut.close();

    return 0;
}
