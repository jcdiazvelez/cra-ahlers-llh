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
#include <hawcnest/HAWCUnits.h>
#include <hawcnest/Logging.h>
#include <hawcnest/CommandLineConfigurator.h>
#else
namespace po = boost::program_options;
#define log_info(args...) cout << args << endl
#endif


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
    std::vector< std::string > input;
    double smoothing;
    int iteration;;
    int nTimesteps;
#ifdef HAWCNEST
    CommandLineConfigurator cl;

    cl.AddPositionalOption<std::vector <std::string> >("input", "Data Map, Differential Relative, Intensity Map, Local Acceptance Map final, Local Acceptance Map initial, All Sky Exposure final, All Sky Exposure initial");
    cl.AddOption<std::string>("outdir","./sample/","Directory of output");
    cl.AddOption<double>("smoothing,s",0.,"Smoothing [degree]");
    cl.AddOption<int>("iteration,i",20,"Iteration number");
    cl.AddOption<int>("timesteps,t",360,"Number of time steps");
    if (!cl.ParseCommandLine(argc, argv)) {
      return 1;
    }

    input = cl.GetArgument<std::vector <std::string> >("input");
    foldername = cl.GetArgument<std::string>("outdir");
    smoothing = cl.GetArgument<double>("smoothing");
    iteration = cl.GetArgument<int>("iteration");
    nTimesteps = cl.GetArgument<int>("timesteps");
#else

    po::options_description desc("Options");
    po::positional_options_description p;
    p.add("input", -1);
    po::variables_map vm;

    try {
       // Define and parse the program options 
       desc.add_options()
             ("help", "produce help message")
             ("input", po::value<std::vector <std::string> >(&input)->multitoken(), "Input files")
             ("outdir", po::value<std::string>(&foldername)->default_value("./sample/"), "Directory of output")
             ("smoothing", po::value<double>(&smoothing)->default_value(0.), "Smoothing [degrees]")
             ;

        //po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw 
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

        /// --help option 
        if ( vm.count("help")  ) {
            std::cout << "Basic Command Line Parameter App"
                      << std::endl << desc << std::endl;
            return 0;
        }
        po::notify(vm); // throws on error, so do after help in case there are any problems

        if (vm.count("input"))
        {
            cout << "Input files: " << "\n";
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

    // output directory
    fs::path dir(foldername);
    if(!(fs::exists(dir)) ) {
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir))
                    log_info("....successfully created !");
    }

    SkyMap dat;
    read_Healpix_map_from_fits(input[0], dat,       1);
    log_info("Loaded data map " << input[0]);
    SkyMap diffCRmap;
    read_Healpix_map_from_fits(input[1], diffCRmap, 1);
    log_info("Loaded diff rel int map " << input[1]);
    int npix = dat.Npix();
    SkyMap Emap; 
    read_Healpix_map_from_fits(input[2], Emap,      1);
    std::vector<double> norm;
    std::ifstream filein(input[3].c_str());
    for (std::string line; std::getline(filein, line); )
    {
            norm.push_back( std::atof( line.c_str() ) );
    }


    // Get large scale map
    SkyMap largeDiffCRmap;
    largeDiffCRmap.SetNside( dat.Nside(), dat.Scheme());
    largeDiffCRmap.fill(0.);
//

//

//    const lon = -97.*HAWCUnits::degree;
//    // calculate signifiance : S_a^(n)
//    for (unsigned int timeidx=0;timeidx < nTimesteps;timeidx++)
//    {
//    
//        double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon;
//        double cb = cos(beta);
//        double sb = sin(beta);
//    
//        vec3 vp;
//    
//        // rotation matrix 
//        matrix<double> rot (3, 3);
//    
//        rot(0,0) = cb*slat;
//        rot(0,1) = sb*slat;
//        rot(0,2) = -clat;
//    
//        rot(1,0) = -sb;
//        rot(1,1) = cb;
//        rot(1,2) = 0;
//    
//        rot(2,0) = cb*clat;
//        rot(2,1) = clat*sb;
//        rot(2,2) = slat;
//    
//        for (unsigned int i=0; i < npix;i++ )
//        {
//            std::vector<int> listOfPixels;
//            dat.query_disc(pt_i, smoothRad, listOfPixels);
//    
//            double sumEmap0 = 0.;
//            double sumEmap0 = 0.;
//            double sumEmap0 = 0.;
//    
//            for (unsigned int p=0; p<listOfPixels.size(); p++){
//                   
//                 int index = listOfPixels[p];
//    
//                 //rotation from Equatorial (ra,dec) to local
//                 vec3 v =  CRmap.pix2vec(index);
//                 vp.x = rot(0,0)*v.x+rot(0,1)*v.y+rot(0,2)*v.z;
//                 vp.y = rot(1,0)*v.x+rot(1,1)*v.y+rot(1,2)*v.z;
//                 vp.z = rot(2,0)*v.x+rot(2,1)*v.y+rot(2,2)*v.z;
//
//                 unsigned int j = dat.vec2pix(vp); // local pixel 
//            
//            }
//        }
//    } 

    SkyMap bkg;
    bkg.SetNside( dat.Nside(), dat.Scheme());
    bkg.fill(0.);
    // Smoothed Bkg : B_a^(n)
    SkyMap bkgSmooth;
    bkgSmooth.SetNside( dat.Nside(), dat.Scheme());
    bkgSmooth.fill(0.);
    // Smoothed Data : D_a^(n)
    SkyMap datSmooth;
    datSmooth.SetNside( dat.Nside(), dat.Scheme());
    datSmooth.fill(0.);
    // Smoothed Diff. Rel Int : dI_a^(n)
    SkyMap diffCRmapSmooth;
    diffCRmapSmooth.SetNside( dat.Nside(), dat.Scheme());
    diffCRmapSmooth.fill(0.);
    // Smoothed Significance : S_a^(n)
    SkyMap significance;
    significance.SetNside( dat.Nside(), RING);
    significance.fill(0.);
    //
    SkyMap simpleSignificance;
    simpleSignificance.SetNside( dat.Nside(), RING);
    simpleSignificance.fill(0.);
  
    // Clean up diffCRmap of unphysical values
    for (unsigned int i=0; i<npix;i++) {
            if ( (diffCRmap[i] < -1.) || diffCRmap[i] > 1.) { 
                    diffCRmap[i] = 0.;
                    dat[i] = 0.;
            }
    }
    log_info("Smoothing " << smoothing << " degrees");
    for (unsigned int i=0; i < npix;i++ )
    {
            pointing pt_i = dat.pix2ang(i);

            std::vector<int> listOfPixels;
            dat.query_disc(pt_i, smoothing*3.14159/180., listOfPixels);

            double sumDat = 0.;
            double sumBkg = 0.;

            for (unsigned int p=0; p<listOfPixels.size(); p++){
                   
                    int index = listOfPixels[p];
                    sumDat += dat[index];
                    sumBkg += dat[index]/( 1. + diffCRmap[index] );
            }
            bkg[i] = dat[i] / (1. + diffCRmap[i]);
            bkgSmooth[i] = sumBkg;
            datSmooth[i] = sumDat;
            if (sumDat > 0. && sumBkg > 0.) {
                diffCRmapSmooth[i] = sumDat/sumBkg - 1.;
                // WRONG
                simpleSignificance[i] = sumDat / sqrt(sumBkg);
            } else {
                diffCRmapSmooth[i] = 0.;
                simpleSignificance[i] = 0.;
            }

    }

    int nsideOut = dat.Nside();

    fitshandle fitsOut;
    // write B_a^(n) smoothed
    stringstream nameBkgfits;
    nameBkgfits << foldername << boost::format("/bkg_HAWC_%d_%d_iteration%02d_S%02d.fits") % nsideOut % nTimesteps % iteration % smoothing;
    if (fs::exists(nameBkgfits.str()) ) {
            fs::remove( nameBkgfits.str() );
    }
    fitsOut.create(nameBkgfits.str().c_str());
    write_Healpix_map_to_fits(fitsOut, bkgSmooth, PLANCK_FLOAT64);
    fitsOut.close();


    // write I_a^(n) smoothed
    stringstream namefits;
    namefits << foldername << boost::format("/CR_HAWC_%d_%d_iteration%02d_S%02d.fits") % nsideOut % nTimesteps % iteration % smoothing;
    if (fs::exists(namefits.str()) ) {
            fs::remove( namefits.str() );
    }
    fitsOut.create(namefits.str().c_str());
    write_Healpix_map_to_fits(fitsOut, diffCRmapSmooth, PLANCK_FLOAT64);
    fitsOut.close();

    // write S_a^(n) smoothed
    stringstream nameSigfits;
    nameSigfits << foldername << boost::format("/Significance_HAWC_%d_%d_iteration%02d_S%02d.fits") % nsideOut % nTimesteps % iteration % smoothing;
    if (fs::exists(nameSigfits.str()) ) {
            fs::remove( nameSigfits.str() );
    }
    fitsOut.create(nameSigfits.str().c_str());
    write_Healpix_map_to_fits(fitsOut, significance, PLANCK_FLOAT64);
    fitsOut.close();


    // calculate signifiance : S_a^(n)
//    for (unsigned int timeidx=0;timeidx < nTimesteps;timeidx++)
//    {
//
//          double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon;
//          double cb = cos(beta);
//          double sb = sin(beta);
//
//          vec3 vp;
//
//          // rotation matrix 
//          matrix<double> rot (3, 3);
//
//          rot(0,0) = cb*slat;
//          rot(0,1) = sb*slat;
//          rot(0,2) = -clat;
//
//          rot(1,0) = -sb;
//          rot(1,1) = cb;
//          rot(1,2) = 0;
//
//          rot(2,0) = cb*clat;
//          rot(2,1) = clat*sb;
//          rot(2,2) = slat;
//
//          for (unsigned int i=0; i < npix;i++ )
//          {
//                  std::vector<int> listOfPixels;
//                  dat.query_disc(pt_i, smoothRad, listOfPixels);
//
//                  double sumEmap0 = 0.;
//                  double sumEmap0 = 0.;
//                  double sumEmap0 = 0.;
//
//                  for (unsigned int p=0; p<listOfPixels.size(); p++){
//                         
//                          int index = listOfPixels[p];
//
//                          //rotation from Equatorial (ra,dec) to local
//                          vec3 v =  CRmap.pix2vec(index);
//                          vp.x = rot(0,0)*v.x+rot(0,1)*v.y+rot(0,2)*v.z;
//                          vp.y = rot(1,0)*v.x+rot(1,1)*v.y+rot(1,2)*v.z;
//                          vp.z = rot(2,0)*v.x+rot(2,1)*v.y+rot(2,2)*v.z;
//
//                          unsigned int j = CRmap.vec2pix(vp); // local pixel 
//                
//                  }
//
//                  // global significance 
//                  if (Emap0[j]*norm0[timeidx] > 0.0 && Emap[j]*norm[timeidx] > 0.0) {
//                          significancemap[i] +=
//                                  -2.0*(diffCRmap[i]+CRmap[i])*Emap[j]*norm[timeidx];
//
//                          significancemap[i] +=
//                                  +2.0*CRmap[i]*Emap0[j]*norm0[timeidx];
//
//                          double temp1 = norm[timeidx]/norm0[timeidx]*Emap[j]/Emap0[j];
//
//                          significancemap[i] +=
//                                  2.0*(*Nmap[timeidx])[j]*log(temp1*(1.0+diffCRmap[i]/CRmap[i]));
//                  }
//                  }
//          }
    return 0;
}
