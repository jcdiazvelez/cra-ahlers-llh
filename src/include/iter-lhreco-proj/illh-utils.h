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

//#include <boost/shared_ptr.hpp>
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

#include <iter-lhreco-proj/pickle.h>
#include <util/multimap.h>


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




//typedef Healpix_Map<double> SkyMap; 
//typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer

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

    bool multifits;
    bool weights;

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




namespace illh
{
    /*
     * Sum over all pixels
     *
     */
    double 
    Sum(const SkyMap& map,unsigned int npix);


#if __cplusplus > 199711L
    void fluctuate(std::vector<SkyMapPtr>& Nmap, boost::mt19937 rng);



    /*
     * Generate isotropic maps with Poisson noise
     *
     */
    void isotropic(std::vector<SkyMapPtr>& Nmap, boost::mt19937 rng);
#endif

    /*
     * loadMap - read local maps and generate vector of maps bined in 
     * siderial time steps
     */
    void loadMap(
        std::vector<SkyMapPtr>& Nmap,
        unsigned nTimesteps, 
        unsigned nsideIn, 
        unsigned nsideOut, 
        std::string filename, 
        bool weights=false);



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
        std::string suffix, 
        bool weights=false);

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
                    unsigned iteration);

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
    loc2eq_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, const SkyMap& CRmap);

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
    eq2loc_idx(unsigned i, unsigned timeidx, double lat, double lon, unsigned nTimesteps, const SkyMap& CRmap);



    void mlh_iteration(
        unsigned int pixmin, 
        unsigned int pixmax, 
        unsigned int nTimesteps,
        const std::vector<DetectorPtr> detectors,
        SkyMap &diffCRmap,
        SkyMap &CRmap,
        SkyMap &bkgMap,
		std::mutex &pmutex);

    void mlh_normalization(
        unsigned int timeidxMin, 
        unsigned int timeidxMax, 
        unsigned int nTimesteps, 
        std::vector<DetectorPtr> detectors,
        SkyMap &diffCRmap,
        SkyMap &CRmap,
		std::mutex &pmutex);

    void mlh_significance(
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
		std::mutex &pmutex);

    void mlh_smooth(
        unsigned int pixmin, 
        unsigned int pixmax, 
        unsigned int nTimesteps,
        double smoothing_radius,
        SkyMap &significancemap,
        SkyMap &muon,
        SkyMap &muoff,
        SkyMap &Na,
		std::mutex &pmutex);

    void mlh_acceptance(
        unsigned int pixmin, 
        unsigned int pixmax, 
        unsigned int nTimesteps,
        std::vector<DetectorPtr> detectors,
        const SkyMap &diffCRmap,
        const SkyMap &CRmap,
		std::mutex &pmutex);


}
