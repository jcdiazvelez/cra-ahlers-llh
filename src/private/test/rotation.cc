/**
 * cra-llh-ahlers
 *
 * @version $Id: $
 *
 * @date: $Date: $
 *
 * @author Juan Carlos Diaz-Velez <juan.diazvelez@alumnos.udg.mx>
 *
*/


#include <iter-lhreco-proj/illh-utils.h>
#include <cmath>

double angular_distance(pointing p0, pointing p1)
{
   return acos(cos(p0.theta)*cos(p1.theta) + sin(p0.theta)*sin(p1.theta)*cos(p0.phi-p1.phi));
}

int main(int argc, char* argv[])
{
  unsigned nside(64);
  if (argc > 1)
    nside = atoi(argv[1]);
  std::cout << "Using NSide = " << nside << std::endl;

  double latitude = 18.0*M_PI/180.;
  double longitude = -97.0*M_PI/180.;

  unsigned npix = 12*nside*nside;
  float pixel_area = 4*M_PI/npix;
  float pixel_size = sqrt(pixel_area);


  SkyMap CRmap;
  CRmap.SetNside(nside, RING);
  CRmap.fill(1.0);

  double crab_ra = (5.+  (34. + 32./60.)/60)/24.*2*M_PI;
  double crab_dec = (22. +52/3600)*M_PI/180.;
  double aries_ra = 0;
  double aries_dec = 0;

  pointing aries(M_PI - aries_ra, aries_ra);
  pointing crab(M_PI - crab_ra, crab_ra);

  int timeidx = 0;
  int nTimesteps = 360;

  //rotation from Equatorial (ra,dec) to local (theta,phi)
  auto integers = {0, 20, 101, 235, 300, 400, 500, 1000};
  int j,k;

  for (auto i : integers) 
  {
    if (i>=npix)
    {
        continue;
    }

    j = illh::eq2loc_idx(i, timeidx, latitude, longitude, nTimesteps, CRmap);
    k = illh::loc2eq_idx(j, timeidx, latitude, longitude, nTimesteps, CRmap);

    // get theta phi from pixels
    pointing pt0 = CRmap.pix2ang(i); 
    pointing pt1 = CRmap.pix2ang(k); 
    auto ang_dist = angular_distance(pt0,pt1);

    if (ang_dist > 2*pixel_size)  // tolerance is one pixel
    {
       std::cout << i << ", " << k << "diff = "<< ang_dist << ": "<< pixel_size <<std::endl;
       return 1;
    }
  }


  for (auto t = 0; t < nTimesteps; t++)
  {
    auto ak = CRmap.ang2pix(aries); 
    auto aj = illh::eq2loc_idx(ak, t, latitude, longitude, nTimesteps, CRmap);
    auto ai = illh::loc2eq_idx(aj, t, latitude, longitude, nTimesteps, CRmap);

    auto ck = CRmap.ang2pix(crab); 
    auto cj = illh::eq2loc_idx(ck, t, latitude, longitude, nTimesteps, CRmap);
    auto ci = illh::loc2eq_idx(cj, t, latitude, longitude, nTimesteps, CRmap);

    // get theta phi from pixels
    pointing pt0 = CRmap.pix2ang(ai); 
    pointing pt1 = CRmap.pix2ang(ci); 

    if (angular_distance(pt0,aries) > 2*pixel_size)  // tolerance is one pixel
    {
       std::cout << boost::format("Aries (%f,%f) -> (%f,%f)") % aries.theta % aries.phi % pt0.theta % pt0.phi << std::endl;
       return 1;
    }
    if (angular_distance(pt0,aries) > 2*pixel_size)  // tolerance is one pixel
    {
       std::cout << boost::format("Crab (%f,%f) -> (%f,%f)") % crab.theta % crab.phi % pt1.theta % pt1.phi << std::endl;
       return 1;
    }
  }
  return 0;
}

