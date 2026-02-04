#ifndef __MULTIMAP_H__

// Healpix
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map_fitsio.h>

// CFITSIO
#include <fitsio.h>
#include <memory>


typedef Healpix_Map<double> SkyMap; 
typedef std::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer



std::vector<SkyMapPtr>
read_group_hdu(fitsfile *fptr, int hdu_index);

void write_map_group(
        fitsfile *fptr,
        const std::vector< SkyMap > &maps,
        const std::string &extname,
        int &status);

void write_map_group(
        fitsfile *fptr,
        const std::vector< SkyMapPtr > aps,
        const std::string &extname,
        int &status);


std::vector<SkyMapPtr>
add_multimaps(
        std::vector<SkyMapPtr>& mm1, 
        std::vector<SkyMapPtr>& mm2);


#define __MULTIMAP_H__
#endif //__MULTIMAP_H__
