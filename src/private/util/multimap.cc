#include <util/multimap.h>


std::vector<SkyMapPtr>
read_group_hdu(fitsfile *fptr, int hdu_index)
{
    int status = 0;
    int hdutype = 0;

    // Move to the requested HDU
    fits_movabs_hdu(fptr, hdu_index, &hdutype, &status);
    if (status) {
        fits_report_error(stderr, status);
        throw std::runtime_error("Could not move to HDU");
    }

    if (hdutype != IMAGE_HDU)
        throw std::runtime_error("HDU is not an image");

    int naxis = 0;
    long naxes[2] = {0, 0};

    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 2, naxes, &status);

    if (status) {
        fits_report_error(stderr, status);
        throw std::runtime_error("Could not read image dims");
    }

    if (naxis != 2)
        throw std::runtime_error("Image is not 2-D (expected layers × pixels)");

    long npix  = naxes[0];
    long nmaps = naxes[1];

    long total = npix * nmaps;

    std::vector<double> buffer(total);

    long fpixel[2] = {1,1};
    fits_read_pix(fptr, TDOUBLE, fpixel, total,
                  nullptr, buffer.data(), nullptr, &status);

    if (status) {
        fits_report_error(stderr, status);
        throw std::runtime_error("Error reading image data");
    }

    // Deduce NSIDE from pixel count
    SkyMap tmpmap;
    int nside = tmpmap.npix2nside(npix);

    // Reconstruct maps
    std::vector<SkyMapPtr> maps(nmaps);
    for (int m = 0; m < nmaps; ++m) {
        maps[m] = SkyMapPtr(new SkyMap);
        maps[m]->SetNside(nside, RING);  // or NESTED
        for (long p = 0; p < npix; ++p)
           (*maps[m])[p] = buffer[m * npix + p];
    }

    return maps;
}

void write_map_group(
        fitsfile *fptr,
        const std::vector< SkyMap > &maps,
        const std::string &extname,
        int &status)
{
    long npix = maps[0].Npix();
    long nmaps = maps.size();

    // Create a new HDU
    long buffer_size[2] = {npix, nmaps};
    fits_create_img(fptr, DOUBLE_IMG, 2,
                    buffer_size, &status);

    // Add EXTNAME keyword
    fits_update_key(fptr, TSTRING, "EXTNAME",
                    (void*)extname.c_str(),
                    "Name of this image extension",
                    &status);

    // Flatten maps into a buffer
    std::vector<double> buffer(nmaps * npix);
    for (int m = 0; m < nmaps; ++m)
        for (long p = 0; p < npix; ++p)
            buffer[m * npix + p] = maps[m][p];

    long fpixel[2]= {1,1};  // first pixel index
    fits_write_pix(fptr, TDOUBLE, fpixel,
                   npix * nmaps,
                   buffer.data(), &status);
}


void write_map_group(
        fitsfile *fptr,
        const std::vector< SkyMapPtr > maps,
        const std::string &extname,
        int &status)
{
    long npix = maps[0]->Npix();
    long nmaps = maps.size();

    // Create a new HDU
    long buffer_size[2] = {npix, nmaps};
    fits_create_img(fptr, DOUBLE_IMG, 2,
                    buffer_size, &status);

    // Add EXTNAME keyword
    fits_update_key(fptr, TSTRING, "EXTNAME",
                    (void*)extname.c_str(),
                    "Name of this image extension",
                    &status);

    // Flatten maps into a buffer
    std::vector<double> buffer(nmaps * npix);
    for (int m = 0; m < nmaps; ++m)
        for (long p = 0; p < npix; ++p)
            buffer[m * npix + p] = (*maps[m])[p];

    long fpixel[2]= {1,1};  // first pixel index
    fits_write_pix(fptr, TDOUBLE, fpixel,
                   npix * nmaps,
                   buffer.data(), &status);
}

std::vector<SkyMapPtr>
add_multimaps(
        std::vector<SkyMapPtr>& mm1, 
        std::vector<SkyMapPtr>& mm2)
{
    int nside = mm2[0]->Nside();
    int npix = mm2[0]->Npix();

    unsigned int t_index(0);
    if ( mm1.size()<1 ) {
        for (unsigned int t_index=0;t_index<mm2.size();t_index++) 
        {
            mm1.push_back(SkyMapPtr(new SkyMap));
            mm1.back()->SetNside(nside, RING);     
            mm1.back()->fill(0.);     
        }
    } else if ( mm1.size() != mm2.size() ) {
        throw std::runtime_error(
        "The two multimaps don't have the same number of time bins");
    } 
    if ( mm1[0]->Nside() != mm2[0]->Nside() ) {
        throw std::runtime_error(
        "The two multimaps don't have the same NSide");
    } 
    //std::cout << "Maps are Nside= "<<nside <<std::endl;

    for (unsigned int t_index=0;t_index<mm2.size();t_index++)
    {
        // Add pixel values 
        for (unsigned int ipix=0;ipix<npix;ipix++)
        {
                (*mm1[t_index])[ipix]+=(*mm2[t_index])[ipix];
        }
    }
    return mm1;
}


