from astropy.io import fits
import numpy as np
import healpy

def read_three_multilayer_sets(filename):
    """
    Read back the file written by write_three_multilayer_sets().
    
    Returns:
        local_maps, weighted_maps, w2_maps
    where each set is a list of 1D maps (numpy arrays).
    """

    hdul = fits.open(filename)

    set1_data = hdul["LMAPS"].data
    set2_data = hdul["WMAPS"].data
    set3_data = hdul["W2MAPS"].data

    # Convert 2D array back into list of 1D maps
    set1 = [row.copy() for row in set1_data]
    set2 = [row.copy() for row in set2_data]
    set3 = [row.copy() for row in set3_data]

    hdul.close()

    return set1, set2, set3



def write_three_multilayer_sets(filename,
                                local_maps,
                                weighted_maps,
                                w2_maps):
    """
    Write three sets of multi-layer HEALPix maps into a FITS file.
    
    Each set is a list of maps:
        setX_maps = [map1, map2, ..., mapN]
    Each map must be 1D array of size npix (= 12*nside^2)

    File structure:
        0: PRIMARY (empty)
        1: LMAPS (ImageHDU)
        2: WMAPS (ImageHDU)
        3: W2MAPS (ImageHDU)
    """

    # sanity check: consistent npix
    npix = len(local_maps[0])
    for s in (local_maps, weighted_maps, w2_maps):
        for m in s:
            if len(m) != npix:
                raise ValueError("All maps must have same npix")

    # convert lists → 2D arrays
    img1 = np.vstack(local_maps)
    img2 = np.vstack(weighted_maps)
    img3 = np.vstack(w2_maps)

    # Create HDUs
    primary = fits.PrimaryHDU()
    hdu1 = fits.ImageHDU(data=img1, name="LMAPS")
    hdu2 = fits.ImageHDU(data=img2, name="WMAPS")
    hdu3 = fits.ImageHDU(data=img3, name="W2MAPS")

    # Write
    hdul = fits.HDUList([primary, hdu1, hdu2, hdu3])
    hdul.writeto(filename, overwrite=True)

    print(f"Wrote file: {filename}")

    
def create_multimap(outfile,infile_list):
    """
    Convert set of FITS HEALPix files to multimap.
    
    Each file in list is a HEALPix map :
        setX_maps = [map1, map2, ..., mapN]
    

    Output File structure:
        0: PRIMARY (empty)
        1: LMAPS (ImageHDU)
        2: WMAPS (ImageHDU)
        3: W2MAPS (ImageHDU)
    """
    local_maps = []
    weighted_maps = []
    w2_maps = []
    for f in infile_list:
        m, w, w2 = healpy.read_map(f)
        local_maps.append(m)
        weighted_maps.append(w)
        w2_maps.append(w2)
        
    write_three_multilayer_sets(outfile,local_maps, weighted_maps, w2_maps)
        
        
