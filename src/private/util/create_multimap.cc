#include <iostream>
#include <vector>
#include <string>


// Healpix
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map_fitsio.h>

//Boost
#include <boost/program_options.hpp> 

// CFITSIO
#include <fitsio.h>

#include <util/utils.h>
#include <util/multimap.h>


namespace po = boost::program_options; 



int main(int argc, char **argv)
{


    bool compressed(false);
    bool weights;
    std::string output_file;
    std::vector<std::string> inputFiles;


    po::options_description desc("Options"); 
    po::variables_map vm; 

    try { 
       // Define and parse the program options 
        desc.add_options() 
             ("help,h", "produce help message") 
             ("weights,w", po::bool_switch(&weights)->default_value(false), "Includes weights column")
			 ("output-file,o", po::value<std::string>(&output_file), "Output FITS multimap")
			 ("input-files", po::value<std::vector<std::string>>(&inputFiles), "Input files");

        po::positional_options_description pos;
        pos.add("input-files", -1); // -1 means "consume all remaining args"
     
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm); 
         
        /// --help option 
        if ( vm.count("help")  ) { 
            std::cout << argv[0] << std::endl;
            std::cout << "Create multi-layered FITS image file from multiple time-binned HEALpix FITS files. ";
            std::cout << "Column 0 of each input file should be a time-binned map of counts in local coordinates. ";
            std::cout << "If 'weights' option is used, ";
            std::cout << "it is assumed that two additional columns are included: " << std::endl;
            std::cout << "- column 1 contains weights. " << std::endl;
            std::cout << "- column 2 contains weights squared to calculate statistical errors. ";
            std::cout << std::endl;
            std::cout << "Output file should have either .fits or .fits.gz sufix. ";
            std::cout << "In the latter case, gzip-2 compression will be used.";
            std::cout << std::endl;
            std::cout << "Allowed options" 
                      << std::endl << desc << std::endl; 
            return 0; 
        } 
        po::notify(vm); // throws on error, so do after help in case there are any problems

    } catch(po::error& e) { 
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
            std::cerr << desc << std::endl; 
            return 1; 
        
    } catch(std::exception& e) { 
            std::cerr << "Unhandled Exception reached the top of main: " 
                      << e.what() << ", application will now exit" << std::endl; 
            return 2; 
    } 

    if ( !ends_with(output_file, ".fits") && !ends_with(output_file, ".fits.gz")  ) { 
        po::notify(vm); // throws on error, so do after help in case there are any problems
        std::cerr << "specified output file does not have the correct suffix. ";
        std::cerr << "Must be either '.fits' or '.fits.gz'" << std::endl;
        return 3; 
    } 



    std::vector<Healpix_Map<double>> maps ;
    std::vector<Healpix_Map<double>> wmaps;
    std::vector<Healpix_Map<double>> w2maps;
    int nside = -1, npix = -1;

    //for (int i = 2; i < argc; ++i) {
    for (const auto& fname : inputFiles) {
        //std::string fname = argv[i];
        std::cout << maps.size() << ":" <<fname << std::endl;
        fitshandle fh;
        fh.open(fname.c_str());
        fh.goto_hdu(2);

        Healpix_Map<double> map;
        Healpix_Map<double> wmap;
        Healpix_Map<double> w2map;

        read_Healpix_map_from_fits(fh, map, 1);
        if (weights) {
            read_Healpix_map_from_fits(fh, wmap, 2);
            read_Healpix_map_from_fits(fh, w2map, 3);
        }

        if (nside < 0) {
            nside = map.Nside();
            npix  = map.Npix();
        } else if (map.Nside() != nside) {
            std::cerr << "Error: NSIDE mismatch" << std::endl;
            return 1;
        }
 

        if (nside < 0) {
            nside = map.Nside();
            npix  = map.Npix();
        } else if (map.Nside() != nside) {
            std::cerr << "Error: NSIDE mismatch" << std::endl;
            return 1;
        }

        maps.push_back(map);
        if (weights) {
            wmaps.push_back(wmap);
            w2maps.push_back(w2map);
        }
        fh.close();
    }

    int nmaps = maps.size();

    std::cout << "Read " << nmaps << " maps with NSIDE=" 
              << nside << " NPIX=" << npix << std::endl;

    // -------------------------------------------------------------
    // Create output FITS file
    // -------------------------------------------------------------
    fitsfile *fptr = nullptr;
    int status = 0;

    fits_create_file(&fptr, ("!" + output_file).c_str(), &status);
    if (ends_with(output_file, ".fits.gz") ) {
        compressed = true;
        fits_set_compression_type(fptr, GZIP_2, &status);
    }   

    // Primary HDU (empty image)
    long naxes[2] = {0,0};
    fits_create_img(fptr, DOUBLE_IMG, 0, naxes, &status);
    if (status) { fits_report_error(stderr, status); return status; }

    // -------------------------------------------------------------
    // Write 3 groups as separate FITS image HDUs
    // -------------------------------------------------------------
    unsigned int columns(1);
    write_map_group(fptr, maps, "LMAPS", status);
    if (weights) {
        write_map_group(fptr, wmaps, "WMAPS", status);
        write_map_group(fptr, w2maps, "W2MAPS", status);
        columns = 3;
    }

    fits_close_file(fptr, &status);

    if (status)
        fits_report_error(stderr, status);
    else
        std::cout << "Wrote " << columns << " multilayer images "
                  << (weights ? "(counts, weights, weights^2) " : "(counts) ") 
                  << "each with " << maps.size() << " layers to "
                  << (compressed ? "compressed " : "") << "file "
                  << output_file << std::endl;

    return status;
}

