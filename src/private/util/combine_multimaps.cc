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

    bool init(false);
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
            std::cout << "Combine several multi-layered FITS image files. ";
            std::cout << "If 'weights' option is used, ";
            std::cout << "it is assumed that two additional HDUs included: " << std::endl;
            std::cout << "- 1 contains weights. " << std::endl;
            std::cout << "- 2 contains weights squared to calculate statistical errors. ";
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

    int nside = -1, npix = -1;
    int status = 0;

    std::vector<SkyMapPtr> maps;
    std::vector<SkyMapPtr> wmaps;
    std::vector<SkyMapPtr> w2maps;
    

    for (const auto& fname : inputFiles) {
        //std::string fname = argv[i];
        fitsfile *fptr;

        fits_open_file(&fptr, fname.c_str(), READONLY, &status);
        if (status) { 
            fits_report_error(stderr, status); 
            throw std::runtime_error("Problem reading " + fname+ ".");
        }
        std::vector<SkyMapPtr> in_maps = read_group_hdu(fptr, 2);
        maps = add_multimaps( maps, in_maps);

        if (weights) {
            std::vector<SkyMapPtr> in_maps = read_group_hdu(fptr, 3);
            wmaps = add_multimaps( wmaps, in_maps);
            std::vector<SkyMapPtr> in_w2maps = read_group_hdu(fptr, 4);
            w2maps = add_multimaps( w2maps, in_maps);
        }
        fits_close_file(fptr, &status);
    }

    int nmaps = maps.size();
    if (nmaps>0) {
        npix  = maps[0]->Npix();
        nside  = maps[0]->Nside();
    }

    std::cout << "Read " << nmaps << " maps with NSIDE=" 
              << nside << " NPIX=" << npix << std::endl;

    // -------------------------------------------------------------
    // Create output FITS file
    // -------------------------------------------------------------
    fitsfile *fptr = nullptr;

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

