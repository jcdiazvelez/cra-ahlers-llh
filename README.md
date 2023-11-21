# cra-ahlers-llh
* Iterative log-likelyhood reconstruction for CRA
* C++ implementation of maximum-likelihood technique for reconstructing cosmic-ray anisotropy maps
* Cite: [M. Ahlers et al 2016 ApJ 823 10](http://iopscience.iop.org/article/10.3847/0004-637X/823/1/10)

## scripts
Collection of scripts for driving production of extraction, reconstruction and analysis of cosmic-ray data.


## Installation


**Prerequisites**:

iter-lhreco dependencies:

* CFITSIO: depends on cfitsio: https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html 
* GNU Scientific Library: depends on gsl: https://www.gnu.org/software/gsl/doc/html/index.html
* HEALpix: depends on healpix: http://healpix.sourceforge.net/
* BOOST C++ libraries: https://www.boost.org/

You will also need:

* [CMake](https://cmake.org) >= 3.1
* A C++11-compliant compiler (e.g. [gcc](https://gcc.gnu.org) >= 4.8.1 or [clang](https://clang.llvm.org) >= 3.3)


**Installation**:

The C++ projects are built using CMake.
Each project has a CMakeList.txt that will detect dependencies and generate a MakeFile. To do this cd into the build directory and excecute the commands:

  cmake ../src;
  make


**Running**:
You need to provide a JSON configuration file that lists the parameters for each observatory. There is an ``examples`` directory that contains configurations for IceCube, HAWC and IceCube+HAWC. The config file can contain an arbitrary number of observatories. 

````{verbatim}
  ./bin/multi-llh [options] --config <json config file>  -o <output file> --iterations <number of iterations> 
  Basic Command Line Parameter App
  Options:
    -h [ --help ]                    produce help message
    -o [ --outdir ] arg (=./sample/) Directory of output
    --nsideout arg (=64)             Healpix Nside for output map
    --timesteps arg (=360)           Number of time steps
    --timestepmin arg (=0)           First time step to use
    --timestepmax arg (=0)           Last time step to use
    --iterations arg (=20)           Number of iterations
    -f [ --fluctuate ]               add random fluctuations
    --seed arg (=123)                RNG seed
    --iso                            make isotropic map
    --config arg (=config.json)      JSON config
````
