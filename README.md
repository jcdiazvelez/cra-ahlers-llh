# cra-ahlers-llh
* Iterative log-likelyhood reconstruction for CRA
* C++ implementation of maximum-likelihood technique for reconstructing cosmic-ray anisotropy maps
* Cite: [M. Ahlers et al 2016 ApJ 823 10](http://iopscience.iop.org/article/10.3847/0004-637X/823/1/10)

## Notes: Note the change in format on the output fits files!

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


## Running:

You need to provide a JSON configuration file that lists the parameters for each observatory. There is an ``examples`` directory that contains configurations for IceCube, HAWC and IceCube+HAWC. The config file can contain an arbitrary number of observatories. 

The output will be in the form of a HEALPix map of NSide = 64 but you can specify a different NSide. If you don't specify the number of iterations, it will default to 20. 

Other parameters such as ``--seed``, ``--iso``, and ``--fluctuate`` are only used for systematic studies.

````{verbatim}
  ./bin/multi-llh [options] --config <json config file>  -o <output directory> --iterations <number of iterations> 
Basic Command Line Parameter App
Options:
  -h [ --help ]                    produce help message
  -o [ --outdir ] arg (=./sample/) Directory of output
  --nsideout arg (=64)             Healpix Nside for output map
  --timesteps arg (=360)           Number of time steps
  --timestepmin arg (=0)           First time step to use
  --timestepmax arg (=0)           Last time step to use
  --iterations arg (=20)           Number of iterations
  --save-iter                      save each iteration
  -f [ --fluctuate ]               add random fluctuations
  --seed arg (=123)                RNG seed
  --iso                            make isotropic map
  --smoothing-radius arg (=1)      Smoothing radius (deg) for significance map
  --config arg (=config.json)      JSON config
````

**Input files (new version):**
The input should be specified in the JSON config file (see examples) and consists of data binned into local HEALPix maps for each sidereal time bin (24h/360 bins by default). The config file specifies an input `multimap` file per detector. The file is a multi-layered FITS file:
````{verbatim}
  "multimap": "CR_ICECUBE_LOCAL_NSIDE64.fits.gz",
````
For more information about multimaps see section below^1.


**Input files (original version):**
The input should be specified in the JSON config file (see examples) and consists of data binned into local HEALPix maps for each sidereal time bin (24h/360 bins by default). The config file specifies the prefix and suffix of each file. For eaxmple:

````{verbatim}
  "prefix": "CR_ICECUBE_LOCAL_NSIDE64_degbin-",
  "suffix": ".fits.gz"
````
would assume that the files are named 
````{verbatim}
  CR_ICECUBE_LOCAL_NSIDE64_degbin-000.fits.gz
  CR_ICECUBE_LOCAL_NSIDE64_degbin-001.fits.gz
  ...
  CR_ICECUBE_LOCAL_NSIDE64_degbin-359.fits.gz
````

**Local coordinates**:

The rotation matrices implemented in this code use a conversion of [J2000.0 equatorial coordinates](http://en.wikipedia.org/wiki/Epoch_(astronomy)) to local coordinates.

Local coordinates are expressed in a right-handed sense, with $x=E$ and $y=N$, with $\phi$ = 0 oriented along the x-axis and following 
the [physics convention](https://en.wikipedia.org/wiki/Spherical_coordinate_system) for a unit sphere with coordinates $(\theta,\phi)$, with polar angle $\theta$ (theta) (angle with respect to positive polar axis), and azimuthal angle $\phi$ (phi).


**Output files:**
***Note file formats have changed:***
The output files will be named:
* ```"CR_{detector(s)}__64_360_iteration{iter}.fits.gz"```: ($a_{\ell m}$-smoothed relative intensity, relative error $\sigma_i/\mu_i$, signed Li-Ma significance $\sigma_{i,\mathrm{Li-Ma}}$)
* ```"RAW_{detector(s)}_{nside}_{nTbins}_iteration{iter}.fits.gz"```: (data, background, expectation $\mu$)
* ```"llhratio_{detector(s)}_{nside}_{nTbins}_iteration{iter}.dat"```: ($\mathrm{llh}_n - \mathrm{llh}_0$)
* ```"exposure_{detector}_{nside}_{nTbins}_iteration{iter}.fits.gz"```: ($A_i$: one file per dectector)
* ```"norm_{detector}_{nside}_{nTbins}_iteration{iter}.dat"```: (<math>$N_\tau$</math>: one file per detector)

If ```--save-iter``` option is used, there will be one of each of the above files per iteration. Otherwise only the initial and final results will be saved.

The code will terminate after convergence or after reaching the number of
iterations (whichever comes first). 


## $^1$ MultiMaps:
Instead of dealing with N time-binned FITS files, cra-ahlers now supports a
multi-layered image FITS format consisting of N layers. If this is an
unweighted datasets, the multimap will have a single HDU called "MAPS".
Additional HDUs are included for weighted data: "WMAP" is a set of weight sums
per pixel and "W2MAPS" is a set of squared weight sums corresponding to the statistical
variance.

### Multi-map utilities
There are a couple of tools included to deal with multimaps:
1. create_multimap - will combined a set of healpix map files into a single
   multimap.
2. combine_multimaps - will sum a set of multimaps (element wise) to combine
   into a single file.




## JSON Config:
The following is an example configuration for HAWC using a weighted multimap:
````{verbatim}
{
    "detectors": {
        "HAWC": { 
            "nside":64,
            "latitude": 19.0,
            "longitude": -97.0,
            "thetamax": 43.0,
            "weights": true,
            "mutimap": "test_multimap2.fits.gz",
        }
    }
}
````
Here you specified the nside as well as the location of the detector in
latitude and longitud coordinates. The "weights" flag indicates that cra-llh
should work usig the weighs rather than the counts.
