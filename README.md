# DiscretePowerlawFitter
This is a fast discrete power-law fitter based on the paper `Power-law distributions in empirical data` from Aaron Clauset, Cosma Rohilla Shalizi, M. E. J. Newman and the [official implementations](https://aaronclauset.github.io/powerlaws/).

This library exists to address the long waiting time of the bootstrap procedure in the official implementations.

## Features
- Fast fast! Optimizations have been made to be able to obtain a goodness-of-fit p-value with 10,000 replicas in 1 second.
- Precise generation of replicas. Instead of discretizing the values of the random number generator for the continuous case, it uses the pseudo CDF inversion with binary search algorithm as described in `Appendix D` of the paper.
- Multiple interfaces: CLI, Mathlink and shared library.


## CLI Usage
```
INSTRUCTIONS: PowerLawFitterCpp [options]

  -d <data_to_test>,       --data=<data_to_test>            
                    Sample data as a list of comma-separated integers.
  -r <number_of_replicas>, --replicas=<number_of_replicas>  
                    Number of bootstrap replicas. Default is 2000.
  -x <xMin>,               --x_min=<xMin>                   
                    Known value of xMin if there is any.
  -s,                      --single_thread                  
                    Use only one thread for the boot-strapping.
                           --help                           
                    Show instructions.
```

### Examples
Example command:
```
./PowerLawFitterCpp --data=1,2,1,1,1,1,1,1,3,2,1,2,2,1,5,1,2,1,3,1,2,1,2,1,1,1,1,1,1,2,2,1,1,1,1,1,1,3,1,1,1,1,1,1,1,6,5,1,1,2,7,1,2,1,1,3,1,5,1,1,1,1,2,1,1,2,1,1,1,1,1,1,2,1,2,1,3,1,1,1,1,1,1,1,2,5,1,3,1,1,1,1,1,7,5,1,1,2,2,3,1,1,1,4,1,1,1,1,1,1,14,1,1,3,1,1,3,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,2,3,1,1,2,1,1,1,1,1,8,2,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,1,1,1,1,1,1,1,1,1,1,1,7,1,1,1,1,1,1,3,1,1,1,1,1,36,1,1,1,1,1,1,2,1
```

## Mathlink Usage
![MathlinkUsage](https://github.com/CarlosManuelRodr/DiscretePowerlawFitter/raw/main/img/MathlinkUsage.png)

## Building

### Requisites
- [gsl](https://www.gnu.org/software/gsl/) Gnu Scientific Library
- [gsl CBLAS](https://www.gnu.org/software/gsl/doc/html/cblas.html)
- [cmake Version >= 3.21](https://cmake.org/)
- A C++17 compiler.

### Instructions
On Unix-based system just use the `BuildCLI.sh` script. If you wish to build the Mathlink version, use the `BuildMathlink.sh` script. To build as a shared library use `BuildLib.sh`.


## Credits
- [thermal_funcs](https://github.com/andrewfowlie/thermal_funcs) For the implementation of the Hurwitz zeta function.
- [The Lean Mean C++ Option Parser](http://optionparser.sourceforge.net/) Is used for the CLI version.
- [thread-pool](https://github.com/bshoshany/thread-pool) Is used as safe parallelization method.
- [FindMathematica](https://github.com/sakra/FindMathematica) for the Mathlink build.