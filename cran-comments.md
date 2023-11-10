This is a minor bugfix release.

## Changes from 0.1.3 to 0.1.4

* C++ code refactor to work around `std::get<variant>()` lack of support
  for MacOS < 10.14.
* Further C++ refactor towards std::unique_ptr use.

## R CMD check results for 0.1.4

> checking installed package size ... NOTE
    installed size is 11.7Mb
    sub-directories of 1Mb or more:
      libs   9.8Mb
      
This is the size of the C++ library after build&installation.

## CRAN results for 0.1.3

Version: 0.1.3
Check: installed package size
Result: NOTE
     installed size is 8.3Mb
     sub-directories of 1Mb or more:
     libs 6.4Mb
Flavors: r-release-macos-arm64, r-release-macos-x86_64, r-oldrel-macos-arm64

Version: 0.1.3
Check: whether package can be installed
Result: ERROR
    Installation failed.
Flavor: r-oldrel-macos-x86_64
    
# revdepcheck results

We checked 6 reverse dependencies (6 from CRAN), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
