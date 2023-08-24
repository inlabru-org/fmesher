## Changes from 0.1.1 to 0.1.2 (requested submission)

* Various documentation improvements, in particular for INLA compatibility
* Modify mesh refinement tests to directly check the refinement criteria
  instead of the specific mesh result, to account for differences in
  floating point behaviour on M1 processors.
* Modify tests of non-mesh-generation features to use precomputed meshes
  or meshes with stable properties
* Protect against invalid `tv` inputs
* Revert from `\text{}` to `\textrm{}`, as AMS extensions are only supported
  from R 4.2.2 (https://www.stats.bris.ac.uk/R/doc/manuals/r-devel/R-exts.pdf
  2023-08-24, page 90), and CRAN oldrel for macOS is 4.2.0, not 4.2.3

## R CMD check results for 0.1.2

0 errors | 0 warnings | 0 notes

## CRAN results for 0.1.1

### C++ library build size note:

Version: 0.1.1
Check: installed package size
Result: NOTE
     installed size is 8.3Mb
     sub-directories of 1Mb or more:
     libs 6.2Mb
Flavors: r-release-macos-x86_64, r-oldrel-macos-x86_64

### Package test suite error, identified and fixed in 0.1.2:

Version: 0.1.1
Check: tests
Result: ERROR
     Running ‘testthat.R’ [11s/12s]
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
     `expected`: 7.781
     ── Failure ('test-fmesher_integration.R:184:3'): Polygon integration with holes ──
     sum(ipA3$weight) (`actual`) not equal to 8.004363 (`expected`).
    
     `actual`: 8.002
     `expected`: 8.004
     ── Failure ('test-fmesher_integration.R:185:3'): Polygon integration with holes ──
     sum(ipA4$weight) (`actual`) not equal to 8.004363 (`expected`).
    
     `actual`: 8.002
     `expected`: 8.004
    
     [ FAIL 4 | WARN 0 | SKIP 9 | PASS 180 ]
     Error: Test failures
     Execution halted
Flavor: r-release-macos-x86_64

### Manual building error on oldrel-macos, fixed in 0.1.2

AMS support (\text) was apparently introduced in 4.2.2, not 4.2.0

Version: 0.1.1
Check: PDF version of manual
Result: WARN
    LaTeX errors when creating PDF version.
    This typically indicates Rd problems.
    LaTeX errors found:
    ! Undefined control sequence.
    <argument> \log (\sum _{i; \text
...
Flavor: r-oldrel-macos-x86_64


### M1mac package test suite error, fixed in 0.1.2

The problem was replicated with Linux/M1/clang-16 and the tests rewritten
to test the appropriate invariants.

> test_check("fmesher")
  Starting 2 test processes
  [ FAIL 1 | WARN 0 | SKIP 9 | PASS 183 ]
