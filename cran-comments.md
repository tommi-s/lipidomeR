## Test Environments

* CentOS Linux 7, R-3.6.1
* Windows 10, R-3.6.2

## devtools check results

devtools::check()

─  building ‘lipidomeR_0.1.1.tar.gz’
   Warning: invalid uid value replaced by that for user 'nobody'
   Warning: invalid gid value replaced by that for user 'nobody'

✓  checking examples (10.4s)
   ** found \donttest examples: check also with --run-donttest

devtools::check( run_dont_test = TRUE )

─  building ‘lipidomeR_0.1.1.tar.gz’
   Warning: invalid uid value replaced by that for user 'nobody'
   Warning: invalid gid value replaced by that for user 'nobody'
   
✓  checking examples (19.2s)
   Examples with CPU or elapsed time > 5s
                  user system elapsed
   liverlipidome 7.802  0.064  10.123

Duration: 1m 40.6s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## R CMD check results

R CMD check --as-cran ../lipidomeR_0.1.1.tar.gz

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Tommi Suvitaival <TSUV0001@RegionH.DK>’

New submission

* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

** found \donttest examples: check also with --run-donttest

R CMD check --as-cran --run-donttest ../lipidomeR_0.1.1.tar.gz

* checking CRAN incoming feasibility ... NOTE

New submission
Maintainer: ‘Tommi Suvitaival <TSUV0001@RegionH.DK>’

* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

* checking examples ... NOTE
Examples with CPU or elapsed time > 5s
               user system elapsed
liverlipidome 7.801  0.074   8.966
