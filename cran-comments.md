## Test Environments

* CentOS Linux 7, R-3.6.1
* Windows 10, R-3.6.2

## devtools check results

### devtools::check()

─  building ‘lipidomeR_0.1.2.tar.gz’
   Warning: invalid uid value replaced by that for user 'nobody'

✓  checking examples (9.8s)
   ** found \donttest examples: check also with --run-donttest
   
── R CMD check results ─────────────────────────────────────────────────────────
─ lipidomeR 0.1.2 ────
Duration: 2m 51.8s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

### devtools::check( run_dont_test = TRUE )

─  building ‘lipidomeR_0.1.2.tar.gz’
   Warning: invalid uid value replaced by that for user 'nobody'
   
✓  checking examples (39.5s)
   Examples with CPU or elapsed time > 5s
                   user system elapsed
   liverlipidome  7.204  0.084  19.010
   humanlipidome  1.903  0.004   6.087
   cancerlipidome 1.597  0.089   5.203

── R CMD check results ─────────────────────────────────────────────────────────
─ lipidomeR 0.1.2 ────
Duration: 2m 2.8s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

## R CMD check results

### R CMD check --as-cran ../lipidomeR_0.1.2.tar.gz

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Tommi Suvitaival <TSUV0001@RegionH.DK>’

New submission

* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

** found \donttest examples: check also with --run-donttest

Status: 2 NOTEs

### R CMD check --as-cran --run-donttest ../lipidomeR_0.1.2.tar.gz

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Tommi Suvitaival <TSUV0001@RegionH.DK>’

New submission

* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

* checking examples ... NOTE
Examples with CPU or elapsed time > 5s
               user system elapsed
liverlipidome 7.101  0.082   8.296
               
Status: 3 NOTEs
