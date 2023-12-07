
<!-- README.md is generated from README.Rmd. Please edit that file -->

# malariaEquilibriumVivax

[![Travis build
status](https://travis-ci.org/mrc-ide/malariaEquilibriumVivax.svg?branch=master)](https://travis-ci.org/mrc-ide/malariaEquilibriumVivax)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/malariaEquilibriumVivax?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/malariaEquilibriumVivax)

Often we are interested in the state of a given malaria transmission
model at equilibrium. However, some models (e.g. the Griffin et al. 2014
model, the White et al. 2018 model) are quite complex, and can result in
equilibrium solutions that are fairly in-depth. In these situations it
is useful to have a “canonical” equilibrium solution that is tried and
tested, and can be used reliably by multiple users. This package
complements the malariaEquilibrium package to store canonical
equilibrium solutions for Plasmodium vivax malaria models.

## How to interact with this repos

If you have a new model with an equilibrium solution that you would like
to host here, then please go ahead and do so *in a new branch*. The
master branch, on the other hand, should always correspond to the
current most active model used by the malaria group (currently the White
et al. 2018 model). This should generally be left alone unless you know
what you’re doing.

Any changes that you propose to branches that you did not create
yourself (e.g. the master branch) should be done through pull requests,
rather than editing directly. This gives others the opportunity to
review your changes before they are implemented in code.

NB. These stipulations are in place because there are other R packages
that test against the equilibrium solutions from this package, and so
any changes here may cause other packages to fail if they are not know
about beforehand.

## Installation

You can install from the master branch of malariaEquilibriumVivax with:

``` r
# install.packages("devtools")
devtools::install_github("mrc-ide/malariaEquilibriumVivax")
```

Or from a specific branch with:

``` r
# install.packages("devtools")
devtools::install_github("mrc-ide/malariaEquilibriumVivax", ref = "mybranch")
```

Finally, don’t forget to load the package with:

``` r
library(malariaEquilibriumVivax)
library(malariasimulation)
```

## Example

### Working with parameters

It’s possibly easiest to begin with the default vivax parameter set as
found in malariasimulation.

``` r
vivax_params <- get_parameters(parasite = "vivax")
```

We can modify parameter values exactly as we would with a named list:

``` r
vivax_params$eta <- 0.1
head(vivax_params)
#> $dd
#> [1] 5
#> 
#> $dt
#> [1] 1
#> 
#> $da
#> [1] 10
#> 
#> $dpcr_max
#> [1] 70
#> 
#> $dpcr_min
#> [1] 10
#> 
#> $kpcr
#> [1] 4.602
```

These do not have parameter names that are consistent with Michael
White’s model and must be translated for use in the equilibrium
solution.

``` r
vivax_eq_params <- malariasimulation:::translate_vivax_parameters(vivax_params)
```

### Getting the equilibrium solution

Now we can calculate the equilibrium solution using the
`vivax_equilibrium()` function. This assumes a fixed EIR and treatment
rate (`ft`), and takes model parameters as input. It calculates the
equilibrium solution for a defined age range:

``` r
eq <- vivax_equilibrium(EIR = 10, ft = 0.2, p = vivax_eq_params, age = 0:10)
```

We have two equilibrium solutions: a full solution that loops over
hypnozoite batch numbers until reaching a limit (the default limit is
set to 10 hypnozoite batches), and a simplified, approximate, version
that calculates the average for each age and biting heterogeneity group.

``` r
eq <- vivax_equilibrium(EIR = 10, ft = 0.2, p = vivax_eq_params, age = 0:10, v_eq = "simplified")
```

The output elements are “states” and “FOIM”. The former is a matrix with
age brackets in rows, and model states or other downstream calculations
in columns. For example, the column “S” gives the equilibrium proportion
of humans in the susceptible state.

``` r
head(eq$states[[1]])
#>   age       prop            S            U            A            D
#> 1   0 0.04347126 8.761186e-05 0.0001821670 0.0001778841 2.704346e-05
#> 2   1 0.04158151 8.125013e-05 0.0001506895 0.0002141042 1.417579e-05
#> 3   2 0.03977391 6.917851e-05 0.0001696545 0.0001952529 8.791530e-06
#> 4   3 0.03804489 6.131761e-05 0.0001909988 0.0001654373 6.787005e-06
#> 5   4 0.03639103 5.552787e-05 0.0002068250 0.0001384096 5.742089e-06
#> 6   5 0.03480906 5.073966e-05 0.0002167239 0.0001165638 5.050992e-06
#>              T            P       ID          IDM       ICA          ICM
#> 1 1.367281e-06 1.330026e-05 10.57063 3.998968e-01  18.34678 2.714350e+00
#> 2 7.022980e-07 7.178265e-06 22.33676 1.235799e-05  43.05848 8.388141e-05
#> 3 4.368955e-07 4.436999e-06 34.00066 3.818982e-10  70.75929 2.592182e-09
#> 4 3.384367e-07 3.407793e-06 45.14801 1.180178e-14  99.90119 8.010607e-14
#> 5 2.866831e-07 2.877534e-06 55.62886 3.647096e-19 129.66663 2.475513e-18
#> 6 2.523018e-07 2.529270e-06 65.39694 1.127060e-23 159.58217 7.650064e-23
#>          HH EIR       inf       E_M       I_M  inf_rvoir      mv0       psi
#> 1  3.093588  10 0.3404829 0.4828513 0.2377735 0.01254184 258.7219 0.3090514
#> 2  6.138467  10 0.3404829 0.4828513 0.2377735 0.02794163 258.7219 0.4529588
#> 3  8.988822  10 0.3404829 0.4828513 0.2377735 0.04860308 258.7219 0.5799565
#> 4 11.588919  10 0.3404829 0.4828513 0.2377735 0.07346903 258.7219 0.6920317
#> 5 13.927406  10 0.3404829 0.4828513 0.2377735 0.10152487 258.7219 0.7909376
#> 6 16.013886  10 0.3404829 0.4828513 0.2377735 0.13200393 258.7219 0.8782218
#>     phi_clin phi_patent       dPCR
#> 1 0.19089373  0.8042138 0.03026526
#> 2 0.05938605  0.5571232 0.08782813
#> 3 0.03199827  0.3375802 0.09799606
#> 4 0.02450078  0.2112293 0.09944715
#> 5 0.02162766  0.1422464 0.09978762
#> 6 0.02027758  0.1027153 0.09989899
```

The “FOIM” element gives the onward force of infection from humans to
mosquitoes, which is a single number:

``` r
eq$FOIM
#> [1] 0.1423578
```
