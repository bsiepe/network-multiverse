# Network Multiverse
This repository contains code and data for the preprint Network Multiverse (Siepe, Schumacher & Heck, 2023). (ADD OSF LINK) Start the .Rproj var-compare.Rproj prior to running the scripts. We used the `renv` package to create a reproducible environment; more information is contained in the corresponding folder.  

Full results of the multiverse analyses are too large and can be requested from the corresponding author. However, we provide all summaries of results necessary to reproduce our analyses in the `output/` folder. 

## Structure

### `data/`
Contains the raw data to reproduce our empirical analyses

### `output/`
Contains the summarized results of both the multiverse analyses as well as the supplementary simulation study

### `scripts/`
Contains all scripts needed to reproduce the results. The main analysis script is **gimme-multiverse.Rmd**. 
**aux_funs.R** contains all auxiliary functions used. 



## Package
The `mv-gimme` fork of the `gimme` package used for the multiverse analyses can be found in another [GitHub repository](https://github.com/bsiepe/mv-gimme).
Please note that the package is not actively maintained and may not be up to date with the original `gimme` package. I am happy to help if you would like to use the package for your own multiverse analysis or simulation study.



## Shiny App
There are multiple ways to run this shiny app. The easiest way is to use the version that is hosted online at (INSERT LINK).

Another possibility is to use the `shiny` package in R:

```
library(shiny)

runGitHub(repo = "network-multiverse",
          username = "bsiepe",
          ref = "main",
          subdir = "shiny-app")

```

Alternatively, you can clone the GitHub repository, navigate to the file `shiny-app/app.R` and use:

```
shiny::runApp()
```

If you happen to have any problems with the code, or any suggestions for further improvement, feel free to contact me via the email listed in my GitHub account. 