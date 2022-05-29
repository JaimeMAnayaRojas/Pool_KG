# Statistical analyses and Integral Projection Models
In this repository, we provide the full code to replicate our analyses. We parametarized our IPM using experimental data.

## Part I: Statistical model
The statistical models were run in `stan` via `R` using the `rstan` packages. This is done by runing the `R/MainScript.R` script.

## Part II: Integral Projection Models (IPM) using Julia
We used IPM to estimate the demographic effects of species removal in Trinidadian streams. If you want to run this models using the new and amazing programing lenguage, [`Julia`](https://julialang.org/), you can run the `Julia/MainScript.jl` file. For a detail explanation of the construction of the IPM models you can read the [🎈 Notebook_full_analyses.jl — Pluto.jl.pdf ](https://github.com/JaimeMAnayaRojas/KG/blob/main/%F0%9F%8E%88%20Notebook_full_analyses.jl%20%E2%80%94%20Pluto.jl.pdf)


## Part III: Integral Projection Models (IPM) using `R`
At first, we run the IPM models in R using the `R/MainScript.R`, but it takes approximatly three days to run. In Julia, it takes less than five hours.
