# `JDRPM` setup guide
1. Be sure to have Julia installed. You can download it from [her official site](https://julialang.org/downloads/).

2. Install on R the [`JuliaConnectoR`](https://github.com/stefan-m-lenz/JuliaConnectoR) library (see the installation section there on the github link for further instructions or if you have any problem). This is like installing any other library, so on R just do:
```R
install.packages("JuliaConnectoR")
```

3. Download and activate the package relative to `JDRPM`. You can do this step also staying inside R. To do it, create a section that looks like this:
```R
library(JuliaConnectoR) # load the interface library
juliaSetupOk() # check it returns TRUE

juliaEval("using Pkg") # load the Package manager on Julia

juliaEval("Pkg.activate(\"<path/to/where/you/stored/JDRPM>\")") # load the JDRPM pacakge

# the next line has to be executed only the first time, i.e. at the installation
# since it is required to download (and install, only once) all the depdendencies
juliaEval("Pkg.instantiate()") 
# but actually Julia Package manager is smart enough to not re-install everything,
# so you can even leave it there after the first installation
# go here for more details https://pkgdocs.julialang.org/v1/environments/

# about <path/to/where/you/stored/JDRPM>, in my case was for example this:
# juliaEval("Pkg.activate(\"../../JDRPM\")")
# so it just depends on where it is the current R file you are working on

juliaEval("Pkg.status()") # a check to have correctly loaded the package;
# it should print something like this:
# Project JDRPM v0.1.0
# Status `C:\Users\feder\Desktop\Uni magistrale\Tesi\src\JDRPM\Project.toml`
#   [6e4b80f9] BenchmarkTools v1.5.0
# ⌃ [31c24e10] Distributions v0.25.108
#   [efd6af41] ProfileCanvas v0.1.6
#   [92933f4c] ProgressMeter v1.10.2
#   [276daf66] SpecialFunctions v2.4.0
#   [2913bbd2] StatsBase v0.34.3
#   [a759f4b9] TimerOutputs v0.5.24
#   [ade2ca70] Dates
#   [37e2e46d] LinearAlgebra
#   [10745b16] Statistics v1.10.0
# Info Packages marked with ⌃ have new versions available and may be upgradable.

module = normalizePath("<path/to/where/you/stored/JDRPM>/src/JDRPM.jl") # locate the "main" file
module_JDRPM = juliaImport(juliaCall("include", module)) # load the "main" file
```
4. To actually fit the model and get the results, you have to call the fitting function which however returns a Julia object, so then you have to convert it back to an R object. That is, the fitting process implies actually two steps:
```R
#### 1) the fit, with ... being all the parameters
out = module_JDRPM$MCMC_fit(...) 

#### 2) the conversion, and reshape
rout = juliaGet(out)
names(rout)  = c("Si","gamma","alpha", "sigma2h", "muh", "eta1","beta","theta", "tau2", 
    "phi0", "phi1","lambda2","fitted","llike","lpml","waic")

# reshape some stuff to uniform it to DRPM output
rout$Si      = aperm(rout$Si,       c(2, 1, 3))
rout$gamma   = aperm(rout$gamma,    c(2, 1, 3))
rout$sigma2h = aperm(rout$sigma2h,  c(2, 1, 3))
rout$muh     = aperm(rout$muh,      c(2, 1, 3))
rout$fitted  = aperm(rout$fitted,   c(2, 1, 3))
rout$llike   = aperm(rout$llike,    c(2, 1, 3))
rout$alpha   = aperm(rout$alpha,    c(2, 1))
rout$theta   = aperm(rout$theta,    c(2, 1))
rout$tau2    = aperm(rout$tau2,     c(2, 1))
rout$eta1    = aperm(rout$eta1,     c(2, 1))
rout$phi0    = matrix(rout$phi0,    ncol = 1)
rout$phi1    = matrix(rout$phi1,    ncol = 1)
rout$lambda2 = matrix(rout$lambda2, ncol = 1)
```

5. `JDRPM` is finally ready-to-use. The `MCMC_fit` has already some quite self explanatory argument names. However here there is a more detailed list. Otherwise see [here](#tests-and-examples) for some examples, or scan trough the julia code itself, or read the modeling and implementation chapters of the Tesi.pdf document (when it will be finished).
```julia
function MCMC_fit(;
  # Y::Matrix{Float64},                   # n*T matrix, the observed values
  Y::Union{Matrix{Float64},Matrix{Union{Missing, Float64}}},   # n*T matrix, the observed values
  # Y::Matrix,                            # n*T matrix, the observed values
  sp_coords = missing,                  # n*2 matrix, the spatial coordinates
  Xlk_covariates = missing,             # n*p*T matrix, the covariates to include in the likelihood
  Xcl_covariates = missing,             # n*p*T matrix, the covariates to include in the clustering process

  M_dp::Float64,                        # Dirichlet mass parameter
  initial_partition = missing,          # Initial partition (if provided)

  starting_alpha::Float64,              # Starting value for alpha
  unit_specific_alpha::Bool,            # Unit-specific alpha values
  time_specific_alpha::Bool,            # Time-specific alpha values
  update_alpha::Bool,                   # Update alpha?
  
  include_eta1::Bool,                   # Include the autoregressive part of eta1?
  include_phi1::Bool,                   # Include the autoregressive part of phi1?
  update_eta1::Bool,                    # Update the autoregressive part of eta1?
  update_phi1::Bool,                    # Update the autoregressive part of phi1?

  sig2h_priors::Vector{Float64},        # Prior parameters for sig2h ∼ invGamma(a_sigma,b_sigma)
  eta1_priors::Vector{Float64},         # Prior parameters for eta1 ∼ Laplace(0,b) so it's the scale parameter b
                                        # plus the std dev for the Metropolis update trough N(μ=eta1_old,σ=mhsig_eta1)
  beta_priors = missing,                # Prior parameters for beta ∼ the mean vector and the s^2 param in fron of the Id matrix
  tau2_priors::Vector{Float64},         # Prior parameters for tau2 ∼ invGamma(a_tau, b_tau)
  phi0_priors::Vector{Float64},         # Prior parameters for phi0 ∼ N(μ=m0,σ^2=s0^2)
  phi1_priors::Float64,                 # Prior parameters for phi1 ∼ U(-1,1)
                                        # so we just need the std dev of the Metropolis update trough N(μ=phi1_old,σ=mhsig_phi1)
  lambda2_priors::Vector{Float64},      # Prior parameters for lambda2 ∼ invGamma(a_lambda, b_lambda)
  alpha_priors::AbstractArray{Float64}, # Prior parameters for alpha ∼ Beta(a_alpha, b_alpha)
                                        # but possibly that pair for each unit j, that's why the abstract array
  
  spatial_cohesion_idx = missing,       # cohesion choice
  sp_params = missing,                  # Parameters for spatial cohesion functions
  covariate_similarity_idx = missing,   # similarity choice
  cv_params = missing,                  # Parameters for covariates similarity functions
  cv_weight = 1.0,                      # factor to which scale the covariate similarity values

  draws::Real,                          # Number of MCMC draws
  burnin::Real,                         # Number of burn-in
  thin::Real,                           # Thinning interval

  beta_update_threshold = 0,            # if to update beta regressor only after some iterates
  logging = false,                      # Wheter to save execution infos to log file
  seed::Real,                           # Random seed for reproducibility
  simple_return = false,                # Return just the partition Si
  verbose = false,                      # if to print additional info
  perform_diagnostics = false,          # if to perform diagnostics on the fit
  skip_checks = false                   # if to skip initial checks
  )
```

# Tests and examples
Try to run the [`Tesi/src/JDRPM/test/JDRPM_small_example.Rmd`](https://github.com/federicomor/Tesi/blob/main/src/JDRPM/test/JDRPM_small_example.Rmd) file to see if everything works fine.   

Try instead [`Tesi/src/test/1 Assessing correctness and NA/assessing_correctness.Rmd`](<https://github.com/federicomor/Tesi/blob/main/src/test/1 Assessing correctness and NA/assessing_correctness.Rmd>) to see more complete examples of fits. In this file there are fits with only the target variable, with spatial information, with missing data, with covariates in the likelihood, in the prior, etc: all possible combinations of usage of JDRPM.