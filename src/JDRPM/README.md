- For further information you can contact me or also [Alessandro Carminati](https://github.com/AleCarminati).  
- For the complete explanation of this work read the Thesis document, available in `Tesi/upload/2024_12_Mor_Thesis.pdf`, i.e. [here](https://github.com/federicomor/Tesi/blob/main/upload/2024_12_Mor_Thesis.pdf).
- For just a short summary of this work read the slides, available in `Tesi/upload/2024_12_Mor_Slides.pdf`, i.e. [here](https://github.com/federicomor/Tesi/blob/main/upload/2024_12_Mor_Slides.pdf).

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
# it should print the list of packages dependencies, like Distributions, LinearAlgebra, etc.

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

5. `JDRPM` is finally ready-to-use. The `MCMC_fit` has already some quite self explanatory argument names. However here there is a more detailed list. Otherwise see [here](#tests-and-examples) for some examples, or scan trough the julia code itself, or read the modeling and implementation chapters of the Tesi.pdf document.
```julia
function MCMC_fit(;
  Y::Union{Matrix{Float64},Matrix{Union{Missing, Float64}}},   # n*T matrix, the observed values
                                        # this strange type combination to allow missing data entries
  sp_coords = missing,                  # n*2 matrix, the spatial coordinates
  Xlk_covariates = missing,             # n*p*T matrix, the covariates to include in the likelihood
  Xcl_covariates = missing,             # n*p*T matrix, the covariates to include in the clustering process

  M_dp::Float64,                        # Dirichlet mass parameter
  initial_partition = missing,          # Initial partition (if provided)

  starting_alpha::Float64,              # Starting value for alpha
  unit_specific_alpha::Bool,            # Employ a unit-specific alpha?
  time_specific_alpha::Bool,            # Employ a time-specific alpha?
  update_alpha::Bool,                   # Update alpha?
  
  include_eta1::Bool,                   # Include the autoregressive part of eta1?
  include_phi1::Bool,                   # Include the autoregressive part of phi1?
  update_eta1::Bool,                    # Update the autoregressive part of eta1?
  update_phi1::Bool,                    # Update the autoregressive part of phi1?

  sig2h_priors::Vector{Float64},        # Prior parameters for sig2h ∼ invGamma(a_sigma=...,b_sigma=...)
  eta1_priors::Vector{Float64},         # Prior parameters for eta1 ∼ Laplace(0,b=...) so it's the scale parameter b
                                        # plus the std dev for the Metropolis update trough N(μ=eta1_old,σ=...)
  beta_priors = missing,                # Prior parameters for beta ∼ the mean vector and the s^2 param in fron of the Id matrix,
                                        # more precisely, a vector v of length p+1, containing in indexes 1:p the mean vector μ,
                                        # while in last position the s^2 term, i.e. beta ∼ N(μ=v[1:p],Σ=v[end]*Id)  
  tau2_priors::Vector{Float64},         # Prior parameters for tau2 ∼ invGamma(a_tau=..., b_tau=...)
  phi0_priors::Vector{Float64},         # Prior parameters for phi0 ∼ N(μ=...,σ^2=...) so mean and variance
                                        # beware of the differences between variance and std dev, here in these arguments
  phi1_priors::Float64,                 # Prior parameters for phi1 ∼ U(-1,1)
                                        # so we just need the std dev of the Metropolis update trough N(μ=phi1_old,σ=...)
  lambda2_priors::Vector{Float64},      # Prior parameters for lambda2 ∼ invGamma(a_lambda=...,b_lambda=...)
  alpha_priors::AbstractArray{Float64}, # Prior parameters for alpha ∼ Beta(a_alpha=...,b_alpha=...)
                                        # but possibly that pair for each unit j, that's why the abstract array
  
  spatial_cohesion_idx = missing,       # cohesion choice
  sp_params = missing,                  # Parameters for spatial cohesion functions
  covariate_similarity_idx = missing,   # similarity choice
  cv_params = missing,                  # Parameters for covariates similarity functions
  cv_weight = 1.0,                      # factor to which scale the covariate similarity values

  beta_update_threshold = 0,            # from which iterate start to update the beta regressor
  # harmless insertion to maybe let the model focus first on the real relevant parameters, and then move to update beta
  # otherwise I thought that early "bad" samples for beta could damage the more significant clusters' parameters
  
  draws::Real,                          # Number of MCMC draws
  burnin::Real,                         # Number of burn-in
  thin::Real,                           # Thinning interval
  # these variables are Reals and not Ints since integer values on R (like 1000) are automatically casted into floats (1000.0) unless
  # we explicitly write as.int(value), on R, which is tedious, so I just left Real as type, should not be much performance-relevant

  logging = false,                      # Wheter to save execution infos to log file
  seed::Real,                           # Random seed for reproducibility
  simple_return = false,                # Return just the partition Si
  verbose = false,                      # if to print additional info
  perform_diagnostics = false,          # if to compute convergence diagnostics (e.g. Rhat, ess) on the sampled parameters
  skip_checks = false                   # if to skip initial checks
  )
```

# Spatial cohesions syntax
A brief recap on how to call the fitting function with respect to the choice of the cohesion function.

$$C_1(S_{h},\vec{s}_{h}^\star) = \begin{cases}
    \dfrac{M \cdot \Gamma(|S_h|)}{\Gamma(\alpha \mathcal{D}_h) \indicator{\mathcal{D}_h\geq 1} + \mathcal{D}_h \indicator{\mathcal{D}_h<1}} & \text{if } |S_h| > 1\\ 
   \hfil M & \text{if } |S_h| = 1
\end{cases}$$
```
spatial_cohesion_idx = 1, 
sp_params = alpha, # the alpha parameter
```

$$ C_2(S_{h},\vec{s}_{h}^\star) = M \cdot \Gamma(|S_h|) \cdot \prod_{i,j \in S_h} \indicator{\| \vec{s}_i - \vec{s}_j \| \leq a} $$
```
spatial_cohesion_idx = 2, 
sp_params = a, # the a parameter
```


$$ C_3(S_{h},\vec{s}_{h}^\star) = M \cdot \Gamma(|S_h|) \cdot \int \prod_{i \in S_h} q(\vec{s}_i|\vec{\csi}_h) q(\vec{\csi}_h) \, d\vec{\csi}_h $$


$$ C_4(S_{h},\vec{s}_{h}^\star) = M \cdot \Gamma(|S_h|) \cdot \int \prod_{i \in S_h} q(\vec{s}_i|\vec{\csi}_h) q(\vec{\csi}_h|\vec{s}_h^\star) \, d\vec{\csi}_h $$
```
spatial_cohesion_idx = 4, # or 4
sp_params = list(c(mu0,mu0),k0,v0,matrix(c(L0,0.0,0.0,L0),nrow=2)),
# the parameter set of mu0 (a vector), k0 and v0 (scalars), and the matrix L0
```

$$ C_5(S_{h},\vec{s}_{h}^\star) = M\cdot \Gamma(|S_h|)\cdot \expp{-\phi \sum_{i \in S_h} \| \vec{s}_i - \vec{\bar{s}}_h\| } $$
$$ C_6(S_{h},\vec{s}_{h}^\star) = M\cdot \Gamma(|S_h|)\cdot \expp{-\phi \log( \sum_{i \in S_h} \| \vec{s}_i - \vec{\bar{s}}_h\| )} $$ 
```
spatial_cohesion_idx = 5, # or 6 
sp_params = phi, # the phi parameter
```

# Covariates similarities syntax




# Tests and examples
Try to run the [`Tesi/src/JDRPM/test/JDRPM_small_example.Rmd`](https://github.com/federicomor/Tesi/blob/main/src/JDRPM/test/JDRPM_small_example.Rmd) file to see if everything works fine.   

Try instead [`Tesi/src/test/1 Assessing correctness and NA/assessing_correctness.Rmd`](<https://github.com/federicomor/Tesi/blob/main/src/test/1 Assessing correctness and NA/assessing_correctness.Rmd>) to see more complete examples of fits. In this file there are fits with only the target variable, with spatial information, with missing data, with covariates in the likelihood, in the prior, etc: all possible combinations of usage of JDRPM.
