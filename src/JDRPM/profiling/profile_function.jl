using ProfileCanvas
using Profile

# run in the REPL, not as a new process
include("../MCMC_fit.jl")
N = 80
T = 12
y = rand(N,T)
sp = rand(N,2)

# params
m0_phi0 = 0.
s20_phi0 = 1.
a_sigma  = 2.; b_sigma  = 2.
a_tau    = 2.; b_tau    = 2.
a_lambda = 2.; b_lambda = 2.
eta1_scale = 0.9
sig_mh_eta1 = 0.1
sig_mh_phi1 = 0.1
update_eta1 = true
update_phi1 = true
a_alpha = 1.; b_alpha = 1.
time_specific_alpha = true
# now space
spatial_cohesion_idx = 3.
mu0 = 0.
k0 = 1.
v0 = 5.
L0 = 1.

niter = 1000.
burnin = 0.
thin = 1.
seed = 123.0

################### with ProfileCanvas (better)

# ProfileCanvas.@profview_allocs MCMC_fit(
ProfileCanvas.@profview MCMC_fit(
    Y=y,              
    sp_coords = sp,
    M_dp = 1.0,                     
    initial_partition = missing,
    Xlk_covariates = missing,
    Xcl_covariates = missing,
    starting_alpha = 0.5,         
    unit_specific_alpha = false,       
    time_specific_alpha = time_specific_alpha,     

    update_alpha = true,             
    include_eta1 = true,                    
    include_phi1 = true,
    update_eta1 = update_eta1,                    
    update_phi1 = update_phi1,

    sig2h_priors = [a_sigma,b_sigma],
    eta1_priors = [eta1_scale,sig_mh_eta1^2],
    # beta_priors = c(rep(1,p),2),
    beta_priors = missing,
    tau2_priors = [a_tau,b_tau],
    phi0_priors = [m0_phi0,s20_phi0],
    phi1_priors = sig_mh_phi1^2,
    lambda2_priors = [a_lambda,b_lambda],
    alpha_priors = [a_alpha,b_alpha],
    
    spatial_cohesion_idx = spatial_cohesion_idx,
    sp_params = [[mu0,mu0],k0,v0,[L0 0.0; 0.0 L0]],
    
    # covariate_similarity_idx = NA,  
    draws = niter,                    
    burnin = burnin,                   
    thin = thin,                     
    logging = true,
    seed = seed
)