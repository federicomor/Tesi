using ProfileCanvas
using Profile
using Cthulhu

# run in the REPL, not as a new process
include("../MCMC_fit.jl")
N = 40
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

niter = 10000.
burnin = 9000.
thin = 5.
seed = 123.0

################### with ProfileCanvas (better)

# ProfileCanvas.@profview_allocs MCMC_fit(
# ProfileCanvas.@profview MCMC_fit(
# @descend MCMC_fit(
out = MCMC_fit(
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
);

using MCMCChains

out[1]
out[2]
out[3]
out[4]
out[5]
out[6]

alpha = out[3][:,801:end]
alpha_res = reshape(alpha,200,12,1)
alpha_1 = alpha_res[:,1,1]
alpha_2 = alpha_res[:,2,1]
alpha_3 = alpha_res[:,3,1]
vals = [alpha_1  alpha_2 alpha_3]
vals_res = reshape(vals,200,3,1)

chn = Chains(vals_res, [:alpha_1,:alpha_2,:alpha_3])
traceplot(chn)