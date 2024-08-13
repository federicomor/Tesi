n = 5
T = 3
p = 2
Y = rand(n,T)
sp_coords = rand(n,2)

include("MCMC_fit.jl")

seed_choice = rand(Uniform(0,1000))
test_draws  = 4.0
test_burnin = 2.0
test_thin   = 1.0

MCMC_fit(
	Y=Y,              
	sp_coords = sp_coords,             
	X_covariates = missing,
	M_dp = 2.0,                     
	initial_partition = missing,     
	starting_alpha = 0.5,         
	unit_specific_alpha = false,       
	time_specific_alpha = true,       
	update_alpha = true,                   
	include_eta1 = false,                    
	include_phi1 = false,
	sigma2h_priors = rand(2),
	eta1_priors = rand(2),
	beta_priors = repeat([rand()],p),
	tau2_priors = rand(2),
	phi0_priors = rand(2),
	phi1_priors = rand(2),
	lambda2_priors = rand(2),
	alpha_priors = rand(2),

	spatial_cohesion_idx = 1,
	sp_params = 0.5,
	# sp_params = c(0.5,1),
	# spatial_cohesion_idx = 3,
	# sp_params = list(c(1,2),1,2,matrix(c(1,2,2,4),nrow=2)),
	
	# covariate_similarity_idx = NA,  
	mh = [0.1,0.1,0.7, 0.1, 0.1],             
	draws = test_draws,                    
	burnin = test_burnin,                   
	thin = test_thin,                     
	verbose = false,
	seed = seed_choice
)

