struct params
	y::Matrix{Float64} # n*T
	sp_coords::Matrix{Float64} # n*2
	X_covariates::Array{Float64, 3} # n*p*T, p covariates for each unit and for all times
	M_dp::Int64 # Dirichlet mass parameter
	initial_parition # la vogliamo davvero?
	starting_alpha
	unit_specific_alpha
	time_specific_alpha
	alpha_0 # maybe better to call it "update_alpha"?
	eta1_0 # maybe better to call it "update_alpha"?
	phi1_0 # maybe better to call it "update_alpha"?
	model_priors::Vector{Float64}
	alpha_priors::Vector{Float64}
	spatial_cohesion # maybe we can use a Dict here to maps integers/strings to functions
	sp_params::Vector{Float64} # spatial cohesion params
	covariate_similarity # as above
	mh::Vector{Float64} # for Metropolis updates
	draws::Int64
	burnin::Int64
	thin::Int64
	verbose::Bool
end