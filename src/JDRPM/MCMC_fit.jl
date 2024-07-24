struct Params{T<:Real}                  # T to allow to choose Float32 or Float64
	# Data fields
	Y::Matrix{T}                        # n*T, observed values
	sp_coords::Matrix{T}                # n*2, spatial coordinates
	X_covariates::Array{T, 3}           # n*p*T, covariates for each unit and all times

	# # Model parameters
	# M_dp::Int64                         # Dirichlet mass parameter
	# initial_partition::Vector{Int64}    # Initial partition (if provided)
	# starting_alpha::T                   # Starting value for alpha
	# unit_specific_alpha::Vector{T}      # Unit-specific alpha values
	# time_specific_alpha::Vector{T}      # Time-specific alpha values
	# alpha_0::Bool                       # Do we want to update alpha?
	# eta1_0::Bool                        # Do we want to update eta1?
	# phi1_0::Bool                        # Do we want to update phi1?
	# model_priors::Vector{T}             # Prior distributions for the model parameters
	# alpha_priors::Vector{T}             # Prior distributions for alpha parameters

	# # Function mappings for spatial and covariate cohesion
	# spatial_cohesion::Dict{String, Function}
	# sp_params::Vector{T}                # Parameters for spatial cohesion functions
	# covariate_similarity::Dict{String, Function}

	# # MCMC control parameters
	# mh::Vector{T}                       # Metropolis-Hastings parameters
	# draws::Int64                        # Number of MCMC draws
	# burnin::Int64                       # Number of burn-in samples
	# thin::Int64                         # Thinning interval
	# verbose::Bool                       # Verbosity flag
end


import Base: show
function show(io::IO, params::Params)
    println(io, "Params{", typeof(params.Y[1,1]), "} struct:")
    println(io, "### Data fields:")
    println(io, "    Y: $(size(params.Y)) Matrix")
    println(io, "    sp_coords: $(size(params.sp_coords)) Matrix")
    println(io, "    X_covariates: $(size(params.X_covariates)) Array")
    # println(io, "### Model parameters:")
    # println(io, "    M_dp: $(params.M_dp)")
    # println(io, "    initial_partition: Vector of length $(length(params.initial_partition))")
    # println(io, "    starting_alpha: ", params.starting_alpha)
    # println(io, "    unit_specific_alpha: Vector of length $(length(params.unit_specific_alpha))")
    # println(io, "    time_specific_alpha: Vector of length $(length(params.time_specific_alpha))")
    # println(io, "    alpha_0: $(params.alpha_0)")
    # println(io, "    eta1_0: $(params.eta1_0)")
    # println(io, "    phi1_0: $(params.phi1_0)")
    # println(io, "    model_priors: Vector of length $(length(params.model_priors))")
    # println(io, "    alpha_priors: Vector of length $(length(params.alpha_priors))")
    # println(io, "### Function mappings:")
    # println(io, "    spatial_cohesion: Dict with $(length(params.spatial_cohesion)) entries")
    # println(io, "    sp_params: Vector of length $(length(params.sp_params))")
    # println(io, "    covariate_similarity: Dict with $(length(params.covariate_similarity)) entries")
    # println(io, "### MCMC control parameters:")
    # println(io, "    mh: Vector of length $(length(params.mh))")
    # println(io, "    draws: $(params.draws)")
    # println(io, "    burnin: $(params.burnin)")
    # println(io, "    thin: $(params.thin)")
    # println(io, "    verbose: $(params.verbose)")
end

# Example usage
# params = Params{Float64}(
#     randn(100, 50),
#     rand(100, 2),
#     randn(100, 3, 50),
#     1,
#     collect(1:100),
#     0.5,
#     rand(100),
#     rand(50),
#     true,
#     true,
#     false,
#     [1.0, 2.0, 3.0],
#     [1.0, 2.0],
#     Dict("gaussian" => (x, y) -> exp(-norm(x - y)^2)),
#     [1.0, 0.5],
#     Dict("linear" => (x, y) -> dot(x, y)),
#     [0.1, 0.2, 0.3],
#     1000,
#     200,
#     2,
#     true
# )
# println(params)


# function MCMC_fit(params::Params)
# 	println(size(params.Y))
# 	println(size(params.sp_coords))
# 	println(size(params.X_covariates))
# 	# for i in 1:params.draws
# 	# 	for t in 1:params.T
# 	# 	end
# 	# end
# end

function MCMC_fit(Y,sp_coords,X_covariates)
	println(size(Y))
	println(size(sp_coords))
	println(size(X_covariates))
	# for i in 1:params.draws
	# 	for t in 1:params.T
	# 	end
	# end
end