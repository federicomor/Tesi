using Distributions
using LinearAlgebra
using Random

log_file = open("output.txt", "w")
function debug(str...)
	println(log_file,str...)
end

include("utils.jl")

function MCMC_fit(;
	Y::Matrix{Float64},               # n*T, observed values
	sp_coords = missing,              # n*2, spatial coordinates
	X_covariates = missing,           # n*p*T, covariates for each unit and all times
	M_dp::Float64,                    # Dirichlet mass parameter
	initial_partition = missing,      # Initial partition (if provided)
	starting_alpha::Float64,          # Starting value for alpha
	unit_specific_alpha::Bool,        # Unit-specific alpha values
	time_specific_alpha::Bool,        # Time-specific alpha values
	update_alpha::Bool,               # Update alpha?
	include_eta1::Bool,               # Include the autoregressive part of eta1?
	include_phi1::Bool,               # Include the autoregressive part of phi1?
	sigma2h_priors::Vector{Float64},  # Prior parameters for sigma2h
	eta_priors::Vector{Float64},      # Prior parameters for eta
	beta_priors::Vector{Float64},     # Prior parameters for beta
	tau2_priors::Vector{Float64},     # Prior parameters for tau2
	phi0_priors::Vector{Float64},     # Prior parameters for phi0
	phi1_priors::Vector{Float64},     # Prior parameters for phi1
	lambda2_priors::Vector{Float64},  # Prior parameters for lambda2
	alpha_priors::Matrix{Float64},    # Prior parameters for alpha
	spatial_cohesion = missing,       # cohesion choice
	sp_params = missing,              # Parameters for spatial cohesion functions
	covariate_similarity = missing,   # similarity choice
	mh::Vector{Float64},              # Metropolis-Hastings 
	draws::Float64,                   # Number of MCMC draws
	burnin::Float64,                  # Number of burn-in sa
	thin::Float64,                    # Thinning interval
	verbose::Bool,                    # Verbosity flag	
	seed::Float64,                    # Random seed for reproducibility
	io=log_file
	)

	############# define auxiliary variables #############
	sPPM = !ismissing(sp_coords)
	xPPM = !ismissing(X_covariates)
	n, T = size(Y)
	T_star = T+1
	p = xPPM ? size(X_covariates)[2] : 0
	nout = Int64((draws - burnin)/(thin))
	if sPPM
		sp1 = vec(sp_coords[:,1])
		sp2 = vec(sp_coords[:,2])
	end

	############# send feedback #############
	println("fitting $draws iterates ($nout actual iterates)")
	println("on n=$n subjects\nfor T=$T time instants")
	println("with space? $sPPM\nwith covariates? $xPPM")
	println()

	############# allocate output variables #############
	Si_out = zeros(Int64,nout,n,T)
	gamma_out = zeros(nout,n,T)
	alpha_out = zeros(nout,T) 
	sigma2h_out = zeros(nout,n,T) 
	muh_out = zeros(nout,n,T) 
	eta1_out = zeros(nout,n) 
	if xPPM
		beta_out = zeros(nout,p,T) # or betareg as name?
	end
	theta_out = zeros(nout,T) 
	tau2_out = zeros(nout,T) 
	phi0_out = zeros(nout) 
	phi1_out = zeros(nout) 
	lambda2_out = zeros(nout) 
	fitted = zeros(nout,T,n)
	llike = zeros(nout,T,n)
	lpml = 0
	waic = 0

	############# allocate working variables #############
	Si_iter = ones(Int64,n,T_star)
	Si_iter[:,end] .= 0
	gamma_iter = zeros(n,T_star)
	alpha_iter = ones(n,T_star)*starting_alpha
	sigma2h = ones(n,T_star)
	muh = zeros(n,T_star)
	eta1_iter = zeros(n)
	if xPPM
		beta_iter = rand(MvNormal(beta_priors[1:end-1],beta_priors[end]*Matrix(I, p, p))) # or betareg as name?
	end
	theta_iter = rand(Normal(),T_star)
	tau2_iter = rand(InverseGamma(tau2_priors...),T_star)
	phi0_iter = (T==2) ? mean(Normal(phi0_priors...)) : rand(Normal(phi0_priors...))
	phi1_iter = rand(Uniform(-1,1))*include_phi1
	lambda2_iter = (T==2) ? mean(InverseGamma(lambda2_priors...)) : rand(InverseGamma(lambda2_priors...))
	
	nh = zeros(n,T_star)
	nh[1,:] .= n
	nclus_iter = ones(T_star)


	############# auxiliary working variables #############


	############# testing on random values #############
	Si_iter = rand(collect(1:5),n,T_star)
	gamma_iter = rand((0,1),n,T_star)
	

	
	println("setting seed $seed")
	Random.seed!(round(Int64,seed))
	
	############# start MCMC algorithm #############
	for i in 1:draws
		print("iteration $i\r")
		for t in 1:T
			
			debug("##### iteration $i - time $t #####")
			
			############# update gamma #############
			for j in 1:n
				if t==1 
					gamma_iter[:,t] .= 0
				else
					# we want to find ρ_t^{R_t(-j)} ...
					indexes = findall(jj -> jj != j && gamma_iter[jj, t] == 1, 1:n)
					Si_red = Si_iter[indexes, t]
					Si_red1 = copy(Si_red)
					push!(Si_red1,Si_iter[j,t]) # ... and ρ_t^R_t(+j)}

					if sPPM
						sp1_red = sp1[indexes]
						sp2_red = sp2[indexes]
					end

					n_red = length(Si_red)
					n_red1 = length(Si_red1)
					relabel!(Si_red,n_red)
					relabel!(Si_red1,n_red1)
					# @info "Si_red = $Si_red", "Si_red1 = $Si_red1" _file=""
					debug("Si_red = $Si_red", "\nSi_red1 = $Si_red1")
					nclus_red = isempty(Si_red) ? 0 : maximum(Si_red)
					nclus_red1 = maximum(Si_red1)

					j_label = Si_red1[end]

					nh_red = zeros(Int64,n); nh_red1 = zeros(Int64,n)
					for jj in 1:n_red
						nh_red[Si_red[jj]] += 1
						nh_red1[Si_red1[jj]] += 1
					end
					nh_red1[Si_red1[end]] += 1 # account for the last unit j


					lCo = 0.0; lCn = 0.0

					# unit j can enter an existing cluster
					for k in 1:nclus_red
						if sPPM
							sp_idxs = findall(jj -> Si_red[jj] == k, 1:n_red)
							s1o = sp1[sp_idxs]
							s2o = sp2[sp_idxs]
							s1n = push!(copy(sp1),sp1[j])
							s2n = push!(copy(sp2),sp2[j])

							# lCo = spatial_cohesion(s1o, s2o,)
							# lCn = spatial_cohesion(s1n, s2n,)
						end
					end

				end
			end

			############# update rho #############
			
			############# update muh #############
			
			############# update sigma2h #############
			
			############# update theta #############
			
			############# update tau2 #############

		end
		
		############# update eta1 #############
		
		############# update alpha #############
		
		############# update phi0 #############
		
		############# update phi1 #############
		
		############# update lambda2 #############
	
		############# save MCMC iterates #############

	end
	println("\ndone!")
	
close(io)
# global_logger(ConsoleLogger())	
end

