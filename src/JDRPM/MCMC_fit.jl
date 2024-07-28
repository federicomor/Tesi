using Distributions
using LinearAlgebra
using Random
using Logging

############# debug helpers #############
log_file = open("log.txt", "w+")
function debug(str)
	print(log_file,str * (str[end]=='\n' ? "" : "\n"))
end
macro showd(exs...)
	args = []
	for ex in exs
		push!(args, sprint(Base.show_unquoted, ex), " = ", esc(ex), '\n')
	end
	return :(string($(args...)))
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
	spatial_cohesion_idx = missing,       # cohesion choice
	sp_params = missing,              # Parameters for spatial cohesion functions
	covariate_similarity_idx = missing,   # similarity choice
	mh::Vector{Float64},              # Metropolis-Hastings 
	draws::Float64,                   # Number of MCMC draws
	burnin::Float64,                  # Number of burn-in sa
	thin::Float64,                    # Thinning interval
	verbose::Bool,                    # Verbosity flag	
	seed::Float64,                    # Random seed for reproducibility
	io=log_file
	)
	println("setting seed $seed")
	Random.seed!(round(Int64,seed))

	############# check some stuff #############
	#=
	function cohesion1(s1::Vector{Float64}, s2::Vector{Float64}, alpha::Real, lg::Bool, M::Real=1.0)
	function cohesion2(s1::Vector{Float64}, s2::Vector{Float64}, a::Real, lg::Bool, M::Real=1.0)
	function cohesion3_4(s1::Vector{Float64}, s2::Vector{Float64}, mu_0::Vector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}, Cohesion::Int, lg::Bool, M::Real=1.0)
	function cohesion5(s1::Vector{Float64}, s2::Vector{Float64}, phi::Real, lg::Bool, M::Real=1.0)
	function cohesion6(s1::Vector{Float64}, s2::Vector{Float64}, phi::Real, lg::Bool, M::Real=1.0)
	sp_params should be 
		- alpha for cohesion 1
		- a for cohesion 2
		- mu_0 (2x1 vec), k0, v0, Psi (2x2 matrix) for cohesion 3_4
		- phi for cohesion 5
		- phi for cohesion 6
	=#
	if spatial_cohesion == 1 && !(sp_params isa Real) 
		@error "wrong params for spatial cohesion 1.\nExpected input form: (Real). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion == 2 && !(sp_params isa Real)
		@error "wrong params for spatial cohesion 2.\nExpected input form: (Real). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion == 3 && !(sp_params isa Vector && length.(sp_params) == [2,1,1,4])
		@error "wrong params for spatial cohesion 3.\nExpected input form: (1x2 Vector, Real, Real, 2x2 Matrix). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion == 4 && !(sp_params isa Vector && length.(sp_params) == [2,1,1,4])
		@error "wrong params for spatial cohesion 4.\nExpected input form: (1x2 Vector, Real, Real, 2x2 Matrix). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion == 5 && !(sp_params isa Real)
		@error "wrong params for spatial cohesion 5.\nExpected input form: (Real). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion == 6 && !(sp_params isa Real)
		@error "wrong params for spatial cohesion 6.\nExpected input form: (Real). Received: $(typeof(sp_params))."
		return
	end

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
	alpha_iter = ones(T_star)*starting_alpha
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
	lg_weights = zeros(n)

	############# testing on random values #############
	Si_iter = rand(collect(1:5),n,T_star)
	gamma_iter = rand((0,1),n,T_star)
	
	
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

					# get also the reduced spatial info if sPPM model
					if sPPM
						sp1_red = sp1[indexes]
						sp2_red = sp2[indexes]
					end

					# compute n_red's and nclus_red's and relabel
					n_red = length(Si_red) # "n" relative to here, i.e. the sub-partition size
					n_red1 = length(Si_red1)
					relabel!(Si_red,n_red)
					relabel!(Si_red1,n_red1)
					nclus_red = isempty(Si_red) ? 0 : maximum(Si_red) # number of clusters
					nclus_red1 = maximum(Si_red1)

					# save the label of the current working-on unit j
					j_label = Si_red1[end]

					# compute also nh_red's
					nh_red = zeros(Int64,n); nh_red1 = zeros(Int64,n)
					for jj in 1:n_red
						nh_red[Si_red[jj]] += 1 # numerosities for each cluster label
						nh_red1[Si_red1[jj]] += 1
					end
					nh_red1[Si_red1[end]] += 1 # account for the last unit j

					# start computing weights
					lCo = 0.0; lCn = 0.0
					# unit j can enter an existing cluster...
					for k in 1:nclus_red
						if sPPM
							# filter the spatial coordinates of the units of label k
							sp_idxs = findall(jj -> Si_red[jj] == k, 1:n_red)
							s1o = sp1[sp_idxs]
							s2o = sp2[sp_idxs]
							s1n = push!(copy(sp1),sp1[j])
							s2n = push!(copy(sp2),sp2[j])
							lCo = spatial_cohesion(spatial_cohesion_idx, s1o, s2o, sp_params, lg=true, M=M_dp)
							lCn = spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
						end
						lg_weights[k] = log(nh_red[k]) + lCn - lCo
					end
					# ... or unit j can enter a singleton
					lCn = 0.0
					if sPPM
						lCn = spatial_cohesion(spatial_cohesion_idx, [sp1[j]], [sp2[j]], sp_params, lg=true, M=M_dp)
					end
					lg_weights[max(nclus_red,1)] = log(M_dp) + lCn 
					# max(nclus_red,1) to account the case of nclus_red=0 (where i.e. we dont enter the above loop)
					debug(@showd lg_weights)

					# now use the weights towards sampling the new gamma_jt
					ph_tmp = copy(lg_weights)
					sort!(ph_tmp)
					max_ph = ph_tmp[end]
					sum_ph = 0.0
					# exponentiate...
					for k in 1:nclus_red
						 # for numeric purposes we subract max_ph
						lg_weights[k] = exp(lg_weights[k] - max_ph)
						sum_ph += lg_weights[k]
					end
					# ... and normalize
					for k in 1:(nclus_red + 1)
						lg_weights[k] /= sum_ph
					end
					# compute probh
					probh = alpha_iter[t] / (alpha_iter[t] + (1 - alpha_iter[t]) * lg_weights[j_label])
					debug(@showd probh alpha_iter[t])

					# compatibility check for gamma transition
					if gamma_iter[j, t] == 0
						# we want to find ρ_(t-1)^{R_t(+j)} ...
						indexes = findall(jj -> jj != j && gamma_iter[jj, t] == 1, 1:n)
						Si_comp1 = Si_iter[indexes, t-1]
						Si_comp2 = Si_iter[indexes, t]
						push!(Si_comp1,Si_iter[j,t-1])
						push!(Si_comp2,Si_iter[j,t]) # ... and ρ_t^R_t(+j)}

						rho_comp = compatibility(Si_comp1, Si_comp2)
						if rho_comp == 0
							probh = 0
						end
					end
					# Sample the new gamma
					gt = rand(Bernoulli(probh))
					gamma_iter[j, t] = gt
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
	
close(log_file)
# global_logger(ConsoleLogger())	
end

