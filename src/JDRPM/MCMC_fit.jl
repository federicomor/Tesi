using Distributions
using Statistics
using LinearAlgebra
using Random
using Logging
using Dates
# using TimerOutputs
using ProgressMeter
using StaticArrays
using Printf
using MCMCChains
using DataFrames
using CSV

log_file = open("log.txt", "w+")
include("debug.jl")
include("utils.jl")

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
	lambda2_priors::Vector{Float64},      # Prior parameters for lambda2 ∼ invGamma(a_lambda=..., b_lambda=...)
	alpha_priors::AbstractArray{Float64}, # Prior parameters for alpha ∼ Beta(a_alpha=..., b_alpha=...)
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
	Random.seed!(round(Int64,seed))

	############# define auxiliary variables #############
	# missing is for when we dont provide the argument when fitting
	# nothing is for when we provide it but is NULL, on R
	sPPM = !ismissing(sp_coords) && !isnothing(sp_coords)
	cl_xPPM = !ismissing(Xcl_covariates) && !isnothing(Xcl_covariates)
	lk_xPPM = !ismissing(Xlk_covariates) && !isnothing(Xcl_covariates)
	n, T = size(Y)
	T_star = T+1
	p_lk = lk_xPPM ? size(Xlk_covariates)[2] : 0
	p_cl = cl_xPPM ? size(Xcl_covariates)[2] : 0
	Y_has_NA = any(ismissing.(Y))
	nout = round(Int64, (draws - burnin)/(thin))

	if sPPM
		sp1::Vector{Float64} = copy(vec(sp_coords[:,1]))
		sp2::Vector{Float64} = copy(vec(sp_coords[:,2]))

		S=@MMatrix zeros(2, 2)
		sp_params_struct = SpParams(
			alpha = spatial_cohesion_idx==1 ? sp_params[1] : 1.,
			a = spatial_cohesion_idx==2 ? sp_params[1] : 1.,
			mu_0 = (spatial_cohesion_idx==3 || spatial_cohesion_idx==4) ? SVector{2}(sp_params[1]...) : @SVector(zeros(2)),
			k0 = (spatial_cohesion_idx==3 || spatial_cohesion_idx==4) ? sp_params[2] : 1.,
			v0 = (spatial_cohesion_idx==3 || spatial_cohesion_idx==4) ? sp_params[3] : 1.,
			Psi = (spatial_cohesion_idx==3 || spatial_cohesion_idx==4) ? SMatrix{2,2}(sp_params[4]...) : @SMatrix(zeros(2,2)),
			phi = (spatial_cohesion_idx==5 || spatial_cohesion_idx==6) ? sp_params[1] : 1.
		)
	end

	if cl_xPPM
		# provide different a's and b's for each p covariate
		cv_params_sim4 = Dict()
		if covariate_similarity_idx == 4
			for i in 1:2:length(cv_params)-2
				cv_params_sim4[i ÷ 2 + 1] = [cv_params[1], cv_params[2], cv_params[i+2], cv_params[i+3]]
			end
		end

		Rs = zeros(p_cl,T)
		if covariate_similarity==2 || covariate_similarity==3
			for p in 1:p_cl
				for t in 1:T
					if isa(Xcl_covariates[1,p,t],Real)
						Rs[p,t] = maximum(Xcl_covariates[:,p,t])-minimum(Xcl_covariates[:,p,t])
					end
				end
			end
		end
	end
	
	if update_eta1 acceptance_ratio_eta1 = 0 end
	if update_phi1 acceptance_ratio_phi1 = 0 end

	if !skip_checks
	if logging
		println("Logging to file:\n", abspath("log.txt"),"\n")
		printlgln(replace(string(now()),"T" => "   "))
		debug("current seed = $seed")
	end


	############# check some stuff #############
	if sPPM
		if ismissing(sp_params)
			@error "Please provide sp_params if you want to fit with spatial information." _file=""
			return
		end
		if !(sp_params isa Vector)
			@error "The sp_params are required to be passed in julia Vector form (i.e. list on R)." _file=""
			return
		end
		if spatial_cohesion_idx == 1 && length.(sp_params) != [1] 
			@error "Wrong params for spatial cohesion 1.\nExpected input form: [Real]." _file=""
			return
		elseif spatial_cohesion_idx == 2 && length.(sp_params) != [1]
			@error "Wrong params for spatial cohesion 2.\nExpected input form: [Real]." _file=""
			return
		elseif spatial_cohesion_idx == 3 && length.(sp_params) != [2,1,1,4]
			@error "Wrong params for spatial cohesion 3.\nExpected input form: [1x2 Vector, Real, Real, 2x2 Matrix]." _file=""
			return
		elseif spatial_cohesion_idx == 4 && length.(sp_params) != [2,1,1,4]
			@error "Wrong params for spatial cohesion 4.\nExpected input form: [1x2 Vector, Real, Real, 2x2 Matrix]." _file=""
			return
		elseif spatial_cohesion_idx == 5 && length.(sp_params) != [1]
			@error "Wrong params for spatial cohesion 5.\nExpected input form: [Real]." _file=""
			return
		elseif spatial_cohesion_idx == 6 && length.(sp_params) != [1]
			@error "Wrong params for spatial cohesion 6.\nExpected input form: [Real]." _file=""
			return
		end
		if (spatial_cohesion_idx == 3 || spatial_cohesion_idx == 4) && !issymmetric(sp_params[4])
			@error "Matrix Psi of the spatial parameters must be symmetric." _file=""
			return
		end
	end

	if cl_xPPM
		if !(cv_params isa Vector)
			@error "The cv_params are required to be passed in julia Vector form (i.e. list on R)." _file=""
			return
		end
		if covariate_similarity == 1 && length.(cv_params) != [1] 
			@error "Wrong params for covariate similarity 1.\nExpected input form: [Real]." _file=""
			return
		elseif covariate_similarity == 2 && length.(cv_params) != [1]
			@error "Wrong params for covariate similarity 2.\nExpected input form: [Real]." _file=""
			return
		elseif covariate_similarity == 3 && length.(cv_params) != [1]
			@error "Wrong params for covariate similarity 3.\nExpected input form: [Real]." _file=""
			return
		elseif covariate_similarity == 4 && length(cv_params) != 2+2p_cl
			@error "Wrong params for covariate similarity 4.\nExpected input form: [Real, Real, (Real, Real) for p_cl times]." _file=""
			return
		end
	end

	if lk_xPPM
		if beta_update_threshold>=draws
			@error "Cannot use such a high value of beta beta_update_threshold." _file=""
			return
		end
		if ismissing(beta_priors)
			@error "Cannot use covariates in the likelihood if beta_priors is not defined." _file=""
			return
		end
	end

	if (time_specific_alpha==false && unit_specific_alpha==false) || (time_specific_alpha==true && unit_specific_alpha==false)
		# cases of alpha being a scalar or a vector in time
		# so we should get as priors the two params of the Beta 
		if !(typeof(alpha_priors) <: Vector) || size(alpha_priors) != (2,)
			@error "Wrong params for alpha priors. \nExpected input form: (Vector) of size 2.\nReceived: $(typeof(alpha_priors)) of size $(size(alpha_priors))." _file=""
			return
		end
	elseif (time_specific_alpha==false && unit_specific_alpha==true) || (time_specific_alpha==true && unit_specific_alpha==true)
		# cases of alpha being a vector in units or a matrix
		# so we should get as priors the two params of the Beta for each unit
		if !(typeof(alpha_priors) <: Matrix) || size(alpha_priors) != (2,n)
			@error "Wrong params for alpha priors. \nExpected input form: (Matrix) of size 2*n.\nReceived: $(typeof(alpha_priors)) of size $(size(alpha_priors))." _file=""
			return
		end
	end

	if update_eta1==true && include_eta1==false 
		@error "Be coherent! I can't have 'update eta1' true and simultaneously not including it.\nCheck what you assigned to update_eta1 and include_eta1." _file=""
		return
	end
	if update_phi1==true && include_phi1==false
		@error "Be coherent! I can't have 'update phi1' true and simultaneously not including it.\nCheck what you assigned to update_phi1 and include_phi1." _file=""
		return
	end


	############# some more checks #############
	if !ismissing(initial_partition)
		if length(initial_partition) != n
			@error "The initial parition provided does not cover all the $n units of data Y." _file=""
			return
		elseif any(initial_partition .<= 0)
			@error "Labels for the initial partition provided should start from 1." _file=""
			return
		end
	end

	if (draws-burnin)%thin != 0
		@error "Please define draws, thin and burnin in an integer division-friendly way;\ni.e. such that (draws-burnin) % thin = 0." _file=""
		return
	end

	if nout <= 0 
		@error "Wrong iterations parameters. Ensure that the number of draws is high enough." _file=""
		return
	end


	############# send feedback #############
	if verbose
		println("Parameters:")
		println("sig2h ∼ InverseGamma($(sig2h_priors[1]), $(sig2h_priors[2]))")
		println("Logit(1/2(eta1+1)) ∼ Laplace(0, $(eta1_priors[1]))")
		if lk_xPPM
			println("beta ∼ MvNormal(μ=$(beta_priors[1:end-1]), Σ=$(beta_priors[end])*I)")
			println("updating beta after iteration $beta_update_threshold")
		end
		println("tau2 ∼ InverseGamma($(tau2_priors[1]), $(tau2_priors[2]))")
		println("phi0 ∼ Normal(μ=$(phi0_priors[1]), σ²=$(phi0_priors[2]))")
		println("lambda2 ∼ InverseGamma($(lambda2_priors[1]), $(lambda2_priors[2]))")
		println("alpha ∼ Beta($(alpha_priors[1]), $(alpha_priors[2]))")
		println()
	end


	println("- using seed $seed -")
	println("fitting $(Int(draws)) total iterates (with burnin=$(Int(burnin)), thinning=$(Int(thin)))")
	println("thus producing $nout valid iterates in the end\n")
	println("on n=$n subjects\nfor T=$T time instants\n")
	println(sPPM ? "[✓]" : "[✗]", " with space? $sPPM", sPPM ? " (cohesion $spatial_cohesion_idx)" : "")
	# if sPPM println("- params: ", sp_params_real) end
	println(lk_xPPM ? "[✓]" : "[✗]", " with covariates in the likelihood? $lk_xPPM", lk_xPPM ? " (p=$p_lk)" : "")
	println(cl_xPPM ? "[✓]" : "[✗]", " with covariates in the clustering process? $cl_xPPM", cl_xPPM ? " (p=$p_cl, similarity $covariate_similarity_idx)" : "")
	# if cl_xPPM println("- params: ", cv_params) end
	println(Y_has_NA ? "[✓]" : "[✗]", " are there missing data in Y? $Y_has_NA")
	println()
	end # of the skip_checks if



	############# update to handle missing data #############
	# to remember which units and at which times had a missing value
	missing_map = ismissing.(Y)
	# eg will be something like this
	# julia> Y
	# 6×4 Matrix{Union{Missing, Float64}}:
	#  1.2        missing  5.0  6.0
	#  3.4       5.5       1.0  1.0
	#   missing  1.0       1.0  3.0
	#   missing  5.0       5.0  7.0
	#  4.0       4.0       4.0  1.0
	#  1.0       1.0       1.0  0.0
	# 
	# julia> ismissing.(Y) # the missing_map
	# 6×4 BitMatrix:
	#  0  1  0  0
	#  0  0  0  0
	#  1  0  0  0
	#  1  0  0  0
	#  0  0  0  0
	#  0  0  0  0
	# so with the map we note that if we encounter a 1 when working on [j,t] we have to simulate the Y value
	# The map is needed since when we simulate it it will no longer be missing, and therefore we could "lose track" of him
	missing_idxs = Tuple.(findall(missing_map))

	############# allocate output variables #############
	i_out = 1
	Si_out = zeros(Int64,n,T,nout)
	gamma_out = zeros(Bool,n,T,nout)
	if time_specific_alpha==false && unit_specific_alpha==false
		# for each iterate, a scalar
		alpha_out = zeros(nout)
	elseif time_specific_alpha==true && unit_specific_alpha==false
		# for each iterate, a vector in time
		alpha_out = zeros(T,nout)
	elseif time_specific_alpha==false && unit_specific_alpha==true
		# for each iterate, a vector in units
		alpha_out = zeros(n,nout)
	elseif time_specific_alpha==true && unit_specific_alpha==true
		# for each iterate, a matrix
		alpha_out = zeros(n,T,nout)
	end
	sigma2h_out = zeros(n,T,nout)
	muh_out = zeros(n,T,nout)
	eta1_out = zeros(n,nout)
	if lk_xPPM
		beta_out = zeros(T,p_lk,nout)
		# so that eg we can write beta_out[t,:,i] = beta_iter[t]
	end
	theta_out = zeros(T,nout)
	tau2_out = zeros(T,nout)
	phi0_out = zeros(nout)
	phi1_out = zeros(nout)
	lambda2_out = zeros(nout)
	fitted = zeros(n,T,nout)
	llike = zeros(n,T,nout)
	LPML = 0.0
	WAIC = 0.0

	CPO = zeros(n,T)
	mean_likelhd = zeros(n,T)
	mean_loglikelhd = zeros(n,T)


	############# allocate auxiliary working variables #############
	aux1_relab = zeros(Int64,n)
	aux2_relab = zeros(Int64,n)
	Si_relab = zeros(Int64,n)
	nh_reorder = zeros(Int64,n)
	old_lab = zeros(Int64,n)
	nh_red = zeros(Int64,n)
	nh_red1 = zeros(Int64,n)
	lg_weights = zeros(n)
	ph = zeros(n)
	muh_iter_copy = zeros(n,T_star)
	sig2h_iter_copy = ones(n,T_star)
	log_Mdp = log(M_dp)
	lC = @MVector zeros(2)
	lS = @MVector zeros(2)
	lPP = @MVector zeros(1)
	Si_red1 = zeros(Int64,n)
	s1n = zeros(n)
	s1o = zeros(n)
	s2n = zeros(n)
	s2o = zeros(n)
	nh_tmp = zeros(Int64,n)
	rho_tmp = zeros(Int64,n)
	Xo = zeros(n)
	Xn = zeros(n)
	Xo_cat = Vector{String}(undef, n)
	Xn_cat = Vector{String}(undef, n)

	a_star_tau = 0
	b_star_tau = 0
	a_star_lambda2 = 0
	b_star_lambda2 = 0


	############# allocate and initialize working variables #############
	Si_iter = ones(Int64,n,T_star) # label assignements for units j at time t
	Si_iter[:,end] .= 0
	gamma_iter = zeros(Bool,n,T_star)
	if !ismissing(initial_partition)
		relabel!(initial_partition,n)
		Si_iter[:,1] = initial_partition
		gamma_iter[:,1] .= 1
	end
	if time_specific_alpha==false && unit_specific_alpha==false
		# a scalar
		alpha_iter = starting_alpha
	elseif time_specific_alpha==true && unit_specific_alpha==false
		# a vector in time
		alpha_iter = ones(T_star)*starting_alpha
	elseif time_specific_alpha==false && unit_specific_alpha==true
		# a vector in units
		alpha_iter = ones(n)*starting_alpha
	elseif time_specific_alpha==true && unit_specific_alpha==true
		# a matrix
		alpha_iter = ones(n,T_star)*starting_alpha
	end

	# hierarchy level 3
	phi0_iter = rand(Normal(phi0_priors[1], sqrt(phi0_priors[2])))
	phi1_iter = 0.0 
	if include_phi1
		phi1_iter = rand(Uniform(-1,1))
	end
	# lambda2_iter = rand(InverseGamma(lambda2_priors...))
	lambda2_iter = 1

	# hierarchy level 2
	theta_iter = zeros(T_star)
	theta_iter[1] = rand(Normal(phi0_iter,sqrt(lambda2_iter)))
	@inbounds for t in 2:T_star
		theta_iter[t] = rand(Normal((1-phi1_iter)*phi0_iter + phi1_iter*theta_iter[t-1],sqrt(lambda2_iter*(1-phi1_iter^2))))
	end
	tau2_iter = zeros(T_star)
	@inbounds for t in 1:T_star
		# tau2_iter[t] = rand(InverseGamma(tau2_priors...))
		tau2_iter[t] = 1
	end

	# hierarchy level 1
	sig2h_iter = ones(n,T_star)
	muh_iter = zeros(n,T_star)
	if lk_xPPM
		beta_iter = Vector{Vector{Float64}}(undef,T)
		beta0 = beta_priors[1:end-1]
		s2_beta = beta_priors[end]
		@inbounds for t in 1:T
			# beta_iter[t] = rand(MvNormal(beta0, s2_beta*I(p_lk)))
			beta_iter[t] = beta0 # initialize with the mean? maybe it's better
		end
	end

	eta1_iter = zeros(n)
	if include_eta1
		@inbounds for j in 1:n
			eta1_iter[j] = rand(Uniform(-1,1))
			# eta1_iter[j] = 0.
		end
	end

	nh = zeros(Int,n,T_star) # numerosity of each cluster k at time t
	# dimension n since at most there are n clusters (all singletons)
	# nh[1,:] .= n
	for j in 1:n
		for t in 1:T
			nh[Si_iter[j,t],t] += 1
		end
	end
	nclus_iter = zeros(Int,T_star) # how many clusters there are at time t
	for t in 1:T
		nclus_iter[t] = maximum(@view Si_iter[:,t])
	end



	############# start MCMC algorithm #############
	println(replace(string(now()),"T" => " ")[1:end-4])
	println("Starting MCMC algorithm")

	t_start = now()
	progresso = Progress(round(Int64(draws)),
			showspeed=true,
			output=stdout, # default is stderr, which turns out in orange color on R
			dt=1, # every how many seconds update the feedback
			barlen=0 # no progress bar
			)

	@inbounds for i in 1:draws
	# for i in 1:draws
		# print("iteration $i of $draws\r") # use only this when all finished, for a shorter feedback
		# there is also the ProgressMeter option, maybe is cooler

		############# sample the missing values #############
		# from the "update rho" section onwards also the Y[j,t] will be needed (to compute weights, update laws, etc)
		# so we need now to simulate the values for the data which are missing (from their full conditional)
		if Y_has_NA
			# we have to use the missing_idxs to remember which units and at which times had a missing value,
			# in order to simulate just them and instead use the given value for the other units and times
			for (j,t) in missing_idxs
				# Problem: if when filling a NA we occur in another NA value? eg when we also need Y[j,t±1]
				# I decided here to set that value to 0, if occurs, since anyway target should be centered
				# so it seems a reasonable patch
				# We could have decided to ignore this computation and just use the likelihood as proposal
				# filling distribution, but this would have just worked in the Y[j,t+1] case so the general
				# problem would have still been there

				# if i==1 println("(j=$j,t=$t)") end
				c_it = Si_iter[j,t]
				Xlk_term_t = (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0)
				aux1 = eta1_iter[j]^2


				if t==1	
					c_itp1 = Si_iter[j,t+1]
					Xlk_term_tp1 = (lk_xPPM ? dot(view(Xlk_covariates,j,:,t+1), beta_iter[t+1]) : 0)

					sig2_post = 1 / (1/sig2h_iter[c_it,t] + aux1/(sig2h_iter[c_itp1,t+1]*(1-aux1)))
					mu_post = sig2_post * ( 
						(1/sig2h_iter[c_it,t])*(muh_iter[c_it,t] + Xlk_term_t) +
						(eta1_iter[j]/(sig2h_iter[c_itp1,t+1]*(1-aux1)))*((ismissing(Y[j,t+1]) ? 0 : Y[j,t+1]) - muh_iter[c_itp1,t+1] - Xlk_term_tp1)
						)

					Y[j,t] = rand(Normal(mu_post,sqrt(sig2_post)))

				elseif 1<t<T
					c_itp1 = Si_iter[j,t+1]
					Xlk_term_tp1 = (lk_xPPM ? dot(view(Xlk_covariates,j,:,t+1), beta_iter[t+1]) : 0)

					sig2_post = (1-aux1) / (1/sig2h_iter[c_it,t] + aux1/sig2h_iter[c_itp1,t+1])
					mu_post = sig2_post * ( 
						(1/(sig2h_iter[c_it,t]*(1-aux1)))*(muh_iter[c_it,t] + eta1_iter[j]*(ismissing(Y[j,t-1]) ? 0 : Y[j,t-1]) + Xlk_term_t) +
						(eta1_iter[j]/(sig2h_iter[c_itp1,t+1]*(1-aux1)))*((ismissing(Y[j,t+1]) ? 0 : Y[j,t+1]) - muh_iter[c_itp1,t+1] - Xlk_term_tp1)
						)

					Y[j,t] = rand(Normal(mu_post,sqrt(sig2_post)))

				else # t==T
					Y[j,t] = rand(Normal(
						muh_iter[c_it,t] + eta1_iter[j]*(ismissing(Y[j,t-1]) ? 0 : Y[j,t-1]) + Xlk_term_t,
						sqrt(sig2h_iter[c_it,t]*(1-aux1))
						))
				end
			end
		end


		for t in 1:T
			############# update gamma #############
			for j in 1:n
				if t==1 
					gamma_iter[j,t] = 0
					# at the first time units get reallocated
				else
					# we want to find ρ_t^{R_t(-j)} ...
					indexes = findall_faster(jj -> jj != j && gamma_iter[jj, t] == 1, 1:n)
					Si_red = Si_iter[indexes, t]
					copy!(Si_red1, Si_red)
					push!(Si_red1, Si_iter[j,t]) # ... and ρ_t^R_t(+j)}

					# get also the reduced spatial info if sPPM model
					if sPPM
						sp1_red = @view sp1[indexes]
						sp2_red = @view sp2[indexes]
					end
					# and the reduced covariates info if cl_xPPM model
					if cl_xPPM
						Xcl_covariates_red = @view Xcl_covariates[indexes,:,t]
					end

					# compute n_red's and nclus_red's and relabel
					n_red = length(Si_red) # = "n" relative to here, i.e. the sub-partition size
					n_red1 = length(Si_red1)
					# @assert n_red1 == n_red+1
					relabel!(Si_red,n_red)
					relabel!(Si_red1,n_red1)
					nclus_red = isempty(Si_red) ? 0 : maximum(Si_red) # = number of clusters
					nclus_red1 = maximum(Si_red1)

					# save the label of the current working-on unit j
					j_label = Si_red1[end]

					# compute also nh_red's
					nh_red .= 0
					nh_red1 .= 0
					for jj in 1:n_red
						nh_red[Si_red[jj]] += 1 # = numerosities for each cluster label
						nh_red1[Si_red1[jj]] += 1
					end
					nh_red1[Si_red1[end]] += 1 # account for the last added unit j, not included in the above loop

					# start computing weights
					lg_weights .= 0

					# unit j can enter an existing cluster...
					for k in 1:nclus_red

						# filter the covariates of the units of label k
						aux_idxs = findall(Si_red .== k) # fast
						lC .= 0.
						if sPPM
							# filter the spatial coordinates of the units of label k
							copy!(s1o, sp1_red[aux_idxs])
							copy!(s2o, sp2_red[aux_idxs])
							copy!(s1n,s1o); push!(s1n, sp1[j])
							copy!(s2n,s2o); push!(s2n, sp2[j])
							spatial_cohesion!(spatial_cohesion_idx, s1o, s2o, sp_params_struct, true, M_dp, S,1,false,lC)
							spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_struct, true, M_dp, S,2,false,lC)
						end

						lS .= 0.
						if cl_xPPM
						# Xcl_covariates is a n*p*T matrix
							for p in 1:p_cl
								if isa(first(Xcl_covariates[j,p,t]),Real) # numerical covariate
									copy!(Xo, @view Xcl_covariates_red[aux_idxs,p])
									copy!(Xn, Xo); push!(Xn,Xcl_covariates[j,p,t])
									
									if covariate_similarity_idx == 4
										covariate_similarity!(covariate_similarity_idx, Xo, cv_params_sim4[p], Rs[p,t], true,1,true,lS,cv_weight)
										covariate_similarity!(covariate_similarity_idx, Xn, cv_params_sim4[p], Rs[p,t], true,2,true,lS,cv_weight)
									else 
										covariate_similarity!(covariate_similarity_idx, Xo, cv_params, Rs[p,t], true,1,true,lS,cv_weight)
										covariate_similarity!(covariate_similarity_idx, Xn, cv_params, Rs[p,t], true,2,true,lS,cv_weight)
									end
								else # categorical covariate
									copy!(Xo_cat, @view Xcl_covariates_red[aux_idxs,p])
									copy!(Xn_cat, Xo_cat); push!(Xn_cat,Xcl_covariates[j,p,t])
									covariate_similarity!(covariate_similarity_idx, Xo_cat, cv_params, Rs[p,t], true,1,true,lS,cv_weight)
									covariate_similarity!(covariate_similarity_idx, Xn_cat, cv_params, Rs[p,t], true,2,true,lS,cv_weight)
								end
							end
						end

						lg_weights[k] = log(nh_red[k]) + lC[2] - lC[1] + lS[2] - lS[1]
					end
					
					# ... or unit j can create a singleton
					lC .= 0.
					if sPPM
						spatial_cohesion!(spatial_cohesion_idx, SVector(sp1[j]), SVector(sp2[j]), sp_params_struct, true, M_dp, S,2,false,lC)
					end
					lS .= 0.
					if cl_xPPM
						for p in 1:p_cl
							if covariate_similarity_idx == 4
								covariate_similarity!(covariate_similarity_idx, SVector(Xcl_covariates[j,p,t]), cv_params_sim4[p], Rs[p,t], true, 2,true,lS,cv_weight)
							else
								covariate_similarity!(covariate_similarity_idx, SVector(Xcl_covariates[j,p,t]), cv_params, Rs[p,t], true, 2,true,lS,cv_weight)
							end
						end
					end
					lg_weights[nclus_red+1] = log_Mdp + lC[2] + lS[2]


					# now use the weights towards sampling the new gamma_jt
					max_ph = maximum(@view lg_weights[1:(nclus_red+1)])
					sum_ph = 0.0

					# exponentiate...
					for k in 1:(nclus_red+1)
						 # for numerical purposes we subract max_ph
						lg_weights[k] = exp(lg_weights[k] - max_ph)
						sum_ph += lg_weights[k]
					end
					# ... and normalize
					lg_weights ./= sum_ph

					# compute probh
					probh::Float64 = 0.0
					if time_specific_alpha==false && unit_specific_alpha==false
						probh = alpha_iter / (alpha_iter + (1 - alpha_iter) * lg_weights[j_label])
					elseif time_specific_alpha==true && unit_specific_alpha==false
						probh = alpha_iter[t] / (alpha_iter[t] + (1 - alpha_iter[t]) * lg_weights[j_label])
					elseif time_specific_alpha==false && unit_specific_alpha==true
						probh = alpha_iter[j] / (alpha_iter[j] + (1 - alpha_iter[j]) * lg_weights[j_label])
					elseif time_specific_alpha==true && unit_specific_alpha==true
						probh = alpha_iter[j,t] / (alpha_iter[j,t] + (1 - alpha_iter[j,t]) * lg_weights[j_label])
					end


					# compatibility check for gamma transition
					if gamma_iter[j, t] == 0
						# we want to find ρ_(t-1)^{R_t(+j)} ...
						indexes = findall_faster(jj -> jj==j || gamma_iter[jj, t]==1, 1:n)
						Si_comp1 = @view Si_iter[indexes, t-1]
						Si_comp2 = @view Si_iter[indexes, t] # ... and ρ_t^R_t(+j)}

						rho_comp = compatibility(Si_comp1, Si_comp2)
						if rho_comp == 0
							probh = 0.0
						end
					end
					# sample the new gamma
					gt = rand(Bernoulli(probh))
					gamma_iter[j, t] = gt
				end
			end # for j in 1:n


			############# update rho #############
			# we only update the partition for the units which can move (i.e. with gamma_jt=0)
			movable_units = findall(gamma_iter[:,t] .== 0) # fast
		
			for j in movable_units
				# remove unit j from the cluster she is currently in

				if nh[Si_iter[j,t],t] > 1 # unit j does not belong to a singleton cluster
					nh[Si_iter[j,t],t] -= 1
					# no nclus_iter[t] change since j's cluster is still alive
				else # unit j belongs to a singleton cluster
					j_label = Si_iter[j,t]
					last_label = nclus_iter[t]

					if j_label < last_label
						# here we enter if j_label is not the last label, so we need to
						# relabel clusters in order to then remove j's cluster
						# eg: units 1 2 3 4 5 j 7 -> units 1 2 3 4 5 j 7
						#     label 1 1 2 2 2 3 4    label 1 1 2 2 2 4 3

						# swap cluster labels...
						for jj in 1:n
							if Si_iter[jj, t] == last_label
								Si_iter[jj, t] = j_label
							end
						end
						Si_iter[j, t] = last_label
						# ... and cluster-specific parameters
						sig2h_iter[j_label, t], sig2h_iter[last_label, t] = sig2h_iter[last_label, t], sig2h_iter[j_label, t]
						muh_iter[j_label, t], muh_iter[last_label, t] = muh_iter[last_label, t], muh_iter[j_label, t]
						nh[j_label, t] = nh[last_label, t]
						nh[last_label, t] = 1

					end
					# remove the j-th observation and the last cluster (being j in a singleton)
					nh[last_label, t] -= 1
					nclus_iter[t] -= 1
				end

				# setup probability weights towards the sampling of rho_jt
				ph .= 0.0
				resize!(ph,nclus_iter[t]+1)
				copy!(rho_tmp, @view Si_iter[:,t])

				# compute nh_tmp (numerosities for each cluster label)
				copy!(nh_tmp, @view nh[:,t])
				# unit j contribute is already absent from the change we did above
				nclus_temp = 0

					
				# we now simulate the unit j to be assigned to one of the existing clusters...
				for k in 1:nclus_iter[t]
					rho_tmp[j] = k
					indexes = findall(gamma_iter[:,t+1] .== 1) # fast
					# we check the compatibility between ρ_t^{h=k,R_(t+1)} ...
					Si_comp1 = @view rho_tmp[indexes]
					Si_comp2 = @view Si_iter[indexes,t+1] # and ρ_(t+1)^{R_(t+1)}
					rho_comp = compatibility(Si_comp1, Si_comp2)
					
					if rho_comp != 1
						ph[k] = log(0) # assignment to cluster k is not compatible
					else
						# update params for "rho_jt = k" simulation
						nh_tmp[k] += 1
						# nclus_temp = sum(nh_tmp .> 0) # = number of clusters
						nclus_temp = count(a->(a>0), nh_tmp) # similarly fast, not clear
						# count is a bit slower but does not allocate

						lPP .= 0.
						for kk in 1:nclus_temp
							aux_idxs = findall(rho_tmp .== kk) # fast
							if sPPM
								copy!(s1n, @view sp1[aux_idxs])
								copy!(s2n, @view sp2[aux_idxs])

								spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_struct, true, M_dp, S,1,true,lPP)
								# LPP += spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_struct, true, M_dp, S,1,true)
							end
							if cl_xPPM
								for p in 1:p_cl
									Xn_view = @view Xcl_covariates[aux_idxs,p,t]
									if covariate_similarity_idx == 4
										covariate_similarity!(covariate_similarity_idx, Xn_view, cv_params_sim4[p], Rs[p,t], true,1,true,lPP,cv_weight)
									else
										covariate_similarity!(covariate_similarity_idx, Xn_view, cv_params, Rs[p,t], true,1,true,lPP,cv_weight)
									end
								end
							end
							lPP[1] += log_Mdp + lgamma(nh_tmp[kk])
						end
						
						# lpp += log(M_dp) + lgamma(length(indexes)) # same
						### debug case
						# ph[k] = 0.5
						### real case
						if t==1
							ph[k] = loglikelihood(Normal(
								muh_iter[k,t] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
								sqrt(sig2h_iter[k,t])),
								Y[j,t]) + lPP[1]
						else
							ph[k] = loglikelihood(Normal(
								muh_iter[k,t] + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
								sqrt(sig2h_iter[k,t]*(1-eta1_iter[j]^2))),
								Y[j,t]) + lPP[1]
						end

						# restore params after "rho_jt = k" simulation
						nh_tmp[k] -= 1
					end
				end
				# ... plus the case of being assigned to a new (singleton for now) cluster 
				k = nclus_iter[t]+1
				rho_tmp[j] = k
				# declare (for later scope accessibility) the new params here
				muh_draw = 0.0; sig2h_draw = 0.0

				indexes = findall(gamma_iter[:,t+1] .== 1)
				Si_comp1 = @view rho_tmp[indexes]
				Si_comp2 = @view Si_iter[indexes,t+1]
				rho_comp = compatibility(Si_comp1, Si_comp2)

				if rho_comp != 1
					ph[k] = log(0) # assignment to a new cluster is not compatible
				else
					# sample new params for this new cluster
					muh_draw = rand(Normal(theta_iter[t], sqrt(tau2_iter[t])))
					sig2h_draw = rand(InverseGamma(sig2h_priors[1],sig2h_priors[2]))

					# update params for "rho_jt = k" simulation
					nh_tmp[k] += 1
					# nclus_temp = sum(nh_tmp .> 0)
					nclus_temp = count(a->(a>0), nh_tmp) # similarly fast, not clear
					# count is a bit slower but does not allocate

					lPP .= 0.
					for kk in 1:nclus_temp
						aux_idxs = findall(rho_tmp .== kk) # fast
						if sPPM
							copy!(s1n, @view sp1[aux_idxs])
							copy!(s2n, @view sp2[aux_idxs])
							spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_struct, true, M_dp, S,1,true,lPP)
						end
						if cl_xPPM
							for p in 1:p_cl
								Xn_view = @view Xcl_covariates[aux_idxs,p,t]
								if covariate_similarity_idx == 4
									covariate_similarity!(covariate_similarity_idx, Xn_view, cv_params_sim4[p], Rs[p,t], true,1,true,lPP,cv_weight)
								else
									covariate_similarity!(covariate_similarity_idx, Xn_view, cv_params, Rs[p,t], true,1,true,lPP,cv_weight)
								end
							end
						end
						lPP[1] += log_Mdp + lgamma(nh_tmp[kk])
						# lpp += log(M_dp) + lgamma(length(indexes)) # same
					end
					### debug case
					# ph[k] = 0.5
					### real case
					if t==1
				ph[k] = loglikelihood(Normal(
					muh_draw + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
					sqrt(sig2h_draw)),
					Y[j,t]) + lPP[1]
					else
						ph[k] = loglikelihood(Normal(
							muh_draw + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_draw*(1-eta1_iter[j]^2))),
							Y[j,t]) + lPP[1]
					end

					# restore params after "rho_jt = k" simulation
					nh_tmp[k] -= 1
				end

				# now exponentiate the weights...
				max_ph = maximum(ph)
				sum_ph = 0.0
				for k in eachindex(ph)
					# for numerical purposes we subract max_ph
					ph[k] = exp(ph[k] - max_ph)
					sum_ph += ph[k]
				end
				# ... and normalize them
				ph ./= sum_ph
				

				# now sample the new label Si_iter[j,t]
				u = rand(Uniform(0,1))
				cph = cumsum(ph)
				cph[end] = 1 # fix numerical problems of having sums like 0.999999etc
				new_label = 0
				for k in eachindex(ph)
					if u <= cph[k]
						new_label = k
						break
					end
				end
				
				if new_label <= nclus_iter[t]
					# we enter an existing cluster
					Si_iter[j, t] = new_label
					nh[new_label, t] += 1
					# rho_tmp[j] = new_label # useless
				else
					# we create a new singleton cluster
					nclus_iter[t] += 1
					cl_new = nclus_iter[t]
					Si_iter[j, t] = cl_new
					nh[cl_new, t] = 1
					# rho_tmp[j] = cl_new # useless
					muh_iter[cl_new, t] = muh_draw
					sig2h_iter[cl_new, t] = sig2h_draw
				end

				# now we need to relabel after the possible mess created by the sampling
				# eg: (before sampling)  (after sampling)
				#     units j 2 3 4 5 ->  units j 2 3 4 5
				#    labels 1 1 1 2 2    labels 3 1 1 2 2
				# the after case has to be relabelled
				Si_tmp = @view Si_iter[:,t]

				# Si_relab, nh_reorder, old_lab = relabel_full(Si_tmp,n)				
				relabel_full!(Si_tmp,n,Si_relab, nh_reorder, old_lab)				
				# - Si_relab gives the relabelled partition
				# - nh_reorder gives the numerosities of the relabelled partition, ie "nh_reorder[k] = #(units of new cluster k)"
				# - old_lab tells "the index in position i (which before was cluster i) is now called cluster old_lab[i]"
				# eg:             Original labels (Si): 4 2 1 1 1 3 1 4 5 
				#          Relabeled groups (Si_relab): 1 2 3 3 3 4 3 1 5
				# Reordered cluster sizes (nh_reorder): 2 1 4 1 1 0 0 0 0
				# 	              Old labels (old_lab): 4 2 1 3 5 0 0 0 0 

				# now fix everything (morally permute params)
				Si_iter[:,t] = Si_relab
				copy!(muh_iter_copy, muh_iter) # copy!(dst,src)
				copy!(sig2h_iter_copy, sig2h_iter) # copy!(dst,src)
				len = findlast(x -> x != 0, nh_reorder)
				for k in 1:nclus_iter[t]
					muh_iter[k,t] = muh_iter_copy[old_lab[k],t]
					sig2h_iter[k,t] = sig2h_iter_copy[old_lab[k],t]
					nh[k,t] = nh_reorder[k]
				end

			end # for j in movable_units


			############# update muh #############
			if t==1
				for k in 1:nclus_iter[t]
					sum_Y = 0.0
					for j in 1:n
						if Si_iter[j,t]==k 
							sum_Y += Y[j,t] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0.0)
						end
					end
					sig2_star = 1 / (1/tau2_iter[t] + nh[k,t]/sig2h_iter[k,t])
					mu_star = sig2_star * (theta_iter[t]/tau2_iter[t] + sum_Y/sig2h_iter[k,t])

					muh_iter[k,t] = rand(Normal(mu_star,sqrt(sig2_star)))
				end

			else # t>1
				for k in 1:nclus_iter[t]
					sum_Y = 0.0
					sum_e2 = 0.0
					for j in 1:n
						if Si_iter[j,t]==k 
							aux1 = 1 / (1-eta1_iter[j]^2)
							sum_e2 += aux1
							sum_Y += (Y[j,t] - eta1_iter[j]*Y[j,t-1] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0.0)) * aux1
						end
					end
					sig2_star = 1 / (1/tau2_iter[t] + sum_e2/sig2h_iter[k,t]) 
					mu_star = sig2_star * (theta_iter[t]/tau2_iter[t] + sum_Y/sig2h_iter[k,t])

					muh_iter[k,t] = rand(Normal(mu_star,sqrt(sig2_star)))
				end
			end
			
			############# update sigma2h #############
			if t==1
				for k in 1:nclus_iter[t]
					a_star = sig2h_priors[1] + nh[k,t]/2 # should be the same
					sum_Y = 0.0
					S_kt = findall(Si_iter[:,t] .== k) # fast
					for j in S_kt
						sum_Y += (Y[j,t] - muh_iter[k,t] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0.0))^2
					end

					b_star = sig2h_priors[2] + sum_Y/2
					sig2h_iter[k,t] = rand(InverseGamma(a_star, b_star))
				end

			else # t>1
				for k in 1:nclus_iter[t]
					a_star = sig2h_priors[1] + nh[k,t]/2
					sum_Y = 0.0
					S_kt = findall(Si_iter[:,t] .== k) # fast
					for j in S_kt
						sum_Y += (Y[j,t] - muh_iter[k,t] - eta1_iter[j]*Y[j,t-1] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0.0))^2
					end

					b_star = sig2h_priors[2] + sum_Y/2
					sig2h_iter[k,t] = rand(InverseGamma(a_star, b_star))
				end
			end
			############# update beta #############
			if lk_xPPM && i>=beta_update_threshold
				if t==1
					sum_Y = zeros(p_lk)
					A_star = I(p_lk)/s2_beta
					for j in 1:n
						X_jt = @view Xlk_covariates[j,:,t]
						sum_Y += (Y[j,t] - muh_iter[Si_iter[j,t],t]) * X_jt / sig2h_iter[Si_iter[j,t],t]
						A_star += (X_jt * X_jt') / sig2h_iter[Si_iter[j,t],t]
					end
					b_star = beta0/s2_beta + sum_Y
					A_star = Symmetric(A_star)
					# Symmetric is needed for numerical problems
					# https://discourse.julialang.org/t/isposdef-and-eigvals-do-not-agree/118191/2
					# but A_star is indeed symm and pos def (by construction) so there is no problem
					beta_iter[t] = rand(MvNormalCanon(b_star, A_star)) # quicker and more accurate method
					# Am1_star = inv(A_star) # old method with the MvNormal and the inversion required
					# beta_iter[t] = rand(MvNormal(inv(Symmetric(A_star))*b_star, inv(Symmetric(A_star))))
				else
					sum_Y = zeros(p_lk)
					A_star = I(p_lk)/s2_beta
					for j in 1:n
						X_jt = @view Xlk_covariates[j,:,t]
						sum_Y += (Y[j,t] - muh_iter[Si_iter[j,t],t] - eta1_iter[j]*Y[j,t-1]) * X_jt / sig2h_iter[Si_iter[j,t],t]
						A_star += (X_jt * X_jt') / sig2h_iter[Si_iter[j,t],t]
					end
					b_star = beta0/s2_beta + sum_Y
					A_star = Symmetric(A_star)
					beta_iter[t] = rand(MvNormalCanon(b_star, A_star)) # quicker and more accurate method
					# Am1_star = inv(A_star) # old method with the MvNormal and the inversion required
					# beta_iter[t] = rand(MvNormal(inv(Symmetric(A_star))*b_star, inv(Symmetric(A_star))))
				end
			end
						
			############# update theta #############
			aux1::Float64 = 1 / (lambda2_iter*(1-phi1_iter^2))
			kt = nclus_iter[t]
			sum_mu=0.0
			for k in 1:kt
				sum_mu += muh_iter[k,t]
			end

			if t==1
				sig2_post = 1 / (1/lambda2_iter + phi1_iter^2*aux1 + kt/tau2_iter[t])
				mu_post = sig2_post * (phi0_iter/lambda2_iter + sum_mu/tau2_iter[t] + (phi1_iter*(theta_iter[t+1] - (1-phi1_iter)*phi0_iter))*aux1)

				theta_iter[t] = rand(Normal(mu_post, sqrt(sig2_post)))

			elseif t==T
				sig2_post = 1 / (aux1 + kt/tau2_iter[t])
				mu_post = sig2_post * (sum_mu/tau2_iter[t] + ((1- phi1_iter)*phi0_iter + phi1_iter*theta_iter[t-1])*aux1)

				theta_iter[t] = rand(Normal(mu_post, sqrt(sig2_post)))

			else # 1<t<T
				sig2_post = 1 / ((1+phi1_iter^2)*aux1 + kt/tau2_iter[t])
				mu_post = sig2_post * (sum_mu/tau2_iter[t] + (phi1_iter*(theta_iter[t-1]+theta_iter[t+1]) + phi0_iter*(1-phi1_iter)^2)*aux1)

				theta_iter[t] = rand(Normal(mu_post, sqrt(sig2_post)))
			end 
			
			############# update tau2 #############
			kt = nclus_iter[t]
			aux1 = 0.0
			for k in 1:kt
				aux1 += (muh_iter[k,t] - theta_iter[t])^2 
			end
			a_star_tau = tau2_priors[1] + kt/2
			b_star_tau = tau2_priors[2] + aux1/2
			tau2_iter[t] = rand(InverseGamma(a_star_tau, b_star_tau))
			# tau2_iter[t] = rand(truncated(InverseGamma(a_star_tau, b_star_tau),0,10))

		end # for t in 1:T
		
		############# update eta1 #############
		# eta1_priors[2] = sqrt(eta1_priors[2]) # from variance to std dev
		# no, the input argument is already the std dev
		if update_eta1
			for j in 1:n
				eta1_old = eta1_iter[j]
				eta1_new = rand(Normal(eta1_old,eta1_priors[2])) # proposal value

				if (-1 <= eta1_new <= 1)
					ll_old = 0.0
					ll_new = 0.0
					for t in 2:T
						# likelihood contribution
						ll_old += loglikelihood(Normal(
							muh_iter[Si_iter[j,t],t] + eta1_old*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_iter[Si_iter[j,t],t]*(1-eta1_old^2))
							), Y[j,t]) 
						ll_new += loglikelihood(Normal(
							muh_iter[Si_iter[j,t],t] + eta1_new*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_iter[Si_iter[j,t],t]*(1-eta1_new^2))
							), Y[j,t]) 
					end
					logit_old = aux_logit(eta1_old) 
					logit_new = aux_logit(eta1_new) 

					# prior contribution
					ll_old += -log(2*eta1_priors[1]) -1/eta1_priors[1]*abs(logit_old)
					ll_new += -log(2*eta1_priors[1]) -1/eta1_priors[1]*abs(logit_new)

					# ll_ratio = min(ll_new-ll_old, 0)
					ll_ratio = ll_new-ll_old
					u = rand(Uniform(0,1))
					if (ll_ratio > log(u))
						eta1_iter[j] = eta1_new # accept the candidate
						acceptance_ratio_eta1 += 1
					end
				end
			end
		end

		############# update alpha #############
		if update_alpha
			if time_specific_alpha==false && unit_specific_alpha==false
				# a scalar
				sumg = sum(@view gamma_iter[:,1:T])
				a_star = alpha_priors[1] + sumg
				b_star = alpha_priors[2] + n*T - sumg
				alpha_iter = rand(Beta(a_star, b_star))

			elseif time_specific_alpha==true && unit_specific_alpha==false
				# a vector in time
				for t in 1:T
					sumg = sum(@view gamma_iter[:,t])
					a_star = alpha_priors[1] + sumg
					b_star = alpha_priors[2] + n - sumg
					alpha_iter[t] = rand(Beta(a_star, b_star))
				end

			elseif time_specific_alpha==false && unit_specific_alpha==true
				# a vector in units
				for j in 1:n
					sumg = sum(@view gamma_iter[j,1:T])
					a_star = alpha_priors[1,j] + sumg
					b_star = alpha_priors[2,j] + T - sumg
					alpha_iter[j] = rand(Beta(a_star, b_star))
				end
			elseif time_specific_alpha==true && unit_specific_alpha==true
				# a matrix
				for j in 1:n
					for t in 1:T
						sumg = gamma_iter[j,t] # nothing to sum in this case
						a_star = alpha_priors[1,j] + sumg
						b_star = alpha_priors[2,j] + 1 - sumg
						alpha_iter[j,t] = rand(Beta(a_star, b_star))
					end
				end
			end
		end

		############# update phi0 #############
		aux1 = 1/lambda2_iter
		aux2 = 0.0
		for t in 2:T
			aux2 += theta_iter[t] - phi1_iter*theta_iter[t-1]
		end
		sig2_post = 1 / ( 1/phi0_priors[2] + aux1 * (1 + (T-1)*(1-phi1_iter)/(1+phi1_iter)) )
		mu_post = sig2_post * ( phi0_priors[1]/phi0_priors[2] + theta_iter[1]*aux1 + aux1/(1+phi1_iter)*aux2 )
		phi0_iter = rand(Normal(mu_post, sqrt(sig2_post)))
		
		############# update phi1 #############
		# phi1_priors = sqrt(phi1_priors) # from variance to std dev
		# no, the input argument is already the std dev
		if update_phi1
			phi1_old = phi1_iter
			phi1_new = rand(Normal(phi1_old, phi1_priors)) # proposal value

			if (-1 <= phi1_new <= 1)
				ll_old = 0.0; ll_new = 0.0
				for t in 2:T
					# likelihood contribution
					ll_old += loglikelihood(Normal(
						(1-phi1_old)*phi0_iter + phi1_old*theta_iter[t-1],
						sqrt(lambda2_iter*(1-phi1_old^2))
						), theta_iter[t]) 
					ll_new += loglikelihood(Normal(
						(1-phi1_new)*phi0_iter + phi1_new*theta_iter[t-1],
						sqrt(lambda2_iter*(1-phi1_new^2))
						), theta_iter[t]) 
				end

				# prior contribution
				ll_old += loglikelihood(Uniform(-1,1), phi1_old)
				ll_new += loglikelihood(Uniform(-1,1), phi1_new)

				# ll_ratio = min(ll_new-ll_old, 0)
				ll_ratio = ll_new-ll_old
				u = rand(Uniform(0,1))
				if (ll_ratio > log(u))
					phi1_iter = phi1_new # accept the candidate
					acceptance_ratio_phi1 += 1
				end
			end
		end

		############# update lambda2 #############
		aux1 = 0.0
		for t in 2:T
			aux1 += (theta_iter[t] - (1-phi1_iter)*phi0_iter - phi1_iter*theta_iter[t-1])^2
		end
		a_star_lambda2 = lambda2_priors[1] + T/2
		b_star_lambda2 = lambda2_priors[2] + ((theta_iter[1] - phi0_iter)^2 + aux1) / 2
		lambda2_iter = rand(InverseGamma(a_star_lambda2,b_star_lambda2))	
		# lambda2_iter = rand(truncated(InverseGamma(a_star_lambda2,b_star_lambda2),0,10))	

		############# save MCMC iterates #############
		if i>burnin && i%thin==0 
			Si_out[:,:,i_out] = Si_iter[:,1:T]
			gamma_out[:,:,i_out] = gamma_iter[:,1:T]
			if time_specific_alpha==false && unit_specific_alpha==false
				# for each iterate, a scalar
				alpha_out[i_out] = alpha_iter
			elseif time_specific_alpha==true && unit_specific_alpha==false
				# for each iterate, a vector in time
				alpha_out[:,i_out] = alpha_iter[1:T]
			elseif time_specific_alpha==false && unit_specific_alpha==true
				# for each iterate, a vector in units
				alpha_out[:,i_out] = alpha_iter
			elseif time_specific_alpha==true && unit_specific_alpha==true
				# for each iterate, a matrix
				alpha_out[:,:,i_out] = alpha_iter[:,1:T]
			end
			for t in 1:T
				for j in 1:n
					sigma2h_out[j,t,i_out] = sig2h_iter[Si_iter[j,t],t]
					muh_out[j,t,i_out] = muh_iter[Si_iter[j,t],t]
				end
			end
			eta1_out[:,i_out] = eta1_iter
			if lk_xPPM
				for t in 1:T
					beta_out[t,:,i_out] = beta_iter[t]
				end
			end
			theta_out[:,i_out] = theta_iter[1:T]
			tau2_out[:,i_out] = tau2_iter[1:T]
			phi0_out[i_out] = phi0_iter
			phi1_out[i_out] = phi1_iter
			lambda2_out[i_out] = lambda2_iter


			############# save fitted values and metrics #############
			for j in 1:n
				for t in 1:T
					muh_jt = muh_iter[Si_iter[j,t],t]
					sig2h_jt = sig2h_iter[Si_iter[j,t],t]
					X_lk_term = lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0.0

					if t==1
						llike[j,t,i_out] = logpdf(Normal(
							muh_jt + X_lk_term,
							sqrt(sig2h_jt)
							), Y[j,t])
						fitted[j,t,i_out] = muh_jt + X_lk_term
					else # t>1
						llike[j,t,i_out] = logpdf(Normal(
							muh_jt + eta1_iter[j]*Y[j,t-1] + X_lk_term,
							sqrt(sig2h_jt*(1-eta1_iter[j]^2))
							), Y[j,t])
						fitted[j,t,i_out] = muh_jt + eta1_iter[j]*Y[j,t-1] + X_lk_term
					end

					mean_likelhd[j,t] += exp(llike[j,t,i_out])
					mean_loglikelhd[j,t] += llike[j,t,i_out]
					CPO[j,t] += exp(-llike[j,t,i_out])
				end
			end

			i_out += 1
		end

	next!(progresso)

	end # for i in 1:draws

	println("\ndone!")
	t_end = now()
	println("Elapsed time: ", Dates.canonicalize(Dates.CompoundPeriod(t_end-t_start)))

	############# compute LPML and WAIC #############
	for j in 1:n
		for t in 1:T
			LPML += log(CPO[j,t])
		end
	end
	LPML -= n*T*log(nout) # scaling factor
	LPML = -LPML # fix sign

	println("LPML: ", LPML, " (the higher the better)")
	# println("LPML: $LPML - the ↑ the :)")
	
	# adjust mean variables
	mean_likelhd ./= nout
	mean_loglikelhd./= nout
	for j in 1:n
		for t in 1:T
			WAIC += 2*mean_loglikelhd[j,t] - log(mean_likelhd[j,t])
		end
	end
	WAIC *= -2
	println("WAIC: ", WAIC, " (the lower the better)")
	# println("WAIC: $WAIC - the ↓ the :)")
	
	if update_eta1 @printf "acceptance ratio eta1: %.2f%%\n" acceptance_ratio_eta1/(n*draws) *100 end
	if update_phi1 @printf "acceptance ratio phi1: %.2f%%" acceptance_ratio_phi1/draws*100 end
	println()

	if perform_diagnostics
		chn = Chains(
			hcat(lambda2_out,phi0_out,tau2_out',theta_out',eta1_out',alpha_out'),
			["lambda2","phi0",
			[string("tau2_t", i) for i in 1:T]...,
			[string("theta_t", i) for i in 1:T]...,
			[string("eta1_j", i) for i in 1:n]...,
			[string("alpha_t", i) for i in 1:T]...,
			]
		)
		ss = DataFrame(summarystats(chn))
		println("\nDiagnostics:")
		@show ss[!,[1,4,5,6,7]];
		if logging CSV.write(log_file,ss[!,[1,4,5,6,7]]) end
	end

	close(log_file)

	if simple_return
		return Si_out, LPML, WAIC
	else
		return Si_out, Int.(gamma_out), alpha_out, sigma2h_out, muh_out, include_eta1 ? eta1_out : NaN,
			lk_xPPM ? beta_out : NaN, theta_out, tau2_out, phi0_out, include_phi1 ? phi1_out : NaN, lambda2_out,
			fitted, llike, LPML, WAIC
	end

end

function close_log_file()
	close(log_file)
end