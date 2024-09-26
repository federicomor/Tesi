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

include("utils.jl")

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
										  # plus the std dev for the Metropolis update trough N(eta1_old,mhsig_eta1^2)
	beta_priors = missing,                # Prior parameters for beta ∼ 
	tau2_priors::Vector{Float64},         # Prior parameters for tau2 ∼ invGamma(a_tau, b_tau)
	phi0_priors::Vector{Float64},         # Prior parameters for phi0 ∼ N(m0, s0^2)
	phi1_priors::Float64,                 # Prior parameters for phi1 ∼ U(-1,1)
										  # so we just need the std dev of the Metropolis update trough N(phi1_old,mhsig_phi1^2)
	lambda2_priors::Vector{Float64},      # Prior parameters for lambda2 ∼ invGamma(a_lambda, b_lambda)
	alpha_priors::AbstractArray{Float64}, # Prior parameters for alpha ∼ Beta(a_alpha, b_alpha)
										  # but possibly that pair for each unit j, that's why the abstract array
	
	spatial_cohesion_idx = missing,       # cohesion choice
	sp_params = missing,                  # Parameters for spatial cohesion functions
	covariate_similarity_idx = missing,   # similarity choice
	cv_params = missing,                  # Parameters for covariates similarity functions

	draws::Float64,                       # Number of MCMC draws
	burnin::Float64,                      # Number of burn-in
	thin::Float64,                        # Thinning interval

	logging = false,                      # Wheter to save execution infos to log file
	seed::Float64,                        # Random seed for reproducibility
	simple_return = false                 # Return just the partition Si
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
	beta_update_threshold = round(Int64, burnin/2)

	if sPPM
		sp1::Vector{Float64} = copy(vec(sp_coords[:,1]))
		sp2::Vector{Float64} = copy(vec(sp_coords[:,2]))

		S=@MMatrix zeros(2, 2)
		if spatial_cohesion_idx==3 || spatial_cohesion_idx==4
			sp_params_real = [SVector{2}(sp_params[1]...), sp_params[2], sp_params[3], SMatrix{2,2}(sp_params[4]...)]
		else 
			sp_params_real = sp_params
		end
	end
	
	if update_eta1 acceptance_ratio_eta1 = 0 end
	if update_phi1 acceptance_ratio_phi1 = 0 end


	if logging
		log_file = open("log.txt", "w+")
		include("debug.jl")
		println("Logging to file: ", abspath("log.txt"))

		printlgln(replace(string(now()),"T" => "   "))
		debug("current seed = $seed")
		# to = TimerOutput()
	end

# try

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
		elseif covariate_similarity == 4 && length.(cv_params) != [1,1,1,1]
			@error "Wrong params for covariate similarity 4.\nExpected input form: [Real, Real, Real, Real]." _file=""
			return
		end
	end

	if lk_xPPM
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
	println("- using seed $seed -")
	println("fitting $(Int(draws)) total iterates (with burnin=$(Int(burnin)), thinning=$(Int(thin)))")
	println("thus producing $nout valid iterates in the end")
	println("\non n=$n subjects\nfor T=$T time instants")
	println("\nwith space? $sPPM")
	println("with covariates in the likelihood? $lk_xPPM", lk_xPPM ? " (p=$p_lk)" : "")
	println("with covariates in the clustering process? $cl_xPPM", cl_xPPM ? " (p=$p_cl)" : "")
	println("are there missing data in Y? $Y_has_NA")
	println()


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
	# debug(@showd missing_map)
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
	lPP = @MVector zeros(1)
	Si_red1 = zeros(Int64,n)
	s1n = zeros(n)
	s1o = zeros(n)
	s2n = zeros(n)
	s2o = zeros(n)
	nh_tmp = zeros(Int64,n)
	rho_tmp = zeros(Int64,n)


	############# allocate and initialize working variables #############
	Si_iter = ones(Int64,n,T_star) # label assignements for units j at time t
	Si_iter[:,end] .= 0
	gamma_iter = zeros(Bool,n,T_star)
	if !ismissing(initial_partition)
		relabel!(initial_partition,n,aux1_relab)
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
	lambda2_iter = rand(InverseGamma(lambda2_priors...))

	# hierarchy level 2
	theta_iter = zeros(T_star)
	theta_iter[1] = rand(Normal(phi0_iter,sqrt(lambda2_iter)))
	@inbounds for t in 2:T_star
		theta_iter[t] = rand(Normal((1-phi1_iter)*phi0_iter + phi1_iter*theta_iter[t-1],sqrt(lambda2_iter*(1-phi1_iter^2))))
	end
	tau2_iter = zeros(T_star)
	@inbounds for t in 1:T_star
		tau2_iter[t] = rand(InverseGamma(tau2_priors...))
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
		# debug(@showd beta_iter)
	end

	eta1_iter = zeros(n)
	if include_eta1
		@inbounds for j in 1:n
			eta1_iter[j] = rand(Uniform(-1,1))
		end
	end

	nh = zeros(Int,n,T_star) # numerosity of each cluster k at time t
	# dimension n since at most there are n clusters (all singletons)
	nh[1,:] .= n
	nclus_iter = ones(Int,T_star) # how many clusters there are at time t


	############# start MCMC algorithm #############
	println("Starting MCMC algorithm")

	t_start = now()
	progresso = Progress(round(Int64(draws)),
			showspeed=true,
			output=stdout, # default is stderr, which turns out in orange color on R
			dt=1, # every how many seconds update the feedback
			barlen=0 # no progress bar
			)

	@inbounds for i in 1:draws
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

				# printlgln("########### iteration $i ###########\n")
				# debug(@showd c_it Si_iter sig2h_iter muh_iter)

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

					# mu_prior = muh_iter[c_it,t] + eta1_iter[j]*Y[j,t-1] + Xlk_term_t
					# sig2_prior = sig2h_iter[c_it,t]*(1-aux1)
					Y[j,t] = rand(Normal(mu_post,sqrt(sig2_post)))
					# debug(@showd c_itp1 mu_prior mu_post sig2_prior sig2_post)

				else # t==T
					Y[j,t] = rand(Normal(
						muh_iter[c_it,t] + eta1_iter[j]*(ismissing(Y[j,t-1]) ? 0 : Y[j,t-1]) + Xlk_term_t,
						sqrt(sig2h_iter[c_it,t]*(1-aux1))
						))
				end
				# debug(@showd Y[j,t])
			end
		end


		for t in 1:T
			############# update gamma #############
			# @timeit to " gamma " begin # if logging uncomment this line, and the corresponding "end"
			for j in 1:n
				if t==1 
					gamma_iter[j,t] = 0
					# at the first time units get reallocated
				else
					# we want to find ρ_t^{R_t(-j)} ...
					indexes = findall_faster(jj -> jj != j && gamma_iter[jj, t] == 1, 1:n)
					Si_red = Si_iter[indexes, t]
					# Si_red1 = copy(Si_red)
					copy!(Si_red1, Si_red)
					# Si_red1 = Si_iter[indexes, t]
					push!(Si_red1, Si_iter[j,t]) # ... and ρ_t^R_t(+j)}

					# get also the reduced spatial info if sPPM model
					# @timeit to " sPPM 1 " begin # if logging uncomment this line, and the corresponding "end"
					if sPPM
						sp1_red = @view sp1[indexes]
						sp2_red = @view sp2[indexes]
						# sp1_red = sp1[indexes]
						# sp2_red = sp2[indexes]
					end
					# debug(@showd sp1 sp2 sp1_red sp2_red)
					# and the reduced covariates info if cl_xPPM model
					if cl_xPPM
						Xcl_covariates_red = @view Xcl_covariates[indexes,:,t]
						# debug(@showd size(Xcl_covariates_red))
						# debug(@showd Xcl_covariates_red)
						# debug(@showd Xcl_covariates)
					end
					# end # of the @timeit for sPPM

					# compute n_red's and nclus_red's and relabel
					n_red = length(Si_red) # = "n" relative to here, i.e. the sub-partition size
					n_red1 = length(Si_red1)
					# @assert n_red1 == n_red+1
					relabel!(Si_red,n_red,aux1_relab)
					relabel!(Si_red1,n_red1,aux1_relab)
					nclus_red = isempty(Si_red) ? 0 : maximum(Si_red) # = number of clusters
					nclus_red1 = maximum(Si_red1)

					# save the label of the current working-on unit j
					j_label = Si_red1[end]

					# compute also nh_red's
					# nh_red = zeros(Int64,nclus_red)
					# nh_red1 = zeros(Int64,nclus_red1)
					nh_red .= 0
					nh_red1 .= 0
					for jj in 1:n_red
						nh_red[Si_red[jj]] += 1 # = numerosities for each cluster label
						nh_red1[Si_red1[jj]] += 1
					end
					nh_red1[Si_red1[end]] += 1 # account for the last added unit j, not included in the above loop
					# debug(@showd Si_red Si_red1 j_label)

					# start computing weights
					# lg_weights = zeros(nclus_red+1)
					lg_weights .= 0
					# lCo = 0.0; lCn = 0.0 # log cohesions (for space) old and new
					# lC = @MVector zeros(2)
					lSo = 0.0; lSn = 0.0 # log similarities (for covariates) old and new

					# unit j can enter an existing cluster...
					for k in 1:nclus_red
						# @timeit to " sPPM 2 " begin # if logging uncomment this line, and the corresponding "end"

						# filter the covariates of the units of label k
						# aux_idxs = findall(jj -> Si_red[jj] == k, 1:n_red) # slow
						aux_idxs = findall_faster(jj -> Si_red[jj] == k, 1:n_red)
						# aux_idxs = findall(Si_red .== k) # fast
						if sPPM
							# filter the spatial coordinates of the units of label k
							# debug(@showd sp_idxs)
							# debug(@showd k Si_red)

							# pretty_log("sp_coords")
							# pretty_log("sp1")
							# pretty_log("sp2")
							# debug(@showd sp_idxs)
							# printlgln("\n")

							# "forse qui qualcosa da poter fare dovrebbe esserci" con le views?
							# o forse no per via della push, che modificherebbe la view
							# s1o = @view sp1_red[aux_idxs]
							# s2o = @view sp2_red[aux_idxs]
							# s1n = copy(s1o); push!(s1n, sp1[j])
							# s2n = copy(s2o); push!(s2n, sp2[j])

							copy!(s1o, sp1_red[aux_idxs])
							copy!(s2o, sp2_red[aux_idxs])
							# s1o_v = @view sp1_red[aux_idxs]
							# s2o_v = @view sp2_red[aux_idxs]
							copy!(s1n,s1o); push!(s1n, sp1[j])
							copy!(s2n,s2o); push!(s2n, sp2[j])
							# copy!(s1n,sp1_red[aux_idxs]); push!(s1n, sp1[j])
							# copy!(s2n,sp2_red[aux_idxs]); push!(s2n, sp2[j])

							# println(typeof(s1n))
							# resize!(s1n,length(s1o)+1)
							# copy!(s1n,s1o)
							# s1n[end] = sp1[j]
							# resize!(s2n,length(s2o)+1)
							# copy!(s2n,s1o)
							# s2n[end] = sp2[j]
							# copy!(s1n,s1o); push!(s1n, sp1[j])
							# copy!(s2n,s2o); push!(s2n, sp2[j])
							# lCo = spatial_cohesion(spatial_cohesion_idx, s1o, s2o, sp_params_real, true, M_dp, S)
							# lCn = spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S)
							# spatial_cohesion!(spatial_cohesion_idx, s1o_v, s2o_v, sp_params_real, true, M_dp, S,1,false,lC)
							spatial_cohesion!(spatial_cohesion_idx, s1o, s2o, sp_params_real, true, M_dp, S,1,false,lC)
							spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S,2,false,lC)
						end
						# debug(@showd lCn lCo)
						# debug(@showd cl_xPPM)

						# Xcl_covariates is a n*p*T matrix
						if cl_xPPM
							for p in 1:p_cl
								Xo = @view Xcl_covariates_red[aux_idxs,p]
								Xn = copy(Xo); push!(Xn,Xcl_covariates[j,p,t])
								lSo += covariate_similarity(covariate_similarity_idx, Xo, cv_params, lg=true)
								lSn += covariate_similarity(covariate_similarity_idx, Xn, cv_params, lg=true)
							end
						end
						# end # of the @timeit for sPPM
						# debug(@showd lCn lCo lSn lSo)
						# printlgln("\n")
						# lg_weights[k] = log(nh_red[k]) + lCn - lCo + lSn - lSo
						lg_weights[k] = log(nh_red[k]) + lC[2] - lC[1] + lSn - lSo
					end
					
					# ... or unit j can create a singleton
					# lCn = 0.0
					lSn = 0.0 
					# @timeit to " sPPM 3 " begin # if logging uncomment this line, and the corresponding "end"
					if sPPM
						# lCn = spatial_cohesion(spatial_cohesion_idx, [sp1[j]], [sp2[j]], sp_params_real, lg=true, M=M_dp, S=S)
						# lCn = spatial_cohesion(spatial_cohesion_idx, [sp1[j]], [sp2[j]], sp_params_real, true, M_dp, S)
						spatial_cohesion!(spatial_cohesion_idx, [sp1[j]], [sp2[j]], sp_params_real, true, M_dp, S,2,false,lC)
						# lCn = spatial_cohesion(spatial_cohesion_idx, SVector(sp1[j]), SVector(sp2[j]), sp_params_real, lg=true, M=M_dp, S=S)
					end
					if cl_xPPM
						for p in 1:p_cl
							lSn += covariate_similarity(covariate_similarity_idx, [Xcl_covariates[j,p,t]], cv_params, lg=true)
						end
					end
					# end # of the @timeit for sPPM
					# lg_weights[nclus_red+1] = log(M_dp) + lCn + lSn
					# lg_weights[nclus_red+1] = log_Mdp + lCn + lSn
					lg_weights[nclus_red+1] = log_Mdp + lC[2] + lSn

					# printlgln("before exp and normalization:")
					# debug(@showd lg_weights)

					# now use the weights towards sampling the new gamma_jt
					# max_ph = maximum(lg_weights)
					max_ph = maximum(@view lg_weights[1:(nclus_red+1)])
					# max_ph = maximum(lg_weights[1:(nclus_red+1)])

					sum_ph = 0.0

					# exponentiate...
					# for k in eachindex(lg_weights)
					for k in 1:(nclus_red+1)
						 # for numerical purposes we subract max_ph
						lg_weights[k] = exp(lg_weights[k] - max_ph)
						sum_ph += lg_weights[k]
					end
					# ... and normalize
					lg_weights ./= sum_ph
					# printlgln("after exp and normalization:")
					# debug(@showd lg_weights)

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
					# debug(@showd probh)


					# compatibility check for gamma transition
					if gamma_iter[j, t] == 0
						# we want to find ρ_(t-1)^{R_t(+j)} ...
						# indexes = findall(jj -> gamma_iter[jj, t]==1, 1:n) # slow
						indexes = findall_faster(jj -> gamma_iter[jj, t]==1, 1:n)
						# indexes = findall(gamma_iter[:, t] .== 1) # fast
						union!(indexes,j)
						Si_comp1 = @view Si_iter[indexes, t-1]
						Si_comp2 = @view Si_iter[indexes, t] # ... and ρ_t^R_t(+j)}

						# rho_comp = compatibility(Si_comp1, Si_comp2)
						rho_comp = compatibility(Si_comp1, Si_comp2,aux1_relab,aux2_relab)
						if rho_comp == 0
							probh = 0.0
						end
					end
					# sample the new gamma
					gt = rand(Bernoulli(probh))
					gamma_iter[j, t] = gt
					# printlgln("sampled gamma:")
					# debug(@showd probh gt gamma_iter[:,t])
				end
			end # for j in 1:n


			# end # of the @timeit for gamma
			############# update rho #############
			# @timeit to " rho " begin # if logging uncomment this line, and the corresponding "end"
			# we only update the partition for the units which can move (i.e. with gamma_jt=0)
			# movable_units = findall(j -> gamma_iter[j,t]==0, 1:n) # slow
			movable_units = findall_faster(j -> gamma_iter[j,t]==0, 1:n) 
			# movable_units = findall(gamma_iter[:,t] .== 0) # fast
		
			for j in movable_units
				# printlgln("working on unit j=$j at t=$t")
				# remove unit j from the cluster she is currently in
				# debug(@showd j t)
				# debug(@showd movable_units)

				if nh[Si_iter[j,t],t] > 1 # unit j does not belong to a singleton cluster
					nh[Si_iter[j,t],t] -= 1
					# no nclus_iter[t] change since j's cluster is still alive
				else # unit j belongs to a singleton cluster
					j_label = Si_iter[j,t]
					last_label = nclus_iter[t]
					# print(@showd last_label)

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
				# ph = zeros(nclus_iter[t]+1) 
				ph .= 0.0
				# rho_tmp = copy(Si_iter[:,t])
				copy!(rho_tmp, @view Si_iter[:,t])
				# rho_tmp = Si_iter[:,t]

				# compute nh_tmp (numerosities for each cluster label)
				# nh_tmp = copy(nh[:,t])
				copy!(nh_tmp, @view nh[:,t])
				# nh_tmp = zeros(Int,nclus_iter[t]+1)
				# for jj in setdiff(1:n,j)
					# nh_tmp[rho_tmp[jj]] += 1
				# end
				nclus_temp = 0
					
				# we now simulate the unit j to be assigned to one of the existing clusters...
				for k in 1:nclus_iter[t]
					rho_tmp[j] = k
					# indexes = findall(j -> gamma_iter[j,t+1]==1, 1:n) # slow
					indexes = findall_faster(j -> gamma_iter[j,t+1]==1, 1:n)
					# indexes = findall(gamma_iter[:,t+1] .== 1) # fast
					# we check the compatibility between ρ_t^{h=k,R_(t+1)} ...
					Si_comp1 = @view rho_tmp[indexes]
					Si_comp2 = @view Si_iter[indexes,t+1] # and ρ_(t+1)^{R_(t+1)}
					# rho_comp = compatibility(Si_comp1, Si_comp2)
					rho_comp = compatibility(Si_comp1, Si_comp2,aux1_relab,aux2_relab)
					# printlgln("Assigning to cluster k=$k :")
					# pretty_log("gamma_iter"); # pretty_log("Si_iter")
					# debug(@showd rho_tmp Si_comp1 Si_comp2 rho_comp)
					
					if rho_comp != 1
						ph[k] = log(0) # assignment to cluster k is not compatible
					else
						# update params for "rho_jt = k" simulation
						nh_tmp[k] += 1
						# nclus_temp = sum(nh_tmp .> 0) # = number of clusters
						nclus_temp = count(a->(a>0), nh_tmp) # similarly fast, not clear
						# count is a bit slower but does not allocate

						# lpp = 0.0
						lPP .= 0.
						for kk in 1:nclus_temp
							# @timeit to " sPPM 4 " begin # if logging uncomment this line, and the corresponding "end"
							# aux_idxs = findall(jj -> rho_tmp[jj]==kk, 1:n) # slow
							aux_idxs = findall_faster(jj -> rho_tmp[jj]==kk, 1:n)
							# aux_idxs = findall(rho_tmp .== kk) # fast
							if sPPM
								# s1n = @view sp1[aux_idxs]
								# s2n = @view sp2[aux_idxs]

								copy!(s1n, @view sp1[aux_idxs])
								copy!(s2n, @view sp2[aux_idxs])

								# lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params_real, lg=true, M=M_dp, S=S)
								# lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S)
								# spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S,1,true,lPP)
								spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S,1,true,lPP)
							end
							if cl_xPPM
								# debug(@showd cv_idxs)
								for p in 1:p_cl
									# debug(@showd p)
									Xn = @view Xcl_covariates[aux_idxs,p,t]
									# debug(@showd Xn)
									lPP[1] += covariate_similarity(covariate_similarity_idx, Xn, cv_params, lg=true)
								end
							end
							# end # of the @timeit for sPPM
							# lpp += log(M_dp) + lgamma(nh_tmp[kk])
							# lpp += log_Mdp + lgamma(nh_tmp[kk])
							lPP[1] += log_Mdp + lgamma(nh_tmp[kk])
							# lpp += log(M_dp) + lgamma(length(indexes)) # same
						end

						### debug case
						# ph[k] = 0.5
						### real case
						if t==1
							ph[k] = loglikelihood(Normal(
								muh_iter[k,t] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
								sqrt(sig2h_iter[k,t])),
								# Y[j,t]) + lpp
								Y[j,t]) + lPP[1]
						else
							ph[k] = loglikelihood(Normal(
								muh_iter[k,t] + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
								sqrt(sig2h_iter[k,t]*(1-eta1_iter[j]^2))),
								# Y[j,t]) + lpp
								Y[j,t]) + lPP[1]
						end

						# restore params after "rho_jt = k" simulation
						nh_tmp[k] -= 1
					end
				end
					
				# print(@showd ph, "before k+1")
				# ... plus the case of being assigned to a new (singleton for now) cluster 
				k = nclus_iter[t]+1
				rho_tmp[j] = k
				# declare (for later scope accessibility) the new params here
				muh_draw = 0.0; sig2h_draw = 0.0

				# indexes = findall(j -> gamma_iter[j,t+1]==1, 1:n) # slow
				indexes = findall_faster(j -> gamma_iter[j,t+1]==1, 1:n) 
				# indexes = findall(gamma_iter[:,t+1] .== 1)
				Si_comp1 = @view rho_tmp[indexes]
				Si_comp2 = @view Si_iter[indexes,t+1]
				# rho_comp = compatibility(Si_comp1, Si_comp2)
				rho_comp = compatibility(Si_comp1, Si_comp2,aux1_relab,aux2_relab)
				# printlgln("Assigning to NEW SINGLETON cluster (k=$k) :")
				# pretty_log("gamma_iter"); # pretty_log("Si_iter")
				# debug(@showd rho_tmp Si_comp1 Si_comp2 rho_comp)

				if rho_comp != 1
					ph[k] = log(0) # assignment to a new cluster is not compatible
				else
					# sample new params for this new cluster
					muh_draw = rand(Normal(theta_iter[t], sqrt(tau2_iter[t])))
					sig2h_draw = rand(InverseGamma(2,2))

					# update params for "rho_jt = k" simulation
					nh_tmp[k] += 1
					# nclus_temp = sum(nh_tmp .> 0)
					nclus_temp = count(a->(a>0), nh_tmp) # similarly fast, not clear
					# count is a bit slower but does not allocate

					# lpp = 0.0
					lPP .= 0.
					for kk in 1:nclus_temp
						# @timeit to " sPPM 5 " begin # if logging uncomment this line, and the corresponding "end"
						# aux_idxs = findall(jj -> rho_tmp[jj]==kk, 1:n) # slow
						aux_idxs = findall_faster(jj -> rho_tmp[jj]==kk, 1:n)
						# aux_idxs = findall(rho_tmp .== kk) # fast
						if sPPM
							# s1n = @view sp1[aux_idxs]
							# s2n = @view sp2[aux_idxs]

							copy!(s1n, @view sp1[aux_idxs])
							copy!(s2n, @view sp2[aux_idxs])
							# lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params_real, lg=true, M=M_dp, S=S)
							# lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S)
							# spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S,1,true,lPP)
							spatial_cohesion!(spatial_cohesion_idx, s1n, s2n, sp_params_real, true, M_dp, S,1,true,lPP)
						end
						if cl_xPPM
							for p in 1:p_cl
								Xn = @view Xcl_covariates[aux_idxs,p,t]
								lPP[1] += covariate_similarity(covariate_similarity_idx, Xn, cv_params, lg=true)
							end
						end
						# end # of the @timeit for sPPM
						# lpp += log(M_dp) + lgamma(nh_tmp[kk])
						# lpp += log_Mdp + lgamma(nh_tmp[kk])
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
							# Y[j,t]) + lpp
							Y[j,t]) + lPP[1]
					else
						ph[k] = loglikelihood(Normal(
							muh_draw + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_draw*(1-eta1_iter[j]^2))),
							# Y[j,t]) + lpp
							Y[j,t]) + lPP[1]
					end

					# restore params after "rho_jt = k" simulation
					nh_tmp[k] -= 1
				end

				# printlgln("Before exp and normalization:")
				# debug(@showd ph)
				# now exponentiate the weights...
				# max_ph = maximum(ph)
				max_ph = maximum(@view ph[1:(nclus_iter[t]+1)])
				# max_ph = maximum(ph[1:(nclus_iter[t]+1)])
				sum_ph = 0.0
				for k in 1:(nclus_iter[t]+1)
					# for numerical purposes we subract max_ph
					ph[k] = exp(ph[k] - max_ph)
					sum_ph += ph[k]
				end
				# ... and normalize them
				ph ./= sum_ph
				
				# printlgln("After exp and normalization:")
				# debug(@showd ph)

				# now sample the new label Si_iter[j,t]
				u = rand(Uniform(0,1))
				cph = cumsum(ph)
				cph[end] = 1.0 # fix numerical problems of having sums like 0.999999etc
				new_label = 0
				for k in 1:length(ph)
					if u <= cph[k]
						new_label = k
						break
					end
				end
				# debug(@showd cph u)
				# debug(@showd new_label j)
				# printlgln("\n")
				
				# old if for testing, now we can remove to optimize/speed up the computation
				# if all(isinf.(ph) .|| isnan.(ph))
				# 	@warn "No valid weight in ph vector when updating for rho.\n"
				# 	close(log_file); return
				# end

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
					# rho_tmp[j] = cl_new
					muh_iter[cl_new, t] = muh_draw
					sig2h_iter[cl_new, t] = sig2h_draw
				end

				# now we need to relabel after the possible mess created by the sampling
				# eg:       (before)            (after)
				#     units j 2 3 4 5 ->  units j 2 3 4 5
				#    labels 1 1 1 2 2    labels 3 1 1 2 2
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
				# discard the zeros at the end of the auxiliary vectors nh_reorde and old_lab
				# muh_iter_copy = copy(muh_iter)
				copy!(muh_iter_copy, muh_iter) # copy!(dst,src)
				# sig2h_iter_copy = copy(sig2h_iter)
				copy!(sig2h_iter_copy, sig2h_iter) # copy!(dst,src)
				len = findlast(x -> x != 0, nh_reorder)
				for k in 1:nclus_iter[t]
					muh_iter[k,t] = muh_iter_copy[old_lab[k],t]
					sig2h_iter[k,t] = sig2h_iter_copy[old_lab[k],t]
					nh[k,t] = nh_reorder[k]
				end
				# debug(@showd Si_iter[:,t])

			end # for j in movable_units

			# end # of the @timeit for rho
			############# update muh #############
			# @timeit to " muh " begin # if logging uncomment this line, and the corresponding "end"
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
			
			# end # of the @timeit for muh
			############# update sigma2h #############
			# @timeit to " sigma2h " begin # if logging uncomment this line, and the corresponding "end"
			if t==1
				for k in 1:nclus_iter[t]
					# a_star = sig2h_priors[1] + sum(Si_iter[:,t] .== k)/2
					a_star = sig2h_priors[1] + nh[k,t]/2 # should be the same
					sum_Y = 0.0
					# S_kt = findall(j -> Si_iter[j,t] == k, 1:n) # slow
					S_kt = findall_faster(j -> Si_iter[j,t] == k, 1:n)
					# S_kt = findall(Si_iter[:,t] .== k) # fast
					# if lk_xPPM
					# 	for j in S_kt
					# 		sum_Y += (Y[j,t] - muh_iter[k,t] - dot(Xlk_covariates[j,:,t], beta_iter[t]))^2
					# 	end
					# else 
					# 	for j in S_kt
					# 		sum_Y += (Y[j,t] - muh_iter[k,t])^2
					# 	end
					# end
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
					# S_kt = findall(j -> Si_iter[j,t] == k, 1:n) # slow
					S_kt = findall_faster(j -> Si_iter[j,t] == k, 1:n)
					# S_kt = findall(Si_iter[:,t] .== k) # fast
					# if lk_xPPM
					# 	for j in S_kt
					# 		sum_Y += (Y[j,t] - muh_iter[k,t] - eta1_iter[j]*Y[j,t-1] - dot(Xlk_covariates[j,t], beta_iter[t]))^2
					# 	end
					# else 
					# 	for j in S_kt
					# 		sum_Y += (Y[j,t] - muh_iter[k,t] - eta1_iter[j]*Y[j,t-1])^2
					# 	end
					# end
					for j in S_kt
						sum_Y += (Y[j,t] - muh_iter[k,t] - eta1_iter[j]*Y[j,t-1] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0.0))^2
					end

					b_star = sig2h_priors[2] + sum_Y/2
					sig2h_iter[k,t] = rand(InverseGamma(a_star, b_star))
				end
			end
			# end # of the @timeit for sigma2h
			############# update beta #############
			# @timeit to " beta " begin # if logging uncomment this line, and the corresponding "end"
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
						
			# end # of the @timeit for beta
			############# update theta #############
			# @timeit to " theta " begin # if logging uncomment this line, and the corresponding "end"
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
			
			# end # of the @timeit for theta
			############# update tau2 #############
			# @timeit to " tau2 " begin # if logging uncomment this line, and the corresponding "end"
			kt = nclus_iter[t]
			aux1 = 0.0
			for k in 1:kt
				aux1 += (muh_iter[k,t] - theta_iter[t])^2 
			end
			a_star = tau2_priors[1] + kt/2
			b_star = tau2_priors[2] + aux1/2
			tau2_iter[t] = rand(InverseGamma(a_star, b_star))
			# end # of the @timeit for tau2

		end # for t in 1:T
		
		############# update eta1 #############
		# @timeit to " eta1 " begin # if logging uncomment this line, and the corresponding "end"
		# eta1_priors[2] = sqrt(eta1_priors[2]) # from variance to std dev
		# no, the input argument is already the std dev
		if update_eta1
			for j in 1:n
				eta1_old = eta1_iter[j]
				eta1_new = rand(Normal(eta1_old,eta1_priors[2])) # proposal value

				if (-1 <= eta1_new <= 1)
					ll_old::Float64 = 0.0
					ll_new::Float64 = 0.0
					for t in 2:T
						# likelihood contribution
						ll_old += loglikelihood(Normal(
							muh_iter[Si_iter[j,t],t] + eta1_old*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_iter[j,t]*(1-eta1_old^2))
							), Y[j,t]) 
						ll_new += loglikelihood(Normal(
							muh_iter[Si_iter[j,t],t] + eta1_new*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_iter[j,t]*(1-eta1_new^2))
							), Y[j,t]) 
					end
					logit_old = logit(1/2*(eta1_old+1)) 
					logit_new = logit(1/2*(eta1_new+1)) 

					# prior contribution
					ll_old += -log(2*eta1_priors[1]) -1/eta1_priors[1]*abs(logit_old)
					ll_new += -log(2*eta1_priors[1]) -1/eta1_priors[1]*abs(logit_new)

					ll_ratio = min(ll_new-ll_old, 0)
					u = rand(Uniform(0,1))
					if (ll_ratio > log(u))
						eta1_iter[j] = eta1_new # accept the candidate
						acceptance_ratio_eta1 += 1
					end
				end
			end
		end

		# end # of the @timeit for eta1
		############# update alpha #############
		# @timeit to " alpha " begin # if logging uncomment this line, and the corresponding "end"
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

		# end # of the @timeit for alpha
		############# update phi0 #############
		# @timeit to " phi0 " begin # if logging uncomment this line, and the corresponding "end"
		aux1 = 1/lambda2_iter
		aux2 = 0.0
		# i found that looping on t rather than using sum(... for t in ...) seems faster
		for t in 2:T
			aux2 += theta_iter[t] - phi1_iter*theta_iter[t-1]
		end
		sig2_post = 1 / ( 1/phi0_priors[2] + aux1 * (1 + (T-1)*(1-phi1_iter)/(1+phi1_iter)) )
		mu_post = sig2_post * ( phi0_priors[1]/phi0_priors[2] + theta_iter[1]*aux1 + aux1/(1+phi1_iter)*aux2 )
		phi0_iter = rand(Normal(mu_post, sqrt(sig2_post)))
		
		# end # of the @timeit for phi0
		############# update phi1 #############
		# @timeit to " phi1 " begin # if logging uncomment this line, and the corresponding "end"
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

				ll_ratio = min(ll_new-ll_old, 0)
				u = rand(Uniform(0,1))
				if (ll_ratio > log(u))
					phi1_iter = phi1_new # accept the candidate
					acceptance_ratio_phi1 += 1
				end
			end
		end

		# end # of the @timeit for phi1
		############# update lambda2 #############
		# @timeit to " lambda2 " begin # if logging uncomment this line, and the corresponding "end"
		aux1 = 0.0
		# I found that looping on t rather than using sum(... for t in ...) seems faster
		for t in 2:T
			aux1 += (theta_iter[t] - (1-phi1_iter)*phi0_iter - phi1_iter*theta_iter[t-1])^2
		end
		a_star = lambda2_priors[1] + T/2
		b_star = lambda2_priors[2] + ((theta_iter[1] - phi0_iter)^2 + aux1) / 2
		lambda2_iter = rand(InverseGamma(a_star,b_star))	

		# end # of the @timeit for lambda2
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
			# debug(@showd fitted mean_likelhd mean_loglikelhd CPO)

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

	# println("LPML: $LPML (the higher the better)")
	println("LPML: ", LPML, " (the higher the better)") # not interpolating with $ is faster
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
	# println("WAIC: $WAIC (the lower the better)")
	println("WAIC: ", WAIC, " (the lower the better)") # not interpolating with $ is faster
	# println("WAIC: $WAIC - the ↓ the :)")

	# println()
	# println("acceptance ratio for eta1: ", acceptance_ratio_eta1/(n*draws) *100, "%")
	# println("acceptance ratio for phi1: ", acceptance_ratio_phi1/draws*100, "%")
	
	if update_eta1 @printf "acceptance ratio eta1: %.2f%%\n" acceptance_ratio_eta1/(n*draws) *100 end
	if update_phi1 @printf "acceptance ratio phi1: %.2f%%" acceptance_ratio_phi1/draws*100 end
	println()

	if logging
		# debug(@showd to)
		close(log_file)
	end

	if simple_return
		return Si_out, LPML, WAIC
	else
		return Si_out, Int.(gamma_out), alpha_out, sigma2h_out, muh_out, include_eta1 ? eta1_out : NaN,
			lk_xPPM ? beta_out : NaN, theta_out, tau2_out, phi0_out, include_phi1 ? phi1_out : NaN, lambda2_out,
			fitted, llike, LPML, WAIC
	end

# catch e
# 	println(e)
# 	close(log_file)
# end

# global_logger(ConsoleLogger())	
end

function close_log_file()
	close(log_file)
end