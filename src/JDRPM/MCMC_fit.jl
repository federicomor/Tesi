using Distributions
using Statistics
using LinearAlgebra
using Random
using Logging
using Dates
using TimerOutputs
using ProgressMeter
using StaticArrays

log_file = open("log.txt", "w+")
include("debug.jl")
include("utils.jl")

function MCMC_fit(;
	Y::Matrix{Float64},                   # n*T matrix, the observed values
	sp_coords = missing,                  # n*2 matrix, the spatial coordinates
	Xlk_covariates = missing,             # n*p*T matrix, the covariates to include in the likelihood
	Xcl_covariates = missing,             # n*p*T matrix, the covariates to include in the clustering process

	M_dp::Float64,                        # Dirichlet mass parameter
	initial_partition = missing,          # Initial partition (if provided)
	# actually not implemented yet the version with the initial partition (which however i think is useless)
	# però vabe dovrebbe essere facile da includere

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

	draws::Float64,                       # Number of MCMC draws
	burnin::Float64,                      # Number of burn-in sa
	thin::Float64,                        # Thinning interval

	logging::Bool,                        # Wheter to save execution infos to log file
	seed::Float64                         # Random seed for reproducibility
	)
	Random.seed!(round(Int64,seed))

	if logging
		println("Logging to file: ", abspath("log.txt"))

		printlgln(replace(string(now()),"T" => "   "))
		debug("current seed = $seed")
	end
	to = TimerOutput()

# try
	############# check some stuff #############
	if spatial_cohesion_idx == 1 && !(sp_params isa Real) 
		@error "Wrong params for spatial cohesion 1.\nExpected input form: (Real).\nReceived: $(typeof(sp_params))." _file=""
		return
	elseif spatial_cohesion_idx == 2 && !(sp_params isa Real)
		@error "Wrong params for spatial cohesion 2.\nExpected input form: (Real).\nReceived: $(typeof(sp_params))." _file=""
		return
	elseif spatial_cohesion_idx == 3 && !(sp_params isa Vector && length.(sp_params) == [2,1,1,4])
		@error "Wrong params for spatial cohesion 3.\nExpected input form: (1x2 Vector, Real, Real, 2x2 Matrix).\nReceived: $(typeof(sp_params))." _file=""
		return
	elseif spatial_cohesion_idx == 4 && !(sp_params isa Vector && length.(sp_params) == [2,1,1,4])
		@error "Wrong params for spatial cohesion 4.\nExpected input form: (1x2 Vector, Real, Real, 2x2 Matrix).\nReceived: $(typeof(sp_params))." _file=""
		return
	elseif spatial_cohesion_idx == 5 && !(sp_params isa Real)
		@error "Wrong params for spatial cohesion 5.\nExpected input form: (Real).\nReceived: $(typeof(sp_params))." _file=""
		return
	elseif spatial_cohesion_idx == 6 && !(sp_params isa Real)
		@error "Wrong params for spatial cohesion 6.\nExpected input form: (Real).\nReceived: $(typeof(sp_params))." _file=""
		return
	end

	if spatial_cohesion_idx == 3 || spatial_cohesion_idx == 4
		if !issymmetric(sp_params[4])
			@error "Matrix Psi of the spatial parameters must be symmetric." _file=""
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
		@error "Be coherent! I cant have update eta1 and not including it.\nCheck what you assigned to update_eta1 and include_eta1." _file=""
		return
	end
	if update_phi1==true && include_phi1==false
		@error "Be coherent! I cant have update phi1 and not including it.\nCheck what you assigned to update_phi1 and include_phi1." _file=""
		return
	end


	############# define auxiliary variables #############
	sPPM = !ismissing(sp_coords)
	cl_xPPM = !ismissing(Xcl_covariates)
	lk_xPPM = !ismissing(Xlk_covariates)
	n, T = size(Y)
	T_star = T+1
	p = lk_xPPM ? size(Xlk_covariates)[2] : 0

	if (draws-burnin)%thin != 0
		@error "Please define draws, thin and burnin in an integer division-friendly way.\nI.e., such that (draws-burnin) % thin = 0." _file=""
		return
	end
	nout = round(Int64, (draws - burnin)/(thin))
	if nout < 0 
		@error "Wrong iterations parameters." _file=""
		return
	end

	if sPPM
		sp1 = copy(vec(sp_coords[:,1]))
		sp2 = copy(vec(sp_coords[:,2]))
		sp_original = copy(sp_coords)
	end

	if lk_xPPM && ismissing(beta_priors)
		@error "Cannot use covariates in the likelihood if beta_priors is not defined." _file=""
		return
	end

	############# send feedback #############
	println("- using seed $seed -")
	println("fitting $(Int(draws)) total iterates (burnin=$(Int(burnin)), thinning=$(Int(thin)))")
	println("thus producing $nout valid iterates in the end")
	println("\non n=$n subjects\nfor T=$T time instants")
	println("with space? $sPPM")
	println("with covariates in the likelihood? $lk_xPPM", lk_xPPM ? " (p=$p)" : "")
	println("with covariates in the clustering process? $cl_xPPM")
	println()


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
		beta_out = zeros(T,p,nout)
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

	############# allocate and initialize working variables #############
	Si_iter = ones(Int64,n,T_star) # label assignements for units j at time t
	Si_iter[:,end] .= 0
	gamma_iter = zeros(Bool,n,T_star)
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
	for t in 2:T_star
		theta_iter[t] = rand(Normal((1-phi1_iter)*phi0_iter + phi1_iter*theta_iter[t-1],sqrt(lambda2_iter*(1-phi1_iter^2))))
	end
	tau2_iter = zeros(T_star)
	for t in 1:T_star
		tau2_iter[t] = rand(InverseGamma(tau2_priors...))
	end

	# hierarchy level 1
	sig2h_iter = ones(n,T_star)
	muh_iter = zeros(n,T_star)
	if lk_xPPM
		beta_iter = Vector{Vector{Float64}}(undef,T_star)
		beta0 = beta_priors[1:end-1]
		s2_beta = beta_priors[end]
		for t in 1:T_star
			beta_iter[t] = rand(MvNormal(beta0, s2_beta*I(p)))
		end
		# debug(@showd beta0 s2_beta)
	end

	eta1_iter = zeros(n)
	if include_eta1
		for j in 1:n
			eta1_iter[j] = rand(Uniform(-1,1))
		end
	end


	nh = zeros(Int,n,T_star) # numerosity of each cluster k at time t
	# dimension n since at most there are n clusters (all singletons)
	nh[1,:] .= n
	nclus_iter = ones(Int,T_star) # how many clusters there are at time t


	############# pre-allocate auxiliary working variables #############
	n_red = 0
	n_red1 = 0
	nclus_red = 0
	nclus_red1 = 0
	j_label = 0

	############# testing (for now) on random values #############
	# Si_iter = rand(collect(1:n),n,T_star)
	# gamma_iter = rand((0,1),n,T_star)
	

	############# some more logging #############
	function pretty_log(str)
		if str=="Si_iter" println(log_file,"Si_iter\n",tostr(Si_iter)); return; end
		if str=="gamma_iter" println(log_file,"gamma_iter\n",tostr(gamma_iter)); return; end
		if str=="nh" println(log_file,"nh\n",tostr(nh)); return; end
		if str=="nclus_iter" println(log_file,"nclus_iter\n",tostr(nclus_iter)); return; end
		if str=="muh_iter" println(log_file,"muh_iter\n",tostr(muh_iter)); return; end
		if str=="sig2h_iter" println(log_file,"sig2h_iter\n",tostr(sig2h_iter)); return; end
		if str=="alpha_iter" println(log_file,"alpha_iter\n",tostr(alpha_iter)); return; end
		if str=="theta_iter" println(log_file,"theta_iter\n",tostr(theta_iter)); return; end
		if str=="tau2_iter" println(log_file,"tau2_iter\n",tostr(tau2_iter)); return; end
		if str=="phi0_iter" println(log_file,"phi0_iter\n",tostr(phi0_iter)); return; end
		if str=="phi1_iter" println(log_file,"phi1_iter\n",tostr(phi1_iter)); return; end
		if str=="lambda2_iter" println(log_file,"lambda2_iter\n",tostr(lambda2_iter)); return; end
		if str=="sp_coords" println(log_file,"sp_coords\n",tostr(sp_coords)); return; end
		if str=="sp1" println(log_file,"sp1\n",tostr(sp1)); return; end
		if str=="sp2" println(log_file,"sp2\n",tostr(sp2)); return; end
		if str=="beta_iter" println(log_file,"beta_iter\n",tostr(beta_iter)); return; end
	end


	############# start MCMC algorithm #############
	println("Starting MCMC algorithm")
	sleep(1.0) # to let all the prints be printed
	# println("loading...\r")
	# sleep(1.0) # to let all the prints be printed

	t_start = now()
	progresso = Progress(round(Int64(draws)),
			showspeed=true,
			output=stdout, # default is stderr, which turns out in orange on R
			dt=1, # every how many seconds update the feedback
			barlen=0 # no progress bar
			)

	for i in 1:draws
		# print("iteration $i of $draws\r") # use only this when all finished, for a shorter feedback
		# there is also the ProgressMeter option, maybe is cooler

		# debug("\n▶ iteration $i")

		# pretty_log("Si_iter")
		# pretty_log("gamma_iter")
		# pretty_log("muh_iter")
		# pretty_log("sig2h_iter")
		# pretty_log("alpha_iter")

		for t in 1:T
			# debug("► time $t")

			############# update gamma #############
			# @timeit to " gamma " begin # if logging uncomment this line, and the corresponding "end"
			for j in 1:n
				# debug(title*"[update gamma]")
				# debug("▸ subject $j")
				if t==1 
					gamma_iter[j,t] = 0
					# debug("we are at time t=1 so nothing to do")
				else
					# we want to find ρ_t^{R_t(-j)} ...
					indexes = findall(jj -> jj != j && gamma_iter[jj, t] == 1, 1:n)
					Si_red = Si_iter[indexes, t]
					Si_red1 = copy(Si_red)
					push!(Si_red1, Si_iter[j,t]) # ... and ρ_t^R_t(+j)}

					# get also the reduced spatial info if sPPM model
					# @timeit to " sPPM 1 " begin # if logging uncomment this line, and the corresponding "end"
					if sPPM
						sp1_red = @view sp1[indexes]
						sp2_red = @view sp2[indexes]
					end
					# end # of the @timeit for sPPM

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
					nh_red = zeros(Int64,nclus_red)
					nh_red1 = zeros(Int64,nclus_red1)
					for jj in 1:n_red
						nh_red[Si_red[jj]] += 1 # = numerosities for each cluster label
						nh_red1[Si_red1[jj]] += 1
					end
					nh_red1[Si_red1[end]] += 1 # account for the last added unit j
					# debug(@showd Si_red Si_red1 j_label)

					# start computing weights
					lg_weights = zeros(nclus_red+1)
					lCo = 0.0; lCn = 0.0 # log cohesions (so for space) old and new
					lSo = 0.0; lSn = 0.0 # log similarities(so for covariates) old and new
					# TOWARDS COVARIATE DEVELOPMENT
					# unit j can enter an existing cluster...
					for k in 1:nclus_red
						# @timeit to " sPPM 2 " begin # if logging uncomment this line, and the corresponding "end"
						if sPPM
							# filter the spatial coordinates of the units of label k
							sp_idxs = findall(jj -> Si_red[jj] == k, 1:n_red)

							# pretty_log("sp_coords")
							# pretty_log("sp1")
							# pretty_log("sp2")
							# debug(@showd sp_idxs)
							# printlgln("\n")

							# forse qui qualcosa da poter fare dovrebbe esserci
							s1o = sp1_red[sp_idxs]
							s2o = sp2_red[sp_idxs]
							s1n = copy(s1o); push!(s1n, sp1[j])
							s2n = copy(s2o); push!(s2n, sp2[j])
							lCo = spatial_cohesion(spatial_cohesion_idx, s1o, s2o, sp_params, lg=true, M=M_dp)
							lCn = spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
						end
						# end # of the @timeit for sPPM
						lg_weights[k] = log(nh_red[k]) + lCn - lCo
					end
					
					# ... or unit j can create a singleton
					lCn = 0.0
					lSn = 0.0 # TOWARDS COVARIATE DEVELOPMENT
					# @timeit to " sPPM 3 " begin # if logging uncomment this line, and the corresponding "end"
					if sPPM
						lCn = spatial_cohesion(spatial_cohesion_idx, [sp1[j]], [sp2[j]], sp_params, lg=true, M=M_dp)
					end
					# end # of the @timeit for sPPM
					lg_weights[nclus_red+1] = log(M_dp) + lCn

					# printlgln("before exp and normalization:")
					# debug(@showd lg_weights)

					# now use the weights towards sampling the new gamma_jt
					max_ph = maximum(lg_weights)
					sum_ph = 0.0

					# exponentiate...
					for k in eachindex(lg_weights)
						 # for numerical purposes we subract max_ph
						lg_weights[k] = exp(lg_weights[k] - max_ph)
						sum_ph += lg_weights[k]
					end
					# ... and normalize
					lg_weights ./= sum_ph
					# printlgln("after exp and normalization:")
					# debug(@showd lg_weights)

					# compute probh
					probh = 0.0
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
						indexes = findall(jj -> gamma_iter[jj, t]==1, 1:n)
						union!(indexes,j)
						Si_comp1 = Si_iter[indexes, t-1]
						Si_comp2 = Si_iter[indexes, t] # ... and ρ_t^R_t(+j)}

						rho_comp = compatibility(Si_comp1, Si_comp2)
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
			# debug(title*"[update rho]")
			# we only update the partition for the units which can move (i.e. with gamma_jt=0)
			movable_units = findall(j -> gamma_iter[j,t]==0, 1:n)
			# debug(@showd movable_units Si_iter[:,t])
		
			for j in movable_units
				# remove unit j from the cluster she is currently in
				# debug(@showd j t)

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
				ph = zeros(nclus_iter[t]+1) 
				rho_tmp = copy(Si_iter[:,t])

				# compute nh_tmp (numerosities for each cluster label)
				nh_tmp = copy(nh[:,t])
				# nh_tmp = zeros(Int,nclus_iter[t]+1)
				# for jj in setdiff(1:n,j)
					# nh_tmp[rho_tmp[jj]] += 1
				# end
				nclus_temp = 0
					
				# we now simulate the unit j to be assigned to one of the existing clusters...
				for k in 1:nclus_iter[t]
					rho_tmp[j] = k
					indexes = findall(j -> gamma_iter[j,t+1]==1, 1:n)
					# we check the compatibility between ρ_t^{h=k,R_(t+1)} ...
					Si_comp1 = rho_tmp[indexes]
					Si_comp2 = Si_iter[indexes,t+1] # and ρ_(t+1)^{R_(t+1)}
					rho_comp = compatibility(Si_comp1, Si_comp2)
					# printlgln("Assigning to cluster k=$k :")
					# pretty_log("gamma_iter"); # pretty_log("Si_iter")
					# debug(@showd rho_tmp Si_comp1 Si_comp2 rho_comp)
					
					if rho_comp != 1
						ph[k] = log(0) # assignment to cluster k is not compatible
					else
						# update params for "rho_jt = k" simulation
						nh_tmp[k] += 1
						nclus_temp = sum(nh_tmp .> 0) # = number of clusters

						lpp = 0.0
						for kk in 1:nclus_temp
							# @timeit to " sPPM 4 " begin # if logging uncomment this line, and the corresponding "end"
							if sPPM
								indexes = findall(jj -> rho_tmp[jj]==kk, 1:n)
								s1n = @view sp1[indexes]
								s2n = @view sp2[indexes]
								lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
							end
							if cl_xPPM
								# TOWARDS COVARIATE DEVELOPMENT
							end
							# end # of the @timeit for sPPM
							lpp += log(M_dp) + lgamma(nh_tmp[kk])
							# lpp += log(M_dp) + lgamma(length(indexes)) # same
						end

						### debug case
						# ph[k] = 0.5
						### real case
						if t==1
							ph[k] = loglikelihood(Normal(
								muh_iter[k,t] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
								sqrt(sig2h_iter[k,t])),
								Y[j,t]) + lpp
						else
							ph[k] = loglikelihood(Normal(
								muh_iter[k,t] + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
								sqrt(sig2h_iter[k,t]*(1-eta1_iter[j]^2))),
								Y[j,t]) + lpp
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

				indexes = findall(j -> gamma_iter[j,t+1]==1, 1:n)
				Si_comp1 = rho_tmp[indexes]
				Si_comp2 = Si_iter[indexes,t+1]
				rho_comp = compatibility(Si_comp1, Si_comp2)
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
					nclus_temp = sum(nh_tmp .> 0)

					lpp = 0.0
					for kk in 1:nclus_temp
						# @timeit to " sPPM 5 " begin # if logging uncomment this line, and the corresponding "end"
						if sPPM
							indexes = findall(jj -> rho_tmp[jj]==kk, 1:n)
							s1n = @view sp1[indexes]
							s2n = @view sp2[indexes]
							lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
						end
						if cl_xPPM
							# TOWARDS COVARIATE DEVELOPMENT
						end
						# end # of the @timeit for sPPM
						lpp += log(M_dp) + lgamma(nh_tmp[kk])
						# lpp += log(M_dp) + lgamma(length(indexes)) # same
					end
					### debug case
					# ph[k] = 0.5
					### real case
					if t==1
						ph[k] = loglikelihood(Normal(
							muh_draw + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_draw)),
							Y[j,t]) + lpp
					else
						ph[k] = loglikelihood(Normal(
							muh_draw + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							sqrt(sig2h_draw*(1-eta1_iter[j]^2))),
							Y[j,t]) + lpp
					end

					# restore params after "rho_jt = k" simulation
					nh_tmp[k] -= 1
				end

				# printlgln("Before exp and normalization:")
				# debug(@showd ph)
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

				Si_relab, nh_reorder, old_lab = relabel_full(Si_tmp,n)				
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
				muh_iter_copy = copy(muh_iter)
				sig2h_iter_copy = copy(sig2h_iter)
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
							sum_Y += Y[j,t] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0)
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
							sum_Y += (Y[j,t] - eta1_iter[j]*Y[j,t-1] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0)) * aux1
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
					S_kt = findall(j -> Si_iter[j,t] == k, 1:n)
					# if lk_xPPM
					# 	for j in S_kt
					# 		sum_Y += (Y[j,t] - muh_iter[k,t] - dot(Xlk_covariates[j,t], beta_iter[t]))^2
					# 	end
					# else 
					# 	for j in S_kt
					# 		sum_Y += (Y[j,t] - muh_iter[k,t])^2
					# 	end
					# end
					for j in S_kt
						sum_Y += (Y[j,t] - muh_iter[k,t] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0))^2
					end

					b_star = sig2h_priors[2] + sum_Y/2
					sig2h_iter[k,t] = rand(InverseGamma(a_star, b_star))
				end

			else # t>1
				for k in 1:nclus_iter[t]
					a_star = sig2h_priors[1] + nh[k,t]/2
					sum_Y = 0.0
					S_kt = findall(j -> Si_iter[j,t] == k, 1:n)
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
						sum_Y += (Y[j,t] - muh_iter[k,t] - eta1_iter[j]*Y[j,t-1] - (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0))^2
					end

					b_star = sig2h_priors[2] + sum_Y/2
					sig2h_iter[k,t] = rand(InverseGamma(a_star, b_star))
				end
			end
			# end # of the @timeit for sigma2h
			############# update beta #############
			# @timeit to " beta " begin # if logging uncomment this line, and the corresponding "end"
			if lk_xPPM
				debug(@showd t)
				if t==1
					sum_Y = zeros(p)
					A_star = I(p)/s2_beta
					for j in 1:n
						X_jt = @view Xlk_covariates[j,:,t]
						sum_Y += (Y[j,t] - muh_iter[Si_iter[j,t],t]) * X_jt / sig2h_iter[Si_iter[j,t],t]
						A_star += (X_jt * X_jt') / sig2h_iter[Si_iter[j,t],t]
					end
					debug(@showd A_star)
					b_star = beta0/s2_beta + sum_Y
					debug(@showd b_star)
					A_star = Symmetric(A_star)
					beta_iter[t] = rand(MvNormalCanon(A_star*b_star, A_star))
					# beta_iter[t] = rand(MvNormal(b_star, inv(Symmetric(A_star))))
					# Symmetric is needed for numerical problems
					# https://discourse.julialang.org/t/isposdef-and-eigvals-do-not-agree/118191/2
					# but A_star is indeed symm and pos def (by construction) so there is no problem
				else
					sum_Y = zeros(p)
					A_star = I(p)/s2_beta
					for j in 1:n
						X_jt = @view Xlk_covariates[j,:,t]
						sum_Y += (Y[j,t] - muh_iter[Si_iter[j,t],t] - eta1_iter[j]*Y[j,t-1]) * X_jt / sig2h_iter[Si_iter[j,t],t]
						A_star += (X_jt * X_jt') / sig2h_iter[Si_iter[j,t],t]
					end
					debug(@showd A_star)
					b_star = beta0/s2_beta + sum_Y
					debug(@showd b_star)
					A_star = Symmetric(A_star)
					beta_iter[t] = rand(MvNormalCanon(A_star*b_star, A_star))
					# beta_iter[t] = rand(MvNormal(b_star, inv(Symmetric(A_star))))
				end
				pretty_log("beta_iter")
				println(beta_iter)
			end
						
			# end # of the @timeit for beta
			############# update theta #############
			# @timeit to " theta " begin # if logging uncomment this line, and the corresponding "end"
			aux1 = 1 / (lambda2_iter*(1-phi1_iter^2))
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
			a_star = tau2_priors[1] + kt
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
					ll_old = 0.0; ll_new = 0.0
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
				sumg = sum(gamma_iter)
				a_star = alpha_priors[1] + sumg
				b_star = alpha_priors[2] + n*T - sumg
				alpha_iter = rand(Beta(a_star, b_star))

			elseif time_specific_alpha==true && unit_specific_alpha==false
				# a vector in time
				for t in 1:T
					sumg = sum(gamma_iter[:,t])
					a_star = alpha_priors[1] + sumg
					b_star = alpha_priors[2] + n - sumg
					alpha_iter[t] = rand(Beta(a_star, b_star))
				end

			elseif time_specific_alpha==false && unit_specific_alpha==true
				# a vector in units
				for j in 1:n
					sumg = sum(gamma_iter[j,:])
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
		b_star = lambda2_priors[2] + (theta_iter[1] - phi0_iter)^2 / 2 + aux1/2
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
					if t==1
						llike[j,t,i_out] = logpdf(Normal(
							muh_jt + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							# muh_jt + (lk_xPPM ? dot(Xlk_covariates[j,:,t], beta_iter[t]) : 0),
							sqrt(sig2h_jt)
							), Y[j,t])
						fitted[j,t,i_out] = muh_jt + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0)
						# fitted[j,t,i_out] = muh_jt + (lk_xPPM ? dot(Xlk_covariates[j,:,t], beta_iter[t]) : 0)
					else # t>1
						llike[j,t,i_out] = logpdf(Normal(
							muh_jt + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0),
							# muh_jt + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(Xlk_covariates[j,:,t], beta_iter[t]) : 0),
							sqrt(sig2h_jt*(1-eta1_iter[j]^2))
							), Y[j,t])
						fitted[j,t,i_out] = muh_jt + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(view(Xlk_covariates,j,:,t), beta_iter[t]) : 0)
						# fitted[j,t,i_out] = muh_jt + eta1_iter[j]*Y[j,t-1] + (lk_xPPM ? dot(Xlk_covariates[j,:,t], beta_iter[t]) : 0)
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

	println("LPML: $LPML (the higher the better)")
	
	# adjust mean variables
	mean_likelhd ./= nout
	mean_loglikelhd./= nout
	for j in 1:n
		for t in 1:T
			WAIC += 2*mean_loglikelhd[j,t] - log(mean_likelhd[j,t])
		end
	end
	WAIC *= -2
	println("WAIC: $WAIC (the lower the better)")

	if logging
		debug(@showd to)
	end
	close(log_file)

	return Si_out, Int.(gamma_out), alpha_out, sigma2h_out, muh_out, include_eta1 ? eta1_out : NaN,
		lk_xPPM ? beta_out : NaN, theta_out, tau2_out, phi0_out, include_phi1 ? phi1_out : NaN, lambda2_out,
		fitted, llike, LPML, WAIC


# catch e
	# println(e)
	# close(log_file)
# end

# global_logger(ConsoleLogger())	
end

function close_log_file()
	close(log_file)
end