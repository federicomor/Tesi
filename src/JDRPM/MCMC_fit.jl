using Distributions
using Statistics
using LinearAlgebra
using Random
using Logging

############# for debug #############
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
	sigma2h_priors::Vector{Float64},  # Prior parameters for sigma2h ∼ 
	eta_priors::Vector{Float64},      # Prior parameters for eta ∼ 
	beta_priors::Vector{Float64},     # Prior parameters for beta ∼ 
	tau2_priors::Vector{Float64},     # Prior parameters for tau2 ∼ invGamma(a_tau2, b_tau2)
	phi0_priors::Vector{Float64},     # Prior parameters for phi0 ∼ N(m0, s0^2)
	phi1_priors::Vector{Float64},     # Prior parameters for phi1 ∼ 
	lambda2_priors::Vector{Float64},  # Prior parameters for lambda2 ∼ invGamma(a_lambda2, b_lambda2)
	alpha_priors::Matrix{Float64},    # Prior parameters for alpha ∼ 
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
	if spatial_cohesion_idx == 1 && !(sp_params isa Real) 
		@error "wrong params for spatial cohesion 1.\nExpected input form: (Real). Received: $(typeof(sp_params))." _file=""
		return
	elseif spatial_cohesion_idx == 2 && !(sp_params isa Real)
		@error "wrong params for spatial cohesion 2.\nExpected input form: (Real). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion_idx == 3 && !(sp_params isa Vector && length.(sp_params) == [2,1,1,4])
		@error "wrong params for spatial cohesion 3.\nExpected input form: (1x2 Vector, Real, Real, 2x2 Matrix). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion_idx == 4 && !(sp_params isa Vector && length.(sp_params) == [2,1,1,4])
		@error "wrong params for spatial cohesion 4.\nExpected input form: (1x2 Vector, Real, Real, 2x2 Matrix). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion_idx == 5 && !(sp_params isa Real)
		@error "wrong params for spatial cohesion 5.\nExpected input form: (Real). Received: $(typeof(sp_params))."
		return
	elseif spatial_cohesion_idx == 6 && !(sp_params isa Real)
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
	println("fitting $draws iterates \n(for $nout actual iterates)")
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
		beta_out = zeros(nout,p,T)
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
	Si_iter = ones(Int64,n,T_star) # label assignements for units j at time t
	Si_iter[:,end] .= 0
	gamma_iter = zeros(n,T_star)
	alpha_iter = ones(T_star)*starting_alpha
	sig2h_iter = ones(n,T_star)
	muh_iter = zeros(n,T_star)
	eta1_iter = zeros(n)
	if xPPM
		beta_iter = rand(MvNormal(beta_priors[1:end-1],beta_priors[end]*Matrix(I, p, p)))
	end
	theta_iter = rand(Normal(),T_star)
	tau2_iter = rand(InverseGamma(tau2_priors...),T_star)
	phi0_iter = (T==2) ? mean(Normal(phi0_priors...)) : rand(Normal(phi0_priors...))
	phi1_iter = rand(Uniform(-1,1))*include_phi1
	lambda2_iter = (T==2) ? mean(InverseGamma(lambda2_priors...)) : rand(InverseGamma(lambda2_priors...))
	
	nh = zeros(Int,n,T_star) # numerosity of each cluster k at time t
	# dimension n since at most there are n clusters (all singletons)
	nh[1,:] .= n
	nclus_iter = ones(Int,T_star) # how many clusters there are at time t


	############# auxiliary working variables #############


	############# testing (for now) on random values #############
	# Si_iter = rand(collect(1:n),n,T_star)
	# gamma_iter = rand((0,1),n,T_star)
	
	
	############# start MCMC algorithm #############
	debug("LOG FILE\ncurrent seed = $seed") # ▶►▸

	for i in 1:draws
		debug("\n▶ iteration $i")
		# println("\n▶ iteration $i")
		# print("iteration $i\r") # use only this when all finished, for a shorter feedback

		for t in 1:T
			debug("► time $t")

			############# update gamma #############
			for j in 1:n
				debug("[update gamma]")
				debug("▸ subject $j")
				if t==1 
					gamma_iter[j,t] = 0
				else
					# we want to find ρ_t^{R_t(-j)} ...
					indexes = findall(jj -> gamma_iter[jj, t] == 1, setdiff(1:n,j))
					Si_red = Si_iter[indexes, t]
					Si_red1 = copy(Si_red)
					push!(Si_red1, Si_iter[j,t]) # ... and ρ_t^R_t(+j)}

					# get also the reduced spatial info if sPPM model
					if sPPM
						sp1_red = sp1[indexes]
						sp2_red = sp2[indexes]
					end

					# compute n_red's and nclus_red's and relabel
					n_red = length(Si_red) # = "n" relative to here, i.e. the sub-partition size
					n_red1 = length(Si_red1)
					@assert n_red1 == n_red+1
					relabel!(Si_red,n_red)
					relabel!(Si_red1,n_red1)
					nclus_red = isempty(Si_red) ? 0 : maximum(Si_red) # = number of clusters
					nclus_red1 = maximum(Si_red1)

					# save the label of the current working-on unit j
					j_label = Si_red1[end]

					# compute also nh_red's
					nh_red = zeros(Int64,nclus_red); nh_red1 = zeros(Int64,nclus_red1)
					for jj in 1:n_red
						nh_red[Si_red[jj]] += 1 # = numerosities for each cluster label
						nh_red1[Si_red1[jj]] += 1
					end
					nh_red1[Si_red1[end]] += 1 # account for the last added unit j
					debug(@showd Si_red n_red nclus_red nh_red Si_red1 n_red1 nclus_red1 nh_red1 j_label)

					# start computing weights
					lg_weights = zeros(nclus_red+1)
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
					# ... or unit j can create a singleton
					lCn = 0.0
					if sPPM
						lCn = spatial_cohesion(spatial_cohesion_idx, [sp1[j]], [sp2[j]], sp_params, lg=true, M=M_dp)
					end
					lg_weights[nclus_red+1] = log(M_dp) + lCn
					debug(@showd nclus_red lg_weights)


					# now use the weights towards sampling the new gamma_jt
					max_ph = maximum(lg_weights)
					sum_ph = 0.0

					# exponentiate...
					for k in 1:length(lg_weights)
						 # for numerical purposes we subract max_ph
						lg_weights[k] = exp(lg_weights[k] - max_ph)
						sum_ph += lg_weights[k]
					end
					# ... and normalize
					lg_weights ./= sum_ph
					# compute probh
					probh = alpha_iter[t] / (alpha_iter[t] + (1 - alpha_iter[t]) * lg_weights[j_label])
					debug(@showd probh alpha_iter[t])


					# compatibility check for gamma transition
					if gamma_iter[j, t] == 0
						# we want to find ρ_(t-1)^{R_t(+j)} ...
						indexes = findall(jj -> gamma_iter[jj, t]==1, 1:n)
						Si_comp1 = Si_iter[indexes, t-1]
						Si_comp2 = Si_iter[indexes, t] # ... and ρ_t^R_t(+j)}

						rho_comp = compatibility(Si_comp1, Si_comp2)
						if rho_comp == 0
							probh = 0
						end
					end
					# sample the new gamma
					gt = rand(Bernoulli(probh))
					gamma_iter[j, t] = gt
				end
			end

			############# update rho #############
			debug("[update rho]")
			# we only update the partition for the units which can move (i.e. with gamma_jt=0)
			movable_units = findall(j -> gamma_iter[j,t]==0, 1:n)
	
			for j in movable_units
				# remove unit j from the cluster she is currently in

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
				nh_tmp = zeros(Int,nclus_iter[t]+1)
				for jj in setdiff(1:n,j)
					nh_tmp[rho_tmp[jj]] += 1
				end
				nclus_temp = 0

				println(@showd nclus_iter t)
				# we now simulate the unit j to be assigned to one of the existing clusters...
				for k in 1:nclus_iter[t]
					rho_tmp[j] = k
					indexes = findall(j -> gamma_iter[j,t+1]==1, 1:n)
					# we check the compatibility between ρ_t^{h=k,R_(t+1)} ...
					Si_comp1 = rho_tmp[indexes]
					Si_comp2 = Si_iter[indexes,t+1] # and ρ_(t+1)^{R_(t+1)}
					rho_comp = compatibility(Si_comp1, Si_comp2)
					
					if rho_comp != 1
						ph[k] = log(0) # assignment to cluster k is not compatible
					else
						# update params for "rho_jt = k" simulation
						nh_tmp[k] += 1
						nclus_temp = sum(nh_tmp .> 0) # = number of clusters

						lpp = 0.0
						for kk in 1:nclus_temp
							if sPPM
								indexes = findall(jj -> rho_tmp[jj]==kk, 1:n)
								s1n = sp1[indexes]
								s2n = sp2[indexes]
								lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
							end
							lpp += log(M_dp) + lgamma(nh_tmp[kk])
							# lpp += log(M_dp) + lgamma(length(indexes)) # same
						end

						# debug case
						ph[k] = 0.5
						# real case
						# if t==1
						# 	ph[k] = loglikelihood(Normal(
						# 		muh_iter[k,t],
						# 		sqrt(sig2h_iter[k,t])),
						# 		Y[j,t]) + lpp
						# else
						# 	ph[k] = loglikelihood(Normal(
						# 		muh_iter[k,t] + eta1_iter[j]*Y[j,t-1],
						# 		sqrt(sig2h_iter[k,t]*(1-eta1_iter[j]^2))),
						# 		Y[j,t]) + lpp
						# end

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
						if sPPM
							indexes = findall(jj -> rho_tmp[jj]==kk, 1:n)
							s1n = sp1[indexes]
							s2n = sp2[indexes]
							lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
						end
						lpp += log(M_dp) + lgamma(nh_tmp[kk])
						# lpp += log(M_dp) + lgamma(length(indexes)) # same
					end
					# debug case
					ph[k] = 0.5
					# real case
					# if t==1
					# 	ph[k] = loglikelihood(Normal(
					# 		muh_draw + eta1_iter[j]*Y[j,t-1],
					# 		sqrt(sig2h_draw)),
					# 		Y[j,t]) + lpp
					# else
					# 	ph[k] = loglikelihood(Normal(
					# 		muh_draw,
					# 		sqrt(sig2h_draw*(1-eta1_iter[j]^2))),
					# 		Y[j,t]) + lpp
					# end

					# # restore params after "rho_jt = k" simulation
					# nh_tmp[k] -= 1
				end


				println(ph, " (before exp and norm)")
				# now exponentiate the weights...
				max_ph = maximum(ph)
				sum_ph = 0.0
				for k in 1:length(ph)
					# for numerical purposes we subract max_ph
					ph[k] = exp(ph[k] - max_ph)
					sum_ph += ph[k]
				end
				# ... and normalize them
				ph ./= sum_ph
				println(ph, " (after exp and norm)")

				# now sample the new label Si_iter[j,t]
				uu = rand(Uniform(0,1))
				cph = cumsum(ph)
				cph[end] = 1.0 # fix numerical problems of having sums like 0.999999etc
				new_label = 0
				for k in 1:length(ph)
					if uu <= cph[k]
						new_label = k
						break
					end
				end

				if new_label <= nclus_iter[t]
					# we enter an existing cluster
					println("we enter an existing cluster")
					Si_iter[j, t] = new_label
					nh[new_label, t] += 1
					# rho_tmp[j] = new_label # useless
				else
					# we create a new singleton cluster
					println("we create a new singleton cluster")
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
				Si_tmp = copy(Si_iter[:,t])
				Si_relab, nh_reorder, old_lab = relabel_full(Si_tmp,n)
				# eg:             Original labels (Si): 4 2 1 1 1 3 1 4 5 
				#          Relabeled groups (Si_relab): 1 2 3 3 3 4 3 1 5
				# Reordered cluster sizes (nh_reorder): 2 1 4 1 1 0 0 0 0
				# 	              Old labels (old_lab): 4 2 1 3 5 0 0 0 0 
				# - Si_relab gives the relabelled partition
				# - nh_reorder gives the numerosities of the relabelled partition, ie "nh_reorder[k] = #(units of new cluster k)"
				# - old_lab tells "the index in position i (which before was cluster i) is now called cluster old_lab[i]"

				# now fix everything (morally permute params)
				Si_iter[:,t] = Si_relab
				# discard the zeros at the end of the auxiliary vectors nh_reorde and old_lab
				len = findlast(x -> x != 0, nh_reorder)
				muh_iter[1:len,t] = muh_iter[old_lab[1:len],t]
				sig2h_iter[1:len,t] = sig2h_iter[old_lab[1:len],t]
				nh[1:len,t] = nh_reorder[1:len]

			end # for j in movable_units


			############# update muh #############
			if t==1
				for k in 1:nclus_iter[t]
					sum_Y = 0.0
					for j in 1:n
						if Si_iter[j,t]==k 
							sum_Y += Y[j,t]
						end
					end
					sig2_post = 1 / (1/tau2_iter[t] + sum(Si_iter[:,t] .== k)/sig2h_iter[k,t])
					mu_post = sig2_post * (theta_iter[t]/tau2_iter[t] + sum_Y/sig2h_iter[k,t])

					muh_iter[k] = rand(Normal(mu_post,sqrt(sig2_post)))
				end

			else # t>1
				for k in 1:nclus_iter[t]
					sum_Y = 0.0
					sum_e2 = 0.0
					for j in 1:n
						if Si_iter[j,t]==k 
							aux1 = 1 / (1-eta1_iter[j]^2)
							sum_e2 += aux1
							sum_Y += (Y[j,t] - eta1_iter[j]*Y[j,t-1]) * aux1
						end
					end

					sig2_post = 1 / (1/tau2_iter[t] + sum_e2/sig2h_iter[k,t]) 
					mu_post = sig2_post * (theta_iter[t]/tau2_iter[t] + sum_Y/sig2h_iter[k,t])

					muh_iter[k] = rand(Normal(mu_post,sqrt(sig2_post)))
				end
			end
			
			############# update sigma2h #############

			
			############# update theta #############
			aux1 = 1 / (lambda2_iter*(1-phi1_iter^2))
			kt = nclus_iter[t]
			sum_mu = sum(muh_iter[1:kt,t])

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
			a_star = tau2_priors[1] + kt
			b_star = tau2_priors[2] + sum((muh_iter[k,t] - theta_iter[t])^2 for k=1:kt) / 2
			tau2_iter[t] = rand(InverseGamma(a_star, b_star))


		end # for t in 1:T
		
		############# update eta1 #############

		
		############# update alpha #############

		
		############# update phi0 #############
		aux1 = 1/lambda2_iter
		sig2_post = 1 / ( 1/phi0_priors[2] + aux1 * (1 + (T-1)*(1-phi1_iter)/(1+phi1_iter)) )
		mu_post = sig2_post * ( phi0_priors[1]/phi0_priors[2] + theta_iter[1]*aux1 + aux1/(1+phi1_iter)*sum(theta_iter[t] - phi1_iter*theta_iter[t-1] for t=2:T) )

		
		############# update phi1 #############

		
		############# update lambda2 #############
		a_star = lambda2_priors[1] + T/2
		b_star = lambda2_priors[2] + (theta_iter[1] - phi0_iter)^2 / 2 + sum((theta_iter[t] - (1-phi1_iter)*phi0_iter - phi1_iter*theta_iter[t-1])^2 for t=2:T) / 2
	

		############# save MCMC iterates #############


	end # for i in 1:draws
	println("\ndone!")
	
close(log_file)
# global_logger(ConsoleLogger())	
end

