using ProfileCanvas
# using JET
using Cthulhu

begin
using Random
Random.seed!(1)
N = 20; T = 50
# N = 50; T = 18
y = rand(N,T)
sp = rand(N,2)
p = 2
X_cl = rand(N,p,T)

# params
m0_phi0 = 0.
s20_phi0 = 10.
a_sigma  = 2.; b_sigma  = 2.
a_tau    = 2.; b_tau    = 2.
a_lambda = 2.; b_lambda = 2.
eta1_scale = 0.9
sig_mh_eta1 = 0.1
sig_mh_phi1 = 0.1
update_eta1 = true
update_phi1 = true
a_alpha = 2.; b_alpha = 2.
time_specific_alpha = true
# now space
spatial_cohesion_idx = 3
mu0 = 0.
k0 = 1.
v0 = 5.
L0 = 1.

# niter = 5000. # profile canvas
niter = 1000 # real for these tests
# niter = 20. # for tracking allocs
burnin = 500
thin = 5
seed = 123.0
end


# include("../old/MCMC_fit.jl")
# include("../new/MCMC_fit.jl")
# include("../new/utils.jl")
# include("../MCMC_fit.jl")

#=

=#

# include("../old/MCMC_fit.jl")
include("../past versions/prelease/MCMC_fit.jl")
include("../past versions/old/MCMC_fit.jl")
include("../past versions/cova testing time/MCMC_fit.jl")

# include("../new/MCMC_fit.jl")

include("../MCMC_fit.jl")

#=
# LPML: -573.5247999716885 (the higher the better)
# WAIC: 647.080357797959 (the lower the better)
# acceptance ratio eta1: 95.18%
# acceptance ratio phi1: 97.90%
# space+cova
Elapsed time: 2 seconds, 680 milliseconds
LPML: -694.8307767293472 (the higher the better)
WAIC: 672.4901618487447 (the lower the better)
acceptance ratio eta1: 94.94%
acceptance ratio phi1: 97.20%
  2.686072 seconds (22.53 M allocations: 1.715 GiB, 6.53% gc time)
elapsed time (ns):  2686072100
gc time (ns):       175324000
bytes allocated:    1841554000
pool allocs:        22524595
non-pool GC allocs: 1014
malloc() calls:     11
free() calls:       11
minor collections:  39
full collections:   1
=#


# ProfileCanvas.@profview MCMC_fit(
# ProfileCanvas.@profview_allocs MCMC_fit(
# @profview_allocs MCMC_fit(
# @code_warntype MCMC_fit(
# @report_opt MCMC_fit(
# @profview_allocs out = MCMC_fit(
# @descend MCMC_fit(

niter = 10000. # real for these tests
burnin = 5000.
thin = 50.

niter = 6. # real for these tests
burnin = 2.
thin = 1.

@descend MCMC_fit(
# @timev MCMC_fit(
# ProfileCanvas.@profview_allocs MCMC_fit(
# ProfileCanvas.@profview MCMC_fit(
# out = MCMC_fit(
	Y=y,              
	sp_coords = sp,

	# Xcl_covariates = X_cl,
	# covariate_similarity_idx = 4,
	# cv_params = [0.,2.,1.,2.],  
	# cv_params = missing,  

	M_dp = 1.0,                     
	initial_partition = missing,
	# Xlk_covariates = missing,
	starting_alpha = 0.5,         
	unit_specific_alpha = false,       
	time_specific_alpha = time_specific_alpha,     

	update_alpha = true,             
	include_eta1 = true,                    
	include_phi1 = true,
	update_eta1 = update_eta1,                    
	update_phi1 = update_phi1,

	sig2h_priors = [a_sigma,b_sigma],
	eta1_priors = [eta1_scale,sig_mh_eta1],
	# beta_priors = c(rep(1,p),2),
	# beta_priors = missing,
	tau2_priors = [a_tau,b_tau],
	phi0_priors = [m0_phi0,s20_phi0],
	phi1_priors = sig_mh_phi1,
	lambda2_priors = [a_lambda,b_lambda],
	alpha_priors = [a_alpha,b_alpha],
	
	spatial_cohesion_idx = spatial_cohesion_idx,
	sp_params = [[mu0,mu0],k0,v0,[L0 0.0; 0.0 L0]],
	# spatial_cohesion_idx = 1,
	# sp_params = [0.5],

	draws = niter * 1.0,                    
	burnin = burnin * 1.0,                   
	thin = thin * 1.0,                     
	logging = false,
	seed = seed
)


#=
Elapsed time: 3 seconds, 124 milliseconds
LPML: -422.2917598847498 (the higher the better)
WAIC: 829.7457081318109 (the lower the better)
acceptance ratio eta1: 72.00%
acceptance ratio phi1: 64.40%
official version #############################
Elapsed time: 4 seconds, 515 milliseconds
(25.35 M allocations: 2.236 GiB, 4.64%% gc time)
elapsed time (ns):  4431394100
gc time (ns):       214217200
bytes allocated:    2400526000
pool allocs:        25350195
non-pool GC allocs: 112
malloc() calls:     10
free() calls:       0
minor collections:  37
full collections:   0

pre(re)lease version #############################
Elapsed time: 54 seconds, 190 milliseconds
(331.76 M allocations: 21.547 GiB, 5.80% gc time)
elapsed time (ns):  55209293200
gc time (ns):       3203797000
bytes allocated:    23135946288
pool allocs:        330762753
non-pool GC allocs: 1001718
malloc() calls:     10
free() calls:       0
minor collections:  267
full collections:   0

old version #############################
Elapsed time: 17 seconds, 314 milliseconds
(193.25 M allocations: 18.233 GiB, 9.22% gc time)
elapsed time (ns):  17322894900
gc time (ns):       1597300100
bytes allocated:    19577849360
pool allocs:        192243554
non-pool GC allocs: 1001720
malloc() calls:     10
free() calls:       0
minor collections:  427
full collections:   0

cova time version #############################
Elapsed time: 17 seconds, 16 milliseconds
(169.99 M allocations: 17.248 GiB, 9.01% gc time)
elapsed time (ns):  18037863700
gc time (ns):       1625568400
bytes allocated:    18519403520
pool allocs:        168985729
non-pool GC allocs: 1001718
malloc() calls:     10
free() calls:       0
minor collections:  381
full collections:   0
=#



# begin
# lg_weights = [-1000.0, -2000, -1500, -100, -200]
# ph_tmp = copy(lg_weights)
# sort!(ph_tmp)
# max_ph = ph_tmp[end]
# sum_ph = 0.0
# # exponentiate...
# for k in 1:length(lg_weights)
# 	 # for numeric purposes we subract max_ph
# 	# lg_weights[k] = exp(lg_weights[k]) # without
# 	lg_weights[k] = exp(lg_weights[k]-max_ph) # with
# 	global sum_ph += lg_weights[k]
# end
# @show lg_weights
# # ... and normalize
# for k in 1:length(lg_weights)
# 	lg_weights[k] /= sum_ph
# end
# @show lg_weights
# end



# begin 
# 	ph = rand((0:100),5) * 1.0
# 	ph ./= sum(ph)
# 	uu = rand()
# 	# uu = 1
# 	begin cph = cumsum(ph)
# 		@show ph
# 		@show cph
# 			@assert isapprox(cph[end], 1, atol=1e-6)
# 		println(findfirst(x -> x==1, uu .<= cph))
# 		end
# 		@show uu
# 		begin cph = cumsum(ph)
# 		iaux = 0
# 		for k in 1:length(ph)
# 			if uu < cph[k]
# 				global iaux = k
# 				break
# 			end
# 		end
# 		println(iaux)
# 	end
# end


# begin
# A_star = [0.5695141296842557 -0.5958323106341032; -0.5958323106341032 5.607107634272117]
# invA_star = inv(A_star)
# # invA_star = [1.9755086272843259 0.20992496432168561; 0.20992496432168564 0.20065248429953192]
# using LinearAlgebra
# println(eigvals(A_star)) # [0.5000000000009021, 5.67662176395547] all >0
# println(eigvals(invA_star)) # [0.1761611115874664, 1.9999999999963913] still all >0
# println(isposdef(A_star)) # true
# isposdef(invA_star) # false
# end


# begin
# x = 1
# println("(before the loop) x=$x")
# for i in 1:3
# 	global x = 5
# 	y = x+i
# 	println("($i) y=$y x=$x")
# end
# println("(after the loop) x=$x")
# end