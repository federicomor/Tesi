module JDRPM

export my_example, MCMC_fit, close_log_file, test_R_to_J_conversion

function test_R_to_J_conversion(variable)
	println("typeof variable: $(typeof(variable))")
	println("size: $(size(variable))")
	println("content:\n",variable)
end

function my_example(iters::Number)
	result = 0.0
	p = Progress(round(Int64(iters));
		showspeed=true,
		dt=0.5,
        # barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        # barglyphs=BarGlyphs("[=> ]"),
        color=:yellow,
        barlen=0,
        # enabled=false,
        )
	for i in 1:iters
		result += rand() # for example
		sleep(0.05)
		next!(p;
			# showvalues = [("iteration",i), (:result,result)]
			)
	end
	return result
end

include("../MCMC_fit.jl")

# a function to call a quick run to compile the code
# otherwise if we do a "serious" run as first run it will be slow
# since all first call to something in julia are slower as they need to be compiled
# then all subsequent runs will be faster 
function trigger_compilation()
	Random.seed!(1)
	N = 4; T = 4
	y = rand(N,T)
	sp = rand(N,2)
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
	spatial_cohesion_idx = 3
	mu0 = 0.
	k0 = 1.
	v0 = 5.
	L0 = 1.
	niter = 10
	burnin = 5
	thin = 1
	seed = 111.0
	out = MCMC_fit(
		Y=y,              
		sp_coords = sp,
		M_dp = 1.0,
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
		tau2_priors = [a_tau,b_tau],
		phi0_priors = [m0_phi0,s20_phi0],
		phi1_priors = sig_mh_phi1,
		lambda2_priors = [a_lambda,b_lambda],
		alpha_priors = [a_alpha,b_alpha],
		spatial_cohesion_idx = spatial_cohesion_idx,
		sp_params = [[mu0,mu0],k0,v0,[L0 0.0; 0.0 L0]],
		draws = niter,                    
		burnin = burnin,                   
		thin = thin,                     
		logging = false,
		seed = seed
		);
	return nothing
end

end # module JDRPM
