include("MCMC_fit.jl")
include("utils.jl")

begin
n = 50; T=12
n_labels = 10
rho_tmp = rand((0:n_labels),n)
kk = 1
sp1 = rand(n)
sp2 = rand(n)
spatial_cohesion_idx = 3
sp_params = Any[[1.0, 1.0], 5.0, 5.0, [1.0 0.0; 0.0 1.0]]
M_dp = 1.0

function compute_cohesion_v1(n,T,rho_tmp,n_labels,spatial_cohesion_idx,sp_params)
	lpp = 0.0
	for kk in 1:n_labels
		indexes = findall(jj -> rho_tmp[jj]==kk, 1:n)
		s1n = @view sp1[indexes]
		s2n = @view sp2[indexes]
		lpp += spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
	end
	return lpp
end

function compute_cohesion_v2(n,T,rho_tmp,n_labels,spatial_cohesion_idx,sp_params)
	lpp = 0.0
	
	sp = zeros(n,2)
	S = zeros(2,2)

	for kk in 1:n_labels
		indexes = findall(jj -> rho_tmp[jj]==kk, 1:n)
		s1n = @view sp1[indexes]
		s2n = @view sp2[indexes]
		lpp += cohesion3_4_v2(s1n, s2n, sp, S, sp_params..., lg=true, M=M_dp,Cohesion=3) 
	end
	return lpp
end

using BenchmarkTools
@btime compute_cohesion_v1(n,T,rho_tmp,n_labels,spatial_cohesion_idx,sp_params)
@btime compute_cohesion_v2(n,T,rho_tmp,n_labels,spatial_cohesion_idx,sp_params)

@timev compute_cohesion_v1(n,T,rho_tmp,n_labels,spatial_cohesion_idx,sp_params);
@timev compute_cohesion_v2(n,T,rho_tmp,n_labels,spatial_cohesion_idx,sp_params);
end