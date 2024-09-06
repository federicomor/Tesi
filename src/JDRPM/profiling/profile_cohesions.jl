using LinearAlgebra
using BenchmarkTools
include("../utils.jl")

function spatial_cohesion_splat(idx::Real, s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, sp_params::Vector; lg::Bool, M::Real)
	idx==1.0 && return cohesion1(s1,s2,sp_params...,lg=lg,M=M) 
	idx==2.0 && return cohesion2(s1,s2,sp_params...,lg=lg,M=M) 
	idx==3.0 && return cohesion3(s1,s2,sp_params...,lg=lg,M=M) 
	idx==4.0 && return cohesion3(s1,s2,sp_params...,lg=lg,M=M) 
	idx==5.0 && return cohesion5(s1,s2,sp_params...,lg=lg,M=M) 
	idx==6.0 && return cohesion6(s1,s2,sp_params...,lg=lg,M=M) 
end

n = 20
s1 = rand(n)
s2 = rand(n)
alpha = 0.13
a = 10
mu_0 = [0.,0.]
k0 = 0.1
v0 = 0.1
Psi = [1. 0.2; 0.2 1]
phi=0.5

sp_params = [mu_0, k0, v0, Psi]

@btime spatial_cohesion(4,s1,s2,sp_params,lg=true,M=1)
@btime spatial_cohesion_splat(4,s1,s2,sp_params,lg=true,M=1)

@code_llvm spatial_cohesion(4,s1,s2,sp_params,lg=true,M=1)
@code_llvm spatial_cohesion_splat(4,s1,s2,sp_params,lg=true,M=1)




cohesion1(s1,s2,alpha,lg=true)
cohesion2(s1,s2,a,lg=true)

cohesion3(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi,Cohesion=3,lg=true)
cohesion4(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi,Cohesion=4,lg=true)

@timev cohesion5(s1,s2,phi,lg=true)
@timev cohesion6(s1,s2,phi,lg=true)


