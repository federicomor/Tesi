using LinearAlgebra
using BenchmarkTools
using Cthulhu

include("../utils.jl")

n = 20
s1 = rand(n)
s2 = rand(n)
alpha = 0.13
a = 10
mu_0 = [1.,2.]
k0 = 0.5
v0 = 0.5
Psi = [1. 2.; 2 4]
phi=0.5

sp_params = [mu_0, k0, v0, Psi]

cohesion1(s1,s2,alpha,lg=true)
cohesion2(s1,s2,a,lg=true)
cohesion3(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion4(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion5(s1,s2,phi,lg=true)
cohesion6(s1,s2,phi,lg=true)

sp_params = [alpha]
@code_warntype spatial_cohesion(1,s1,s2,sp_params,lg=false,M=1)

println()
y1 = cohesion1(s1,s2,alpha,true)
y2 = cohesion2(s1,s2,a,true)
y3 = cohesion3_4(s1, s2, mu_0, k0, v0, Psi,3,true)
y4 = cohesion3_4(s1, s2, mu_0, k0, v0, Psi,4,true)
y5 = cohesion5(s1,s2,phi,true)
y6 = cohesion6(s1,s2,phi,true)



