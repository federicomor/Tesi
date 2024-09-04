using LinearAlgebra
using BenchmarkTools
using Cthulhu

include("../utils.jl")

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

cohesion1(s1,s2,alpha,lg=true)
cohesion2(s1,s2,a,lg=true)

cohesion3(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi,Cohesion=3,lg=true)
cohesion4(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi,Cohesion=4,lg=true)

cohesion5(s1,s2,phi,lg=true)
cohesion6(s1,s2,phi,lg=true)


