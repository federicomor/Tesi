# Copy-Item "C:/Users/feder/AppData/Local/Temp/jl_cdGNc7W9m9.html" -Destination ".\Tesi\src\JDRPM\profiling\"

using ProfileCanvas
# @profview ECC
# @profview_allocs ECC

using Distributions
using Statistics
using LinearAlgebra
using Random
using SpecialFunctions
using BenchmarkTools

function G2a(a::Real, lg::Bool)
	out = log(π) + lgamma(a) + lgamma(a - 0.5)
	return lg ? out : exp(out)
end

function cohesion3_4(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}, S::Matrix{Float64}, vtemp::Vector{Float64}; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	sp = [s1 s2]
	sbar = vec(mean(sp, dims=1))
	S = zeros(2,2)
	for i in 1:sdim
		vtemp = sp[i,:] - sbar
		S += (vtemp)*(vtemp)'
	end
	kn = k0+sdim
	vn = v0+sdim
	Psi_n = zeros(2,2)
	Psi_n = Psi + S + (k0*sdim)/(k0+sdim)*(sbar-mu_0)*(sbar-mu_0)'
	
	out = -sdim * log(π) + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * logdet(Psi) - 0.5 * vn * logdet(Psi_n) + log(k0) - log(kn)
	return lg ? out : exp(out)
end

function cohesion3_4_v2(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}, sbar::Vector{Float64}, S::Matrix{Float64}, vtemp::Vector{Float64}, Psi_n::Matrix{Float64}; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	sp = [s1 s2]
	sbar = [mean(s1), mean(s2)]
	fill!(S,0)
	for i in 1:sdim
		vtemp = sp[i,:] - sbar
		S += (vtemp)*(vtemp)'
	end
	kn = k0+sdim
	vn = v0+sdim
	vtemp2 = sbar-mu_0
	Psi_n = Psi + S + (k0*sdim)/(k0+sdim)*vtemp2*vtemp2'

	out = -sdim * log(π) + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * logdet(Psi) - 0.5 * vn * logdet(Psi_n) + log(k0) - log(kn)
	return lg ? out : exp(out)
end

spatial_cohesion_idx = 3
n = 20
Random.seed!(1)
s1n = rand(n)
s2n = rand(n)
mu_0 = [0.0, 0.0]
k0 = 1.0
v0 = 5.0
Psi =[1.0 0.0; 0.0 1.0]
M_dp = 1

function to_profile(iterations::Int)
	result = 0.0
	# for _ in 1:iterations 
	# 	s1n = rand(n)
	# 	s2n = rand(n)
	# 	result += cohesion3_4(s1n,s2n,mu_0,k0,v0,Psi, lg=true, M=M_dp)
	# end
	sbar = zeros(2)
	S = zeros(2,2)
	vtemp = zeros(2)
	Psi_n = zeros(2,2)
	for _ in 1:iterations 
		s1n = rand(n)
		s2n = rand(n)
		result += cohesion3_4_v2(s1n,s2n,mu_0,k0,v0,Psi,sbar,S,vtemp,Psi_n,lg=true, M=M_dp)
	end
	return result
end

@benchmark to_profile(1000)
# ProfileCanvas.@profview_allocs to_profile(1000)

# ProfileCanvas.@profview spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)
# ProfileCanvas.@profview_allocs spatial_cohesion(spatial_cohesion_idx, s1n, s2n, sp_params, lg=true, M=M_dp)

# Copy-Item "C:/Users/feder/AppData/Local/Temp/jl_cdGNc7W9m9.html" -Destination ".\Tesi\src\JDRPM\profiling\"



