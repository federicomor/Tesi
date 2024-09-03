# Copy-Item "C:/Users/feder/AppData/Local/Temp/jl_cdGNc7W9m9.html" -Destination ".\Tesi\src\JDRPM\profiling\"

# using ProfileCanvas
# @profview ECC
# @profview_allocs ECC

using Distributions
using Statistics
using LinearAlgebra
using Random
using SpecialFunctions
using BenchmarkTools

const logpi = log(π) 
function G2a(a::Real, lg::Bool)
	out = logpi + lgamma(a) + lgamma(a - 0.5)
	return lg ? out : exp(out)
end

# paper 3 section 3.1
function cohesion3_4(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}; Cohesion::Int, lg::Bool, M::Real=1.0)
	sdim = length(s1)
	sp = [s1 s2]
	sbar = @SVector [mean(s1), mean(s2)]
	# S = sum( (sp[i,:] - sbar)*(sp[i,:] - sbar)' for i in 1:sdim) # sum is slow
	# FIXED: sum is slow because i didnt initialize S
	# S = zeros(2,2); S = sum( (sp[i,:] - sbar)*(sp[i,:] - sbar)' for i in 1:sdim) # this is fast now
	S = zeros(2,2)
	for i in 1:sdim
		vtemp1 = sp[i,:] - sbar
		S += (vtemp1)*(vtemp1)'
	end
	
	# compute updated parameters
	# for cohesion 3
	kn = k0+sdim
	vn = v0+sdim
	vtemp2 = sbar-mu_0
	Psi_n = Psi + S + (k0*sdim)/(k0+sdim)*vtemp2*vtemp2'
	
	if Cohesion == 3
		out = -sdim * logpi + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * logdet(Psi) - 0.5 * vn * logdet(Psi_n) + log(k0) - log(kn)
		return lg ? out : exp(out)
	end

	# for cohesion 4
	knn = kn+sdim
	vnn = vn+sdim
	mu_n = (k0*mu_0 + sdim*sbar)/(k0+sdim)
	Psi_nn = Psi_n + S + (kn*sdim)/(kn+sdim)*(sbar-mu_n)*(sbar-mu_n)'
	
	if Cohesion == 4
		out = -sdim * logpi + G2a(0.5 * vnn, true) - G2a(0.5 * vn, true) + 0.5 * vn * logdet(Psi_n) - 0.5 * vnn * logdet(Psi_nn) + log(kn) - log(knn)
		return lg ? out : exp(out)
	end
end


##################
# below there is the one "translated" from C, but the above is more readable/efficient/easy
# but are the same, test them on these if you want
# (s1, s2, Psi, Psi_vec) = ([0.06541394973925674, 0.1875903839556078, 0.7065742551867602, 0.8223492385591462], [0.6711654699192571, 0.6199278925430733, 0.36880242735326396, 0.9723482028752322], [2.0 1.0; 1.0 3.0], [2.0, 1.0, 1.0, 3.0])
##################
function cohesion3(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	# Compute sample means
	sbar1 = mean(s1)
	sbar2 = mean(s2)
	# Compute deviations from the sample mean
	S1, S2, S3, S4 = 0.0, 0.0, 0.0, 0.0
	for i in 1:sdim
		s_sbar1 = s1[i] - sbar1
		s_sbar2 = s2[i] - sbar2

		S1 += s_sbar1 * s_sbar1
		S4 += s_sbar2 * s_sbar2
		S2 += s_sbar1 * s_sbar2
	end
	S3 = copy(S2) # to avoid repeating computations

	# Updated parameters for cohesion 3
	kn = k0 + sdim
	vn = v0 + sdim

	auxvec1_1 = sbar1 - mu_0[1]
	auxvec1_2 = sbar2 - mu_0[2]

	auxmat1_1 = auxvec1_1^2
	auxmat1_2 = auxvec1_1 * auxvec1_2
	auxmat1_3 = copy(auxmat1_2)
	auxmat1_4 = auxvec1_2^2

	auxconst1 = k0 * sdim
	auxconst2 = k0 + sdim
	Psi_n_1 = Psi[1] + S1 + auxconst1 / (auxconst2) * auxmat1_1
	Psi_n_2 = Psi[2] + S2 + auxconst1 / (auxconst2) * auxmat1_2
	Psi_n_3 = Psi[3] + S3 + auxconst1 / (auxconst2) * auxmat1_3
	Psi_n_4 = Psi[4] + S4 + auxconst1 / (auxconst2) * auxmat1_4

	detPsi_n = Psi_n_1 * Psi_n_4 - Psi_n_2 * Psi_n_3
	detPsi = Psi[1] * Psi[4] - Psi[2] * Psi[3]

	out = -sdim * logpi + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * log(detPsi) - 0.5 * vn * log(detPsi_n) + log(k0) - log(kn)
	return lg ? out : exp(out)

end

function cohesion4(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	# Compute sample means
	sbar1 = mean(s1)
	sbar2 = mean(s2)
	# Compute deviations from the sample mean
	S1, S2, S3, S4 = 0.0, 0.0, 0.0, 0.0
	for i in 1:sdim
		s_sbar1 = s1[i] - sbar1
		s_sbar2 = s2[i] - sbar2

		S1 += s_sbar1 * s_sbar1
		S4 += s_sbar2 * s_sbar2
		S2 += s_sbar1 * s_sbar2
	end
	S3 = copy(S2) # to avoid repeating computations

	# Updated parameters for cohesion 3
	kn = k0 + sdim
	vn = v0 + sdim

	auxvec1_1 = sbar1 - mu_0[1]
	auxvec1_2 = sbar2 - mu_0[2]

	auxmat1_1 = auxvec1_1^2
	auxmat1_2 = auxvec1_1 * auxvec1_2
	auxmat1_3 = copy(auxmat1_2)
	auxmat1_4 = auxvec1_2^2

	auxconst1 = k0 * sdim
	auxconst2 = k0 + sdim
	Psi_n_1 = Psi[1] + S1 + auxconst1 / (auxconst2) * auxmat1_1
	Psi_n_2 = Psi[2] + S2 + auxconst1 / (auxconst2) * auxmat1_2
	Psi_n_3 = Psi[3] + S3 + auxconst1 / (auxconst2) * auxmat1_3
	Psi_n_4 = Psi[4] + S4 + auxconst1 / (auxconst2) * auxmat1_4

	detPsi_n = Psi_n_1 * Psi_n_4 - Psi_n_2 * Psi_n_3
	
	# Updated parameters for cohesion 4
	knn = kn + sdim
	vnn = vn + sdim

	mu_n_1 = (k0 * mu_0[1] + sdim * sbar1) / (auxconst2)
	mu_n_2 = (k0 * mu_0[2] + sdim * sbar2) / (auxconst2)
	sbar_mun1 = sbar1 - mu_n_1
	sbar_mun2 = sbar2 - mu_n_2

	auxmat2_1 = sbar_mun1^2
	auxmat2_2 = sbar_mun1 * sbar_mun2
	auxmat2_3 = copy(auxmat2_2)
	auxmat2_4 = sbar_mun2^2

	auxconst3 = kn * sdim
	auxconst4 = kn + sdim
	Psi_nn_1 = Psi_n_1 + S1 + auxconst3 / (auxconst4) * auxmat2_1
	Psi_nn_2 = Psi_n_2 + S2 + auxconst3 / (auxconst4) * auxmat2_2
	Psi_nn_3 = Psi_n_3 + S3 + auxconst3 / (auxconst4) * auxmat2_3
	Psi_nn_4 = Psi_n_4 + S4 + auxconst3 / (auxconst4) * auxmat2_4
	detPsi_nn = Psi_nn_1 * Psi_nn_4 - Psi_nn_2 * Psi_nn_3
	
	out = -sdim * logpi + G2a(0.5 * vnn, true) - G2a(0.5 * vn, true) + 0.5 * vn * log(detPsi_n) - 0.5 * vnn * log(detPsi_nn) + log(kn) - log(knn)
	return lg ? out : exp(out)
end

n = 20
Random.seed!(1)
s1n = rand(n)
s2n = rand(n)
mu_0 = [0.0, 0.0]
k0 = 1.0
v0 = 5.0
Psi =[1.0 0.3; 0.3 1.0]
M_dp = 1

println(cohesion3_4(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=3,lg=true) == cohesion3(s1n,s2n,mu_0,k0,v0,Psi,lg=true))
println(cohesion3_4(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=3,lg=false) == cohesion3(s1n,s2n,mu_0,k0,v0,Psi,lg=false))
println(cohesion3_4(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=4,lg=true) == cohesion4(s1n,s2n,mu_0,k0,v0,Psi,lg=true))
println(cohesion3_4(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=4,lg=false) == cohesion4(s1n,s2n,mu_0,k0,v0,Psi,lg=false))

@btime cohesion3_4(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=3,lg=true)
@btime cohesion3_4_C(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=3,lg=true)
@btime cohesion3_4(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=4,lg=true)
@btime cohesion3_4_C(s1n,s2n,mu_0,k0,v0,Psi,Cohesion=4,lg=true)

##################################
using ProfileCanvas

function similarity4(X_jt::AbstractVector{<:Real}, mu_c::Real, lambda_c::Real, a_c::Real, b_c::Real; lg::Bool)
	n = length(X_jt)
	nm = n/2
	xbar = mean(X_jt)
	aux1 = b_c + 0.5 * (sum(X_jt .^ 2) - (n*xbar + lambda_c*mu_c)^2/(n+lambda_c) + lambda_c*mu_c^2 )
	out = -nm*log2pi + 0.5*log(lambda_c/(lambda_c+n)) + lgamma(a_c+nm) - lgamma(a_c) + a_c*log(b_c) + (-a_c-nm)*log(aux1)
	return lg ? out : exp(out)
end

function similarity4_v2(X_jt::AbstractVector{<:Real}, mu_c::Real, lambda_c::Real, a_c::Real, b_c::Real; lg::Bool)
	n = length(X_jt)
	nm = n/2
	xbar = mean(X_jt)
	aux2 = 0.0
	# @inbounds @fastmath @simd for x in X_jt
		# aux2 += x^2
	# end	
	@inbounds @fastmath @simd for i in 1:n
		aux2 += X_jt[i]^2
	end
	# @show aux2
	aux1 = b_c + 0.5 * (aux2 - (n*xbar + lambda_c*mu_c)^2/(n+lambda_c) + lambda_c*mu_c^2 )
	out = -nm*log2pi + 0.5*log(lambda_c/(lambda_c+n)) + lgamma(a_c+nm) - lgamma(a_c) + a_c*log(b_c) + (-a_c-nm)*log(aux1)
	return lg ? out : exp(out)
end


X_jt = rand(80)
mu_c = 0
lambda_c = 1
a_c = 2
b_c = 2

@btime similarity4(X_jt,mu_c,lambda_c,a_c,b_c,lg=true)
@btime similarity4_v2(X_jt,mu_c,lambda_c,a_c,b_c,lg=true)

similarity4_v2(X_jt,mu_c,lambda_c,a_c,b_c,lg=true)
# right one: 26.18445410365393
# with simd: 26.184454103653934
#  fastmath: 26.184454103653934

# ProfileCanvas.@profview_allocs similarity4(X_jt,mu_c,lambda_c,a_c,b_c,lg=true)
# ProfileCanvas.@profview begin 
# ProfileCanvas.@profview_allocs begin 
# @btime begin
# 	res = 0
# 	for _ in 1:100
# 		idxs = rand(collect(1:60),10)
# 		Xn = @view X_jt[idxs]
# 		res += similarity4(Xn,mu_c,lambda_c,a_c,b_c,lg=true)
# 	end
# end