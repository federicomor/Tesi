using Random
using BenchmarkTools

@show Threads.nthreads()

Random.seed!(234)

N = 20
a = collect(1:N)
b = zeros(N)
c = zeros(N)
@btime Threads.@threads for i in 1:N
	b[i] = a[i]^2
end

@btime for i in 1:N
	c[i] = a[i]^2
end

@assert b == c


# T = 24
# theta_iter = rand(T)
# phi1_iter = 1.32
# # classical method
# aux2 = 0.0
# @btime for t in 2:T
# 	aux2 += theta_iter[t] - phi1_iter*theta_iter[t-1]
# end

##########################################

using Distributions
using Statistics

p = 4
s2_beta = 1
n = 10
T = 5
Si_iter = rand((1:5),n,T)
muh_iter = rand(n,T)
sig2h_iter = rand(n,T)
Xlk_covariates = rand(n,p,T)
Y = rand(n,T)

sum_Y = zeros(p)
A_star = I(p)/s2_beta
for j in 1:n
	label_j = Si_iter[j,t]
	X_jt = @view Xlk_covariates[j,:,t]
	sum_Y += (Y[j,t] - muh_iter[label_j,t]) * X_jt / sig2h_iter[label_j,t]
	A_star += (X_jt * X_jt') / sig2h_iter[label_j,t]
end
# debug(@showd A_star)
b_star = beta0/s2_beta + sum_Y
# debug(@showd b_star)
beta_iter[t] = rand(MvNormal(b_star, inv(Symmetric(A_star))))



##########################################


p = 5
mu = rand(p)
A = rand(p,p) # auxiliary for sigma construction
sigma = Symmetric(A * A')
isigma = Symmetric(inv(sigma))

Random.seed!(34)
println(rand(MvNormal(sigma*mu, sigma)))
Random.seed!(34)
println(rand(MvNormalCanon(mu, isigma)))

######################################

function similarity1(X_jt::Union{Vector{Float64}, Vector{String}}, alpha::Real; lg::Bool)
	# eg Xstar_jt would be X_cl[indexes,covariate_index,t] 

	if isa(first(X_jt),Float64) # numerical
		xbar_j = mean(X_jt)
		card_Sjt = length(X_jt)
		Hx = 0.0
		for i in 1:card_Sjt
			Hx += (X_jt[i] - xbar_j)^2 
		end
		Hx /= card_Sjt
		
		return lg ? -alpha * Hx : exp(-alpha * Hx) 

	else # categorical
		# X_jt = ["Old", "Young", "Middle", "Young", "Old", "Old"]
		unique_keys = unique(X_jt)
		# @show unique_keys
		counts = Dict(key => 0.0 for key in unique_keys)
		for item in X_jt
		    counts[item] += 1
		end
		# @show counts
		total = length(X_jt)
		for key in keys(counts)
		    counts[key] = counts[key]/total
		end
		Hx = 0.0
		for key in keys(counts)
			Hx -= counts[key] * log(counts[key])
		end
		
		return lg ? -alpha * Hx : exp(-alpha * Hx) 
	end
end

similarity1([repeat(["a"],1000)...,"b"], 5, lg=false)


########################################à
# Gower dissimilarity
function gower_d(x1::Real, x2::Real, R::Real) # numerical case
	return abs(x1-x2)/R
end
function gower_d(x1::String, x2::String, R::Real) # categorical case
	# R is fictitious here, just to allow the same function call later
	return Int(!(x1 == x2))
end

# Total Gower dissimilarity - paper 6 pag 4
function similarity2(X_jt::Union{Vector{<:Real}, Vector{String}}, alpha::Real; lg::Bool)
	H = 0.0
	if isa(first(X_jt),Real) 
		R = (maximum(X_jt) - minimum(X_jt))
		if R==0.0 R=eps() end
		H /= R
	else R=0.0 end
	n_j = length(X_jt)
	for l in 1:(n_j-1)
		for k in (l+1):n_j
			H += gower_d(X_jt[l], X_jt[k], R)
		end
	end 

	out = -alpha * H
	return lg ? out : exp(out)
end

# Average Gower dissimilarity - paper 6 pag 4
function similarity3(X_jt::Union{Vector{<:Real}, Vector{String}}, alpha::Real; lg::Bool)
	H = 0.0
	if isa(first(X_jt),Real) 
		R = (maximum(X_jt) - minimum(X_jt))
		if R==0.0 R=eps() end
		H /= R
	else R=0.0 end
	n_j = length(X_jt)
	for l in 1:(n_j-1)
		for k in l:n_j
			H += gower_d(X_jt[l], X_jt[k], R)
		end
	end 
	out = -2*alpha / (n_j*(n_j-1)) * H
	return lg ? out : exp(out)
end

X_jt = rand(Int8,20)
X_jt = [1]
X_jt = [1,1,1,1,1,2]
X_jt = [repeat([1],200)...,2]

# 0 => completely dissimilar
# 1 => completely similar

X_jt = [1,1,1]
similarity2(X_jt,1,lg=false)
similarity3(X_jt,1,lg=false)
XX_jt = [1,1,2]
similarity2(XX_jt,1,lg=false)
similarity3(XX_jt,1,lg=false)

X_jt = [repeat([1],5)...,2]
similarity2(X_jt,1,lg=false)
similarity3(X_jt,1,lg=false)
similarity3(X_jt,2,lg=false)
X_jt = collect(1:201)
similarity2(X_jt,1,lg=false) # 0 ie completely dissimilar
similarity3(X_jt,1,lg=false) # strangely high, better to increase alpha
similarity3(X_jt,2,lg=false)
similarity3(X_jt,10,lg=false)

XX_jt = [1,1,1,1,1,1,7,8,9,10,11]
similarity2(XX_jt,1,lg=false)
similarity3(XX_jt,3,lg=false)
XX_jt = [1,2,3,4,5,6,7,8,9,10,11]
similarity2(XX_jt,1,lg=false)
similarity3(XX_jt,3,lg=false)


################################################

function auxiliary_sim(X_jt::AbstractVector{<:Real}, mu_c::Real, lambda_c::Real, a_c::Real, b_c::Real; lg::Bool)
	n = length(X_jt)
	nm = n/2
	xbar = mean(X_jt)
	aux1 = b_c + 0.5 * (sum(X_jt .^ 2) - (n*xbar + lambda_c*mu_c)^2/(n+lambda_c) + lambda_c*mu_c^2 )
	# @show aux1
	out = -nm*log2pi + 0.5*log(lambda_c/(lambda_c+n)) + lgamma(a_c+nm) - lgamma(a_c) + a_c*log(b_c) + (-a_c-nm)*log(aux1)
	return lg ? out : exp(out)
end



mu_c = 0.0
lambda_c = 1.0
a_c = 2.0
b_c = 2.0

X_jt = rand(10)*2
X_jt = [1,1,2,2,2,1,2,1,1,1] ./ 10
auxiliary_sim(X_jt,mu_c,lambda_c,a_c,b_c,lg=true)
X_jt = [1,1,2,2,2,1,2,1,1,100] ./ 10
auxiliary_sim(X_jt,mu_c,lambda_c,a_c,b_c,lg=true)


############################################
include("utils.jl")

cv_params = [mu_c,lambda_c,a_c,b_c]

X_jt = [1,1,2,2,2,1,2,1,1,1] ./ 10
@btime covariate_similarity(4,X_jt,cv_params,lg=true)

covariate_similarity(4,[1],cv_params,lg=true)

X_jt = [1,1,2,2,2,1,2,1,1,100] ./ 10
covariate_similarity(4,X_jt,cv_params,lg=true)

covariate_similarity(1,[repeat(["a"],100000)...,"b"], [2], lg=true)
covariate_similarity(1,[repeat(["a"],10)...,"b"], [2], lg=true)

##############################################


# paper 3 section 3.1
function cohesion3_4(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}; Cohesion::Int, lg::Bool, M::Real=1.0, S=@MMatrix zeros(2, 2))::Float64
	sdim = length(s1)
	sp = [s1 s2]
	sbar = SVector(mean(s1),mean(s2))
	# S = sum( (sp[i,:] - sbar)*(sp[i,:] - sbar)' for i in 1:sdim) # sum is slow
	# FIXED: sum is slow because i didnt initialize S
	# S = zeros(2,2); S = sum( (sp[i,:] - sbar)*(sp[i,:] - sbar)' for i in 1:sdim) # this is fast now
	# S = @MMatrix zeros(2,2)
	S .= 0.
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

# below there is the one "translated" from C, which is way more efficient
function cohesion3_old(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	# Compute sample means
	sbar1 = mean(s1)
	sbar2 = mean(s2)
	# Compute deviations from the sample mean
	S1, S2, S3, S4 = 0.0, 0.0, 0.0, 0.0
	@inbounds for i in 1:sdim
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

# here the in between version, which is faster
function cohesion3(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}, lg::Bool, M::Real=1.0, S=@MMatrix zeros(2, 2))::Float64
	sdim = length(s1)
	# Compute sample means
	sbar1 = mean(s1)
	sbar2 = mean(s2)
	# Compute deviations from the sample mean
	S .= 0.
	@inbounds for i in 1:sdim
		s_sbar1 = s1[i] - sbar1
		s_sbar2 = s2[i] - sbar2

		S[1, 1] += s_sbar1 * s_sbar1
		S[2, 2] += s_sbar2 * s_sbar2
		S[2, 1] += s_sbar1 * s_sbar2
	end
	S[1, 2] = S[2, 1] # to avoid repeating computations
	# Updated parameters for cohesion 3
	kn = k0 + sdim
	vn = v0 + sdim

	sbar = SVector((sbar1, sbar2))
	auxvec1 = sbar .- mu_0
	auxmat1 = auxvec1 * auxvec1'

	auxconst1 = k0 * sdim
	auxconst2 = k0 + sdim
	Psi_n = Psi .+ S .+ auxconst1 / (auxconst2) .* auxmat1
	
	out = -sdim * logpi + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * logdet(Psi) - 0.5 * vn * logdet(Psi_n) + log(k0) - log(kn)
	return lg ? out : exp(out)
end

function cohesion3!(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}, lg::Bool,  M::Real=1.0, S=@MMatrix(zeros(2, 2)), case::Int=1, add::Bool=false, lC=@MVector(zeros(2)))::Float64
	sdim = length(s1)
	# Compute sample means
	sbar1 = mean(s1)
	sbar2 = mean(s2)
	# Compute deviations from the sample mean
	S .= 0.
	@inbounds for i in 1:sdim
		s_sbar1 = s1[i] - sbar1
		s_sbar2 = s2[i] - sbar2

		S[1, 1] += s_sbar1 * s_sbar1
		S[2, 2] += s_sbar2 * s_sbar2
		S[2, 1] += s_sbar1 * s_sbar2
	end
	S[1, 2] = S[2, 1] # to avoid repeating computations
	# Updated parameters for cohesion 3
	kn = k0 + sdim
	vn = v0 + sdim

	sbar = SVector((sbar1, sbar2))
	auxvec1 = sbar .- mu_0
	auxmat1 = auxvec1 * auxvec1'

	auxconst1 = k0 * sdim
	auxconst2 = k0 + sdim
	Psi_n = Psi .+ S .+ auxconst1 / (auxconst2) .* auxmat1
	
	out = -sdim * logpi + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * logdet(Psi) - 0.5 * vn * logdet(Psi_n) + log(k0) - log(kn)
	if add
		lC[case] += lg ? out : exp(out)
	else
		lC[case] = lg ? out : exp(out)
	end
end

include("../utils.jl")
# s1 = rand(10)
# s2 = rand(10)
s1 = [1.0]
s2 = [1.1]
mu_0 = SVector((0.,0.))
k0 = 1.
v0 = 5.
Psi = @SMatrix [1. 0;0 1]
M = 1.0
S= @MMatrix(zeros(2,2))

cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=3, lg=true, M=M, S=S)

cohesion3_old(s1, s2, mu_0, k0, v0, Psi; lg=true, M=M)

cohesion3(s1, s2, mu_0, k0, v0, Psi, lg, M, S)

lC = @MVector zeros(2)
case = 1; add = false; 
lC[case]
cohesion3!(s1, s2, mu_0, k0, v0, Psi, true, M, S, case, add, lC)
lC[case]
