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