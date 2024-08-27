using LinearAlgebra
using BenchmarkTools

n = 20
s1 = rand(n)
s2 = rand(n)
alpha = 0.13

function cohesion1(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, alpha::Real; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	if sdim==1 
		return lg ? log(M) : M
	end
	# out = log(M) + lgamma(sdim)
	out = 0.0
	# compute the centroids
	cent1 = mean(s1)
	cent2 = mean(s2)
	# compute the sum of the distances (the D_h in the paper)
	sum_dist = sum(sqrt((s1[i] - cent1)^2 + (s2[i] - cent2)^2) for i in 1:sdim)

	# decide what to return
	if sum_dist >= 1
		out -= lgamma(alpha*sum_dist) # minus since we are at the denominator
	elseif sum_dist != 0
		out -= log(sum_dist) # minus since we are at the denominator
	else
		out = log(1)
	end
	return lg ? out : exp(out)
end


@btime cohesion1(s1,s2,alpha,lg=true)
