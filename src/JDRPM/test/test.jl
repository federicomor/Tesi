nh_tmp = rand((-2:2),10)
using BenchmarkTools

function method1(nh_tmp)
	nclus_temp = count(x->(x>0),nh_tmp)
	return nclus_temp
end
function method2(nh_tmp)
	nclus_temp = sum(nh_tmp .> 0)
	return nclus_temp
end

method2(nh_tmp)
method2(nh_tmp)

@btime method2(nh_tmp)
@btime method1(nh_tmp)

using Random
Random.seed!(1)
function f()
lg_weights = rand(10)
max_ph = maximum(lg_weights)
sum_ph = 0.0
# exponentiate...
@simd for k in eachindex(lg_weights)
	 # for numerical purposes we subract max_ph
	lg_weights[k] = exp(lg_weights[k] - max_ph)
	sum_ph += lg_weights[k]
end
# sum_ph
# lg_weights
end

@btime f()