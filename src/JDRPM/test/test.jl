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

#######################################
using Random
Random.seed!(1)
sp1 = rand(10)
indexes = [1,2,3,4,6,8]

function f_test(indexes,j)
	sp1_red = sp1[indexes]
	# sp1_red = @view sp1[indexes]
	result = 0.
	for k in 1:10
		aux_idxs = randsubseq(1:length(sp1_red),0.5)
		s1o = @view sp1_red[aux_idxs]
		# copy!(s1o,sp1_red[aux_idxs])
		s1n = copy(s1o); push!(s1n, sp1[j])
		# println(typeof(s1n))
		# resize!(s1n,length(s1o)+1)
		# copy!(s1n,s1o)
		# s1n[end] = sp1[j]
		# copy!(s1n,s1o); push!(s1n, sp1[j])
		result += (sum(s1o)-sum(s1n))*k
	end
	return result
end
function f_test_prealloc(indexes,j)
	sp1_red = sp1[indexes]
	s1n = zeros(length(sp1))
	# sp1_red = @view sp1[indexes]
	result = 0.
	for k in 1:10
		aux_idxs = randsubseq(1:length(sp1_red),0.5)
		s1o = @view sp1_red[aux_idxs]
		copy!(s1n,s1o); push!(s1n, sp1[j])
		result += (sum(s1o)-sum(s1n))*k
	end
	return result
end
f_test(indexes,1)
f_test_prealloc(indexes,1)