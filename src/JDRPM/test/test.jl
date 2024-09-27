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



using BenchmarkTools
nh_tmp = rand(100)
@btime sum(nh_tmp .> 0)
@btime count(x->(x>0), nh_tmp)


n = 100; rho_tmp = rand((1:5),n); k = 1
@btime findall(j -> rho_tmp[j]==k, 1:n) 
@btime findall(rho_tmp .== k)
@btime findall_faster(j -> rho_tmp[j]==k, 1:n)

# v = rand(20)
# @btime begin
# 	extra = extrema(v)
# 	R = extra[2] - extra[1]
# end
# @btime begin
# 	R = maximum(v) - minimum(v)
# end


function Gdist(a::Real,b::Real,R::Real)
	a==b && return 1
	return 1-abs(a-b)/R
end
function Gdist(a::String,b::String,R::Real)
	a==b && return 1
	return 0
end
function gower_sim(X_jt::Vector{<:Real})
	S = 0.
	isa(first(X_jt),Real) ? R = maximum(X_jt) - minimum(X_jt) : R = 0.
	cdim = length(v)
	for i in 1:(cdim-1)
		for j in (i+1):cdim
			S += Gdist(X_jt[i],X_jt[j],R)
		end
	end
	return S / (2*(cvdim*(cvdim-1))) # scale to [0,1]

end
v = [1,1000,2,1,0,-1,15,1,1,2] + rand(10)/10
v = rand(10)
v = repeat([1],10)
gower_sim(v)


include("../utils.jl")
lg = true
X_jt = [1,1,1]
similarity2(X_jt,1,lg)
similarity3(X_jt,1,lg)
XX_jt = [1,1,2]
similarity2(XX_jt,1,lg)
similarity3(XX_jt,1,lg)

similarity3([0],2,lg)


X_jt = [repeat([1],5)...,2]
similarity2(X_jt,1,lg)
similarity3(X_jt,1,lg)
similarity3(X_jt,2,lg)
X_jt = collect(1:201)
similarity2(X_jt,1,lg) # 0 ie completely dissimilar
similarity3(X_jt,1,lg) # strangely high, better to increase alpha
similarity3(X_jt,2,lg)
similarity3(X_jt,10,lg)

XX_jt = [1,1,1,1,1,1,7,8,9,10,11]
similarity2(XX_jt,1,lg)
similarity3(XX_jt,3,lg)
XX_jt = [1,2,3,4,5,6,7,8,9,10,11]
similarity2(XX_jt,1,lg)
similarity3(XX_jt,3,lg)

