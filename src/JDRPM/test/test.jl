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

