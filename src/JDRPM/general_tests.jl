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

