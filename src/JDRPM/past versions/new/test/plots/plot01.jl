using LinearAlgebra
using BenchmarkTools
using Plots

using SpecialFunctions: logabsgamma
using Statistics
using StaticArrays
using LinearAlgebra

function lgamma(x::Real)
    first(logabsgamma(x))
end

logit(x::Real) = log(x / (one(x) - x))
const logpi = log(π)
const log2pi = log(2*π)

function G2a(a::Real, lg::Bool)
	out = logpi + lgamma(a) + lgamma(a - 0.5)
	return lg ? out : exp(out)
end

function cohesion4_v1(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}; lg::Bool, M::Real=1.0)::Float64
	sdim = length(s1)
	sp = [s1 s2]
	sbar =[mean(s1),mean(s2)]
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
	
	# for cohesion 4
	knn = kn+sdim
	vnn = vn+sdim
	mu_n = (k0*mu_0 + sdim*sbar)/(k0+sdim)
	Psi_nn = Psi_n + S + (kn*sdim)/(kn+sdim)*(sbar-mu_n)*(sbar-mu_n)'
	
	out = -sdim * logpi + G2a(0.5 * vnn, true) - G2a(0.5 * vn, true) + 0.5 * vn * logdet(Psi_n) - 0.5 * vnn * logdet(Psi_nn) + log(kn) - log(knn)
	return lg ? out : exp(out)
end

# scalar old version
function cohesion4_v2(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}; lg::Bool, M::Real=1.0)::Float64
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

# new in between one equally faster but tidier
function cohesion4_v3(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}; lg::Bool, M::Real=1.0, S=@MMatrix zeros(2, 2))::Float64
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
     # (You could probably also stack s1 and s2 into a matrix, and use row/columnwise mean to get sbar directly.)
    auxvec1 = sbar .- mu_0
    auxmat1 = auxvec1 * auxvec1'

    auxconst1 = k0 * sdim
    auxconst2 = k0 + sdim
    Psi_n = Psi .+ S .+ auxconst1 / (auxconst2) .* auxmat1

	# Updated parameters for cohesion 4
	knn = kn + sdim
	vnn = vn + sdim

	mu_n = (k0 * mu_0 + sdim * sbar) / (auxconst2)
	sbar_mun = sbar - mu_n
	auxmat2 = sbar_mun * sbar_mun' 

	auxconst3 = kn * sdim
	auxconst4 = kn + sdim
	Psi_nn = Psi_n .+ S .+ auxconst3 / (auxconst4) .* auxmat2

	out = -sdim * logpi + G2a(0.5 * vnn, true) - G2a(0.5 * vn, true) + 0.5 * vn * logdet(Psi_n) - 0.5 * vnn * logdet(Psi_nn) + log(kn) - log(knn)
	return lg ? out : exp(out)
end

function benchmark_functions()
    results = Dict(
        "n" => Int[],
        "time_cohes3_v2" => Float64[],
        "allocs_cohes3_v2" => Int[],
        "memory_cohes3_v2" => Float64[],
        "time_cohes3_v1" => Float64[],
        "allocs_cohes3_v1" => Int[],
        "memory_cohes3_v1" => Float64[],
		"time_cohes3_v3" => Float64[],
        "allocs_cohes3_v3" => Int[],
        "memory_cohes3_v3" => Float64[]
    )

	mu_0 = [0.,0.]
	k0 = 1.
	v0 = 5.
	Psi = [1. 0.1; 0.1 1.]	
	mu_0_static = @SVector [0.,0.]
	Psi_static = @SMatrix [1. 0.1; 0.1 1.]

    for n in [5,10,20,50,100,250,500,1000]
		print("analysing case n=$n\r")
        s1 = rand(n)
        s2 = rand(n)
        
        benchmark1 = @benchmark cohesion4_v1($s1, $s2, $mu_0, $k0, $v0, $Psi, lg=true)
        # @btime cohesion4_v1($s1, $s2, $mu_0, $k0, $v0, $Psi, lg=true)
		# display(benchmark1)
        time_cohes3_v1 = minimum(benchmark1).time / 1e6  # Convert to milliseconds
        allocs_cohes3_v1 = minimum(benchmark1).allocs
        memory_cohes3_v1 = minimum(benchmark1).memory / 1024  # Convert to KB

        benchmark2 = @benchmark cohesion4_v2($s1, $s2, $mu_0, $k0, $v0, $Psi, lg=true)
        # @btime cohesion4_v2($s1, $s2, $mu_0, $k0, $v0, $Psi, lg=true)
		# display(benchmark2)
        time_cohes3_v2 = minimum(benchmark2).time / 1e6  # Convert to milliseconds
        allocs_cohes3_v2 = minimum(benchmark2).allocs
        memory_cohes3_v2 = minimum(benchmark2).memory / 1024  # Convert to KB

		S=@MMatrix zeros(2, 2)
		benchmark3 = @benchmark cohesion4_v3($s1, $s2, $mu_0_static, $k0, $v0, $Psi_static, lg=true, S=$S)
		# display(benchmark3)
        time_cohes3_v3 = minimum(benchmark3).time / 1e6  # Convert to milliseconds
        allocs_cohes3_v3 = minimum(benchmark3).allocs
        memory_cohes3_v3 = minimum(benchmark3).memory / 1024  # Convert to KB
        
        # Store results
        push!(results["n"], n)
        push!(results["time_cohes3_v2"], time_cohes3_v2)
        push!(results["allocs_cohes3_v2"], allocs_cohes3_v2)
        push!(results["memory_cohes3_v2"], memory_cohes3_v2)

        push!(results["time_cohes3_v1"], time_cohes3_v1)
        push!(results["allocs_cohes3_v1"], allocs_cohes3_v1)
        push!(results["memory_cohes3_v1"], memory_cohes3_v1)

		push!(results["time_cohes3_v3"], time_cohes3_v3)
        push!(results["allocs_cohes3_v3"], allocs_cohes3_v3)
        push!(results["memory_cohes3_v3"], memory_cohes3_v3)
    end
    
    return results
end
using Random
Random.seed!(4)
results = benchmark_functions()

# Plot execution time
plot(results["n"], results["time_cohes3_v2"], label="scalar form", xlabel="n", ylabel="Time (ms)", title="Execution Time",legend=:bottomright,
yscale=:log10,
marker=:circle,markersize=3
)
plot!(results["n"], results["time_cohes3_v1"], label="vector form",
yscale=:log10,
marker=:circle,markersize=3
)
plot!(results["n"], results["time_cohes3_v3"], label="static form",
yscale=:log10,
marker=:circle,markersize=3
)
savefig("execution_time.pdf")

# Plot memory allocations
plot(results["n"], results["allocs_cohes3_v2"], label="scalar form", xlabel="n", ylabel="Allocations", title="Memory Allocations",
# yscale=:log10,
marker=:circle,markersize=3
)
plot!(results["n"], results["allocs_cohes3_v1"], label="vector form",
# yscale=:log10,
marker=:circle,markersize=3
)
plot!(results["n"], results["allocs_cohes3_v3"], label="static form",
# yscale=:log10,
marker=:circle,markersize=3
)
savefig("memory_allocations.pdf")

# Plot memory usage
plot(results["n"], results["memory_cohes3_v2"], label="scalar form", xlabel="n", ylabel="Memory (KB)", title="Memory Usage",
marker=:circle,markersize=3
)
plot!(results["n"], results["memory_cohes3_v1"], label="vector form",
marker=:circle,markersize=3
)
plot!(results["n"], results["memory_cohes3_v3"], label="static form",
marker=:circle,markersize=3
)
savefig("memory_usage.pdf")


#########################################

X_jt = rand(100)
function test1(X_jt::Vector{Float64})
	aux2 = 0.0
	@inbounds @fastmath @simd for i in eachindex(X_jt)
		aux2 += X_jt[i]^2
	end
	return aux2
end
function test2(X_jt::Vector{Float64})
	aux2 = 0.0
	for i in eachindex(X_jt)
		aux2 += X_jt[i]^2
	end
	return aux2
end
function test3(X_jt::Vector{Float64})
	aux2 = 0.0
	@inbounds for i in eachindex(X_jt)
		@inbounds aux2 += X_jt[i]^2
	end
	return aux2
end
function test4(X_jt::Vector{Float64})
	aux2 = 0.0
	@inbounds @simd for i in eachindex(X_jt)
		aux2 += X_jt[i]^2
	end
	return aux2
end

test1_label = "@inbounds @simd @fastmath"
test2_label = "no annotations"
test3_label = "@inbounds"
test4_label = "@inbounds @simd"


test1(X_jt)
test2(X_jt)
test3(X_jt)
test4(X_jt)

using BenchmarkTools
using Plots
using Random
Random.seed!(42)
ns = [4,5,8,10,16,20,32,50,64,100,128, 250, 256, 500, 512, 1000] 
function benchmark_tests(ns)
    results = Dict(
        "time_test1" => Float64[],
        "allocs_test1" => Int[],
        "memory_test1" => Float64[],
        "time_test2" => Float64[],
        "allocs_test2" => Int[],
        "memory_test2" => Float64[],
        "time_test3" => Float64[],
        "allocs_test3" => Int[],
        "memory_test3" => Float64[],
        "time_test4" => Float64[],
        "allocs_test4" => Int[],
        "memory_test4" => Float64[]
    )

    for n in ns # Vector sizes to benchmark
        println("Benchmarking with n = $n")
        X_jt = rand(n)

        # Benchmark test1
        bench1 = @benchmark test1($X_jt)
		println("\t",test1(X_jt))
        push!(results["time_test1"], minimum(bench1).time / 1e6)
        push!(results["allocs_test1"], minimum(bench1).allocs)
        push!(results["memory_test1"], minimum(bench1).memory / 1024)

        # Benchmark test2
        bench2 = @benchmark test2($X_jt)
		println("\t",test2(X_jt))
        push!(results["time_test2"], minimum(bench2).time / 1e6)
        push!(results["allocs_test2"], minimum(bench2).allocs)
        push!(results["memory_test2"], minimum(bench2).memory / 1024)

        # Benchmark test3
        bench3 = @benchmark test3($X_jt)
		println("\t",test3(X_jt))
        push!(results["time_test3"], minimum(bench3).time / 1e6)
        push!(results["allocs_test3"], minimum(bench3).allocs)
        push!(results["memory_test3"], minimum(bench3).memory / 1024)

        # Benchmark test4
        bench4 = @benchmark test4($X_jt)
		println("\t",test4(X_jt))
        push!(results["time_test4"], minimum(bench4).time / 1e6)
        push!(results["allocs_test4"], minimum(bench4).allocs)
        push!(results["memory_test4"], minimum(bench4).memory / 1024)
    end

    return results
end

results = benchmark_tests(ns)

# Plot execution time
msize=2
plot(ns, results["time_test2"], label=test2_label, marker=:circle,markersize=msize,yscale=:log10)
plot!(ns, results["time_test3"], label=test3_label, marker=:circle,markersize=msize,yscale=:log10)
plot!(ns, results["time_test4"], label=test4_label, marker=:circle,markersize=msize,yscale=:log10)
plot!(xlabel="n", ylabel="Time (ms)", title="Execution Time", legend=:bottomright,
      ns, results["time_test1"], label=test1_label, marker=:circle,markersize=msize,yscale=:log10)
savefig("execution_time_tests_sims_pow2.pdf")

# Plot memory allocations
plot(ns, results["allocs_test1"], label="test1", xlabel="n", ylabel="Allocations", title="Memory Allocations", marker=:circle)
plot!(ns, results["allocs_test2"], label="test2", marker=:circle)
plot!(ns, results["allocs_test3"], label="test3", marker=:circle)
plot!(ns, results["allocs_test4"], label="test4", marker=:circle)
savefig("memory_allocations_tests.pdf")

# Plot memory usage
plot(ns, results["memory_test1"], label="test1", xlabel="n", ylabel="Memory (KB)", title="Memory Usage", marker=:circle)
plot!(ns, results["memory_test2"], label="test2", marker=:circle)
plot!(ns, results["memory_test3"], label="test3", marker=:circle)
plot!(ns, results["memory_test4"], label="test4", marker=:circle)
savefig("memory_usage_tests.pdf")






#########################################
using Random
using ProfileCanvas

mu_0 = [0.,0.]
k0 = 0.1
v0 = 0.1
Psi = [1. 0; 0 1.]

n = 200

Random.seed!(1)
@timed begin
res = 0
for _ in 1:10_000
	s1 = rand(n)
	s2 = rand(n)
	res += cohesion3(s1, s2, mu_0, k0, v0, Psi, lg=true)
end
end
res

Random.seed!(1)
@timed begin
	res = 0
	for _ in 1:10_000
		s1 = rand(n)
		s2 = rand(n)
		res += cohesion3_4(s1, s2, mu_0, k0, v0, Psi, Cohesion=3, lg=true)
	end
end
res

bout = @benchmark cohesion3(s1, s2, mu_0, k0, v0, Psi, lg=true)