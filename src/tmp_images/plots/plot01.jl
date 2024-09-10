using LinearAlgebra
using BenchmarkTools
using Plots

include("../utils.jl")

function benchmark_functions()
    results = Dict(
        "n" => Int[],
        "time_cohes3_scal" => Float64[],
        "allocs_cohes3_scal" => Int[],
        "memory_cohes3_scal" => Float64[],
        "time_cohes3_vect" => Float64[],
        "allocs_cohes3_vect" => Int[],
        "memory_cohes3_vect" => Float64[]
    )

	mu_0 = [0.,0.]
	k0 = 0.1
	v0 = 0.1
	Psi = [1. 0; 0 1.]
    
    for n in [5,10,20,50,100,250,500,1000]
		print("analysing case n=$n\r")
        s1 = rand(n)
        s2 = rand(n)
        
        # Benchmark cohes3
        benchmark1 = @benchmark cohesion3($s1, $s2, $mu_0, $k0, $v0, $Psi, lg=true)
        # benchmark1 = @benchmark cohesion3(s1, s2, mu_0, k0, v0, Psi, lg=true)
        time_cohes3 = minimum(benchmark1).time / 1e6  # Convert to milliseconds
        allocs_cohes3 = minimum(benchmark1).allocs
        memory_cohes3 = minimum(benchmark1).memory / 1024  # Convert to KB
        
        # Benchmark cohes3_vect
        benchmark2 = @benchmark cohesion3_4($s1, $s2, $mu_0, $k0, $v0, $Psi, Cohesion=3, lg=true)
        # benchmark2 = @benchmark cohesion3_4(s1, s2, mu_0, k0, v0, Psi, Cohesion=3, lg=true)
        time_cohes3_vect = minimum(benchmark2).time / 1e6  # Convert to milliseconds
        allocs_cohes3_vect = minimum(benchmark2).allocs
        memory_cohes3_vect = minimum(benchmark2).memory / 1024  # Convert to KB
        
        # Store results
        push!(results["n"], n)
        push!(results["time_cohes3_scal"], time_cohes3)
        push!(results["allocs_cohes3_scal"], allocs_cohes3)
        push!(results["memory_cohes3_scal"], memory_cohes3)
        push!(results["time_cohes3_vect"], time_cohes3_vect)
        push!(results["allocs_cohes3_vect"], allocs_cohes3_vect)
        push!(results["memory_cohes3_vect"], memory_cohes3_vect)
    end
    
    return results
end

results = benchmark_functions()

# Plot execution time
plot(results["n"], results["time_cohes3_scal"], label="scalar form", xlabel="n", ylabel="Time (ms)", title="Execution Time",
# yscale=:log10,
marker=:circle
)
plot!(results["n"], results["time_cohes3_vect"], label="vector form",
# yscale=:log10,
marker=:circle
)
savefig("execution_time.pdf")

# Plot memory allocations
plot(results["n"], results["allocs_cohes3_scal"], label="scalar form", xlabel="n", ylabel="Allocations", title="Memory Allocations",
# yscale=:log10,
marker=:circle
)
plot!(results["n"], results["allocs_cohes3_vect"], label="vector form",
# yscale=:log10,
marker=:circle
)
savefig("memory_allocations.pdf")

# Plot memory usage
plot(results["n"], results["memory_cohes3_scal"], label="scalar form", xlabel="n", ylabel="Memory (KB)", title="Memory Usage",
marker=:circle
)
plot!(results["n"], results["memory_cohes3_vect"], label="vector form",
marker=:circle
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
ns = [5,10,20,50,100, 250, 500, 1000] 
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
plot(ns, results["time_test2"], label=test2_label, marker=:circle,markersize=3,yscale=:log10)
plot!(ns, results["time_test3"], label=test3_label, marker=:circle,markersize=3,yscale=:log10)
plot!(ns, results["time_test4"], label=test4_label, marker=:circle,markersize=3,yscale=:log10)
plot!(xlabel="n", ylabel="Time (ms)", title="Execution Time", legend=:bottomright,
      ns, results["time_test1"], label=test1_label, marker=:circle,markersize=3,yscale=:log10)
savefig("execution_time_tests_sims.pdf")

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