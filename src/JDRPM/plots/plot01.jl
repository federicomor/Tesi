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
        time_cohes3 = minimum(benchmark1).time / 1e6  # Convert to milliseconds
        allocs_cohes3 = minimum(benchmark1).allocs
        memory_cohes3 = minimum(benchmark1).memory / 1024  # Convert to KB
        
        # Benchmark cohes3_vect
        benchmark2 = @benchmark cohesion3_4($s1, $s2, $mu_0, $k0, $v0, $Psi, Cohesion=3, lg=true)
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
yscale=:log10,
marker=:circle
)
plot!(results["n"], results["time_cohes3_vect"], label="vector form",
yscale=:log10,
marker=:circle
)
savefig("execution_time.svg")

# Plot memory allocations
plot(results["n"], results["allocs_cohes3_scal"], label="scalar form", xlabel="n", ylabel="Allocations", title="Memory Allocations",
# yscale=:log10,
marker=:circle
)
plot!(results["n"], results["allocs_cohes3_vect"], label="vector form",
# yscale=:log10,
marker=:circle
)
savefig("memory_allocations.svg")

# Plot memory usage
plot(results["n"], results["memory_cohes3_scal"], label="scalar form", xlabel="n", ylabel="Memory (KB)", title="Memory Usage",
marker=:circle
)
plot!(results["n"], results["memory_cohes3_vect"], label="vector form",
marker=:circle
)
savefig("memory_usage.svg")


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
