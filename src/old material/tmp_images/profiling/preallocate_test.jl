using Random
Random.seed!(1)

function process_data(n)
    results = []
	Random.seed!(1)
    for i in 1:n
        data = rand(100)  # Allocating a new array inside the loop
        processed = data .^ 2  # Another allocation
        push!(results, processed)
    end
    return results
end


function process_data!(data::Vector{Float64}, processed::Vector{Float64}, results::Vector{Vector{Float64}}, n::Int)
	Random.seed!(1)
    for i in 1:n
        rand!(data)  # Reuse the preallocated 'data' array
        @. processed = data ^ 2  # Reuse the preallocated 'processed' array
        results[i] = copy(processed)  # Store the result
    end
end


n = 10_000
data = Vector{Float64}(undef, 100)        # Preallocate 'data' array
processed = Vector{Float64}(undef, 100)   # Preallocate 'processed' array
res2 = Vector{Vector{Float64}}(undef, n)  # Preallocate 'results' array

process_data(1) # compilation
println("---- Without preallocation")
@timev res1 = process_data(n);
println("---- With preallocation")
@timev process_data!(data, processed, res2, n)
println("Same results? ", res1 == res2)

