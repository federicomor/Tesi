- lgamma from SpecialFunctions is faster than definining a function as
```julia
function log_factorial(n::Int64)
	return sum(log.(i for i in 1:n))
end
```