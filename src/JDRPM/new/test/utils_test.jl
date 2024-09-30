# Assuming all cohesion functions are already defined (cohesion1, cohesion2, etc.)
using Test
include("../utils.jl")

# Define helper function to run and compare different versions
function test_spatial_cohesion()
    # Test inputs
    s1 = rand(10)
    s2 = rand(10)
	mu0 = 0.
    sp_params_vector1256 = [0.5]
    sp_params_vector34 = [[mu0,mu0],0.2,0.4,[1. 0.0; 0.0 1.]]  # Example vector
    sp_params_struct = SpParams(alpha=0.5, a=0.5, mu_0=[mu0,mu0], k0=0.2, v0=0.4, Psi=[1.0 0.0; 0.0 1.0], phi=0.5)
    lg = true
    M = 1.0
    S = @MMatrix zeros(2, 2)

    # Collect results for each function
	for i in 1:6
		println("\n[cohesion $i]")
		results_vector = 0
		results_struct = 0
		if i in [1,2,5,6]
			results_vector = spatial_cohesion(i, s1, s2, sp_params_vector1256, lg, M, S) 
			results_struct = spatial_cohesion(i, s1, s2, sp_params_struct, lg, M, S)
		else
			results_vector = spatial_cohesion(i, s1, s2, sp_params_vector34, lg, M, S) 
			results_struct = spatial_cohesion(i, s1, s2, sp_params_struct, lg, M, S)
		end
		println("\nvec: $(results_vector)\nstr: $(results_struct)")
        # @test isapprox(results_vector[i], results_struct[i], atol=1e-6)
    end
    # println("All tests passed!")
end
test_spatial_cohesion()

# Helper function to test the in-place (mutating) versions
function test_spatial_cohesion_mutating()
    # Test inputs
    s1 = rand(10)
    s2 = rand(10)
	mu0 = 0.
    sp_params_vector1256 = [0.5]
    sp_params_vector34 = [[mu0,mu0],0.2,0.4,[1. 0.0; 0.0 1.]]  # Example vector
    sp_params_struct = SpParams(alpha=0.5, a=0.5, mu_0=[mu0,mu0], k0=0.2, v0=0.4, Psi=[1.0 0.0; 0.0 1.0], phi=0.5)
    lg = true
    M = 1.0
    S = @MMatrix zeros(2, 2)
	lC_v = @MVector zeros(2)
	lC_s = @MVector zeros(2)
	case = 1
	add = false

    # Collect results for each function
	for i in 1:6
		println("\n[cohesion $i]")
		results_vector = 0
		results_struct = 0
		if i in [1,2,5,6]
			spatial_cohesion!(i, s1, s2, sp_params_vector1256, lg, M, S, case,add,lC_v) 
			results_vector = lC_v
			spatial_cohesion!(i, s1, s2, sp_params_struct, lg, M, S, case,add,lC_s)
			results_struct = lC_s
		else
			spatial_cohesion!(i, s1, s2, sp_params_vector34, lg, M, S, case,add,lC_v) 
			results_vector = lC_v
			spatial_cohesion!(i, s1, s2, sp_params_struct, lg, M, S, case,add,lC_s)
			results_struct = lC_s
		end
		println("\nvec: $(results_vector)\nstr: $(results_struct)")
        @test isapprox(results_vector, results_struct, atol=1e-6)
    end
end

# Run the tests
test_spatial_cohesion()
test_spatial_cohesion_mutating()
