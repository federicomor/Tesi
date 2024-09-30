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
		# println("\n[cohesion $i]")
		results_vector = 0
		results_struct = 0
		if i in [1,2,5,6]
			results_vector = spatial_cohesion(i, s1, s2, sp_params_vector1256, lg, M, S) 
			results_struct = spatial_cohesion(i, s1, s2, sp_params_struct, lg, M, S)
		else
			results_vector = spatial_cohesion(i, s1, s2, sp_params_vector34, lg, M, S) 
			results_struct = spatial_cohesion(i, s1, s2, sp_params_struct, lg, M, S)
		end
		# println("\nvec: $(results_vector)\nstr: $(results_struct)")
        # @test isapprox(results_vector, results_struct, atol=1e-6)
        @test isequal(results_vector, results_struct)
    end
    println("All tests passed!")
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
		# println("\n[cohesion $i]")
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
		# println("\nvec: $(results_vector)\nstr: $(results_struct)")
        # @test isapprox(results_vector, results_struct, atol=1e-6)
        @test isequal(results_vector, results_struct)
    end
	print("All tests passed!")
end
test_spatial_cohesion_mutating()

# Helper function to compare mutating and non-mutating numerical versions
function test_covariate_similarity_numerical()
    # Test inputs
    X_jt_num = rand(5)  # Example numerical covariates
    cv_params_num = [0.5, 0.1, 0.3, 0.2]  # Example parameters
    lg = true
	add = false
	case = 1

    # Run and collect results for non-mutating numerical similarity function
    results_non_mutating = [covariate_similarity(idx, X_jt_num, cv_params_num, lg) for idx in 1:4]

    # Check if mutating versions produce the same results
    lS = @MVector zeros(2)
    for idx in 1:4
        # Reset lS for each idx
        lS .= 0.0
        covariate_similarity!(idx, X_jt_num, cv_params_num, lg, 1, false, lS)
        
        # Compare the results
        # @test isapprox(results_non_mutating[idx], lS[1], atol=1e-6)
        @test isequal(results_non_mutating[idx], lS[1])
        println("Numerical test passed for idx=$idx")
    end
end
test_covariate_similarity_numerical()

# Helper function to compare mutating and non-mutating categorical versions
function test_covariate_similarity_categorical()
    # Test inputs
    X_jt_cat = ["A", "B", "A", "C", "B"]  # Example categorical covariates
    cv_params_cat = [0.5]  # Example parameters
    lg = true

    # Run and collect results for non-mutating categorical similarity function
    results_non_mutating = [covariate_similarity(idx, X_jt_cat, cv_params_cat, lg) for idx in 1:3]

    # Check if mutating versions produce the same results
    lS = @MVector zeros(2)
    for idx in 1:3
        # Reset lS for each idx
        lS .= 0.0
        covariate_similarity!(idx, X_jt_cat, cv_params_cat, lg, 1, false, lS)
        
        # Compare the results
        # @test isapprox(results_non_mutating[idx], lS[1], atol=1e-6)
        @test isequal(results_non_mutating[idx], lS[1])
        println("Categorical test passed for idx=$idx")
    end
end
test_covariate_similarity_categorical()
