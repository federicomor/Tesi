using LinearAlgebra
using BenchmarkTools
include("../utils.jl")

function spatial_cohesion_splat(idx::Real, s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, sp_params::Vector; lg::Bool, M::Real)
	idx==1.0 && return cohesion1(s1,s2,sp_params...,lg=lg,M=M) 
	idx==2.0 && return cohesion2(s1,s2,sp_params...,lg=lg,M=M) 
	idx==3.0 && return cohesion3(s1,s2,sp_params...,lg=lg,M=M) 
	idx==4.0 && return cohesion3(s1,s2,sp_params...,lg=lg,M=M) 
	idx==5.0 && return cohesion5(s1,s2,sp_params...,lg=lg,M=M) 
	idx==6.0 && return cohesion6(s1,s2,sp_params...,lg=lg,M=M) 
end

n = 20
s1 = rand(n)
s2 = rand(n)
alpha = 0.13
a = 10
mu_0 = [0.,0.]
k0 = 0.1
v0 = 0.1
Psi = [1. 0.2; 0.2 1]
phi=0.5

sp_params = [mu_0, k0, v0, Psi]

@btime spatial_cohesion(4,s1,s2,sp_params,lg=true,M=1)
@btime spatial_cohesion_splat(4,s1,s2,sp_params,lg=true,M=1)

@code_llvm spatial_cohesion(4,s1,s2,sp_params,lg=true,M=1)
@code_llvm spatial_cohesion_splat(4,s1,s2,sp_params,lg=true,M=1)


cohesion1(s1,s2,alpha,lg=true)
cohesion2(s1,s2,a,lg=true)

cohesion3(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi,Cohesion=3,lg=true)
cohesion4(s1, s2, mu_0, k0, v0, Psi,lg=true)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi,Cohesion=4,lg=true)

@timev cohesion5(s1,s2,phi,lg=true)
@timev cohesion6(s1,s2,phi,lg=true)


#############################################


using Plots, Distributions

# Define the target density π(x) as a mixture of two normal distributions
pi_density(x) = 0.5 * pdf(Normal(-1, 0.7), x) + 0.5 * pdf(Normal(1, 0.5), x)

# Define the proposal density q(x) as a normal distribution N(0, 2)
q = Normal(0, 2)
q_density(x) = pdf(q, x)

# Constant c such that c*q(x) covers π(x)
c = 2.8  # Adjust based on trial and error to ensure c*q(x) is above π(x)

# Set up the plot range
x_range = -4:0.01:4

# Initialize list to store accepted samples
accepted_samples = []

# Create animation using @animate macro
anim = @animate for i in 1:10
	print("$i\r")
    # Step 2: Draw X from q(x)
    X = rand(q)
    
    # Step 3: Draw U from a uniform between [0, c*q(X)]
    U = rand() * c * q_density(X)
    
    # Step 4: Accept/reject
    if U <= pi_density(X)
        push!(accepted_samples, X)  # Store the accepted sample
        accepted = true
        status = "Accepted"
        color = :green
    else
        accepted = false
        status = "Rejected"
        color = :red
    end
    
    # Plot the target and proposal distributions
    p = plot(x_range, pi_density.(x_range), label="Target Density π(x)", lw=2,
	# legend=:outerright)
	)
    plot!(x_range, c*q_density.(x_range), label="Scaled Proposal c*q(x)", lw=2, linestyle=:dash)

    # Plot all accepted samples
    scatter!([accepted_samples...], [0.0 for _ in accepted_samples], 
             label="", color=:black, marker=:circle, markersize=3)
            #  label="Accepted Samples", color=:black, marker=:circle, markersize=8)
    
    # Plot the current sampled point
    scatter!([X], [0], label="Sampled X", color=:black, marker=:circle, markersize=6)
    plot!([X], [U], label="Sampled U", color=:black, marker=:circle, markersize=6)
    
    # Add vertical lines for current sample
    plot!([X, X], [0, pi_density(X)], color=:green, linestyle=:solid, label="", lw=3)
    plot!([X, X], [pi_density(X), c*q_density(X)], color=:red, linestyle=:solid, label="", lw=3)

    # Annotate the intermediate steps
    annotate!([(X+0.05, 0, text("Sampled X: $(round(X, digits=2))", :left, 10, :gray)),
               (X+0.05, U, text("Sampled U: $(round(U, digits=2))", :left, 10, color)),
               (X+0.05, pi_density(X), text("Status: $status", :left, 10, color))])
end

# Save as GIF
gif(anim, "rejection_sampling_with_steps.gif", fps=0.2)
