@kwdef struct SpParams
	# cohesion 1
	alpha::Real
	# cohesion 2
	a::Real
	# cohesion 3 and 4
	mu_0::AbstractVector{Float64}
	k0::Real
	v0::Real
	Psi::AbstractMatrix{Float64}
	# cohesion 5 and 6
	phi::Real
end

#### the version without CvParams struct seems faster
@kwdef struct CvParams
	# cohesion 1
	alpha::Real
	# cohesion 2 and 3 (Gower)
	alpha_g::Real
	# cohesion 4
	mu_c::Real
	lambda_c::Real
	a_c::Real
	b_c::Real
end
