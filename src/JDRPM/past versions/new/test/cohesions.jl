include("../utils.jl")

n = 1000
s1 = rand(n)
s2 = rand(n)
phi = 0.5
epsilon = 0.1
a = 5.0
mu_0 = @SVector [1.,2.]
k0 = 0.5
v0 = 0.1
Psi = @SMatrix [1.0 0.5; 0.5 1.0]
lg=false

function cohesion4_propose(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::AbstractMatrix{Float64}; lg::Bool, M::Real=1.0, S=@MMatrix zeros(2, 2))::Float64
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

@timev cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=3, lg=lg, M=1.0)
@timev cohesion3(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0)
@timev cohesion3_propose(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0)

@btime cohesion3_4($s1, $s2, $mu_0, $k0, $v0, $Psi; Cohesion=3, lg=lg, M=1.0)
@btime cohesion3($s1, $s2, $mu_0, $k0, $v0, $Psi, lg=lg, M=1.0)
S = @MMatrix zeros(2, 2)
@btime cohesion3_propose($s1, $s2, $mu_0, $k0, $v0, $Psi, lg=false, S=$S)

@timev cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=4, lg=lg, M=1.0)
@timev cohesion4(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0)
@timev cohesion4_propose(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0,S = S)
@timev cohesion4_v3(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0,S = S)

@btime cohesion3_4($s1, $s2, $mu_0, $k0, $v0, $Psi, Cohesion=4, lg=lg, M=1.0)
@btime cohesion4($s1, $s2, $mu_0, $k0, $v0, $Psi; lg=lg, M=1.0)
@btime cohesion4_v2($s1, $s2, $mu_0, $k0, $v0, $Psi; lg=$lg, M=1.0)
S=@MMatrix zeros(2, 2)
@btime cohesion4_v3($s1, $s2, $mu_0, $k0, $v0, $Psi; lg=lg, M=1.0, S=$S)
@btime cohesion4_propose($s1, $s2, $mu_0, $k0, $v0, $Psi; lg=lg, M=1.0, S=$S)

b1 = @benchmark cohesion4($s1, $s2, $mu_0, $k0, $v0, $Psi; lg=lg, M=1.0)
S=@MMatrix zeros(2, 2)
b2 = @benchmark cohesion4_propose($s1, $s2, $mu_0, $k0, $v0, $Psi; lg=lg, M=1.0, S=$S)




cohesion1(s1, s2, epsilon; lg=lg, M=1.0)
cohesion2(s1, s2, a; lg=lg, M=1.0)
@btime cohesion3_4($s1, $s2, $mu_0, $k0, $v0, $Psi; Cohesion=3, lg=lg, M=1.0)
@btime cohesion3($s1, $s2, $mu_0, $k0, $v0, $Psi; lg=lg, M=1.0)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=4, lg=lg, M=1.0)
cohesion4(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0)
cohesion5(s1, s2, phi; lg=lg, M=1.0)
cohesion6(s1, s2, phi; lg=lg, M=1.0)


a = 20
b = a-0.5
@btime log(gamma($a)*gamma($b))
@btime llgamma($a) + llgamma($b)
