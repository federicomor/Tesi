include("../utils.jl")

s1 = [1.0, 1.0, 1.4]
s2 = [1.0, 1.1, 1.4]
# s1 = rand(10)
# s2 = rand(10)
phi = 0.5
epsilon = 0.1
a = 5.0
mu_0 = @SVector [2.,1.]
k0 = 0.5
v0 = 0.1
Psi = @SMatrix [1.0 0.5; 0.5 1.0]
lg=false
S = @MMatrix zeros(2,2)
M = 1.

cohesion1(s1, s2, epsilon, lg, 1.0)
cohesion2(s1, s2, a, lg, 1.0)
# cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=3, lg=lg, M=1.0)
cohesion3(s1, s2, mu_0, k0, v0, Psi, lg, 1.0, S)
# @btime cohesion3($s1, $s2, $mu_0, $k0, $v0, $Psi, $lg, $M, $S)
# cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=4, lg=lg, M=1.0)
cohesion4(s1, s2, mu_0, k0, v0, Psi, lg,1.0, S)
cohesion5(s1, s2, phi, lg, 1.0)
cohesion6(s1, s2, phi, lg, 1.0)

spatial_cohesion_idx = 3
sp_params_real = [mu_0, k0, v0, Psi]
@btime spatial_cohesion(spatial_cohesion_idx, s1, s2, sp_params_real, lg, M, S)
@btime spatial_cohesion($spatial_cohesion_idx, $s1, $s2, $sp_params_real, $lg, $M, $S)