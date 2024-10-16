include("../utils.jl")

s1 = [10,1,21,1,41,2.]
s2 = [1,11,1,31,1,1.5]
# s1 = rand(10)
# s2 = rand(10)
s1 = [1.0,2]
s2 = [1.1,2]

phi = 0.5
alphaS = 0.1
a = 5.0
mu_0 = @SVector [2.,1.]
k0 = 0.5
v0 = 0.1
Psi = @SMatrix [1.0 0.5; 0.5 1.0]
lg=true
S = @MMatrix zeros(2,2)
M = 1.

lC = @MVector zeros(2)
case = 1; add = false; 
lC[case]

cohesion1(s1, s2, alphaS, lg, M)
cohesion1!(s1, s2, alphaS, lg, M,case,add,lC)
cohesion2(s1, s2, a, lg, M)
cohesion2!(s1, s2, a, lg, M,case,add,lC)
# cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=3, lg=lg, M=M)
cohesion3(s1, s2, mu_0, k0, v0, Psi, lg, M, S)
cohesion3!(s1, s2, mu_0, k0, v0, Psi, lg, M, S,case,add,lC)
# @btime cohesion3($s1, $s2, $mu_0, $k0, $v0, $Psi, $lg, $M, $S)
# cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=4, lg=lg, M=M)
cohesion4(s1, s2, mu_0, k0, v0, Psi, lg,M, S)
cohesion4!(s1, s2, mu_0, k0, v0, Psi, lg,M, S,case,add,lC)
cohesion5(s1, s2, phi, lg, M)
cohesion5!(s1, s2, phi, lg, M,case,add,lC)
cohesion6(s1, s2, phi, lg, M)
cohesion6!(s1, s2, phi, lg, M,case,add,lC)

# spatial_cohesion_idx = 3
# sp_params_real = [mu_0, k0, v0, Psi]
# @btime spatial_cohesion(spatial_cohesion_idx, s1, s2, sp_params_real, lg, M, S)
# @btime spatial_cohesion($spatial_cohesion_idx, $s1, $s2, $sp_params_real, $lg, $M, $S)