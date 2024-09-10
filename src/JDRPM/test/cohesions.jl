include("../utils.jl")

s1 = [1.0, 1.0, 11]
s2 = [1.0, 1.1, 1.1]
phi = 0.5
epsilon = 0.1
a = 5.0
mu_0 = [2.0, 1.0]
k0 = 0.5
v0 = 0.1
Psi = [1.0 0.5; 0.5 1.0]
lg=false

cohesion1(s1, s2, epsilon; lg=lg, M=1.0)
cohesion2(s1, s2, a; lg=lg, M=1.0)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=3, lg=lg, M=1.0)
cohesion3(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0)
cohesion3_4(s1, s2, mu_0, k0, v0, Psi; Cohesion=4, lg=lg, M=1.0)
cohesion4(s1, s2, mu_0, k0, v0, Psi; lg=lg, M=1.0)
cohesion5(s1, s2, phi; lg=lg, M=1.0)
cohesion6(s1, s2, phi; lg=lg, M=1.0)
