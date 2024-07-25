using Distributions
using LinearAlgebra
using Random

n = 10
T = 10
Si_iter = rand(collect(1:5),n,T)
gamma_iter = rand((0,1),n,T)