using MCMCChains
using DataFrames
using CSV
# chn = Chains(rand(100, 2, 1), ["a.1", "b.1"]);

chn = Chains(
	hcat(out[end-4],out[end-6],out[end-7]',out[end-8]',out[6]',out[3]'),
	["lambda2","phi0",[string("tau2_t", i) for i in 1:T]...,
	[string("theta_t", i) for i in 1:T]...,
	[string("eta1_j", i) for i in 1:N]...,
	[string("alpha_t", i) for i in 1:T]...,
	]
)

using PrettyTables
ss = DataFrame(summarystats(chn))
ss = DataFrame(ess_rhat(chn))
ss[!,:mcse] = DataFrame(mcse(chn))[!,2]
@show ss;
@show ss[!,[1,4,5,6,7]];
# CSV.write(log_file,ss)