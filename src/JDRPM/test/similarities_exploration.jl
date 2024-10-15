include("../utils.jl")

using Plots
using Random
using Colors
using Printf

################################################
# COHESIONS
################################################
n = 105
m=[ 
	0.17257268  1.77336056
	-0.27246619 -0.23385607
	-0.07427892 -0.39462921
	-0.14965488  0.57109899
	0.23213730 -0.58727116
	-0.52771176  0.45863350
	-1.15200111  0.28535072
	-0.62015132  0.13869821
	-0.46710297  0.18576242
	-0.66498812  0.59908107
	-0.98541550  0.13653825
	-0.61767590  0.15529367
	-0.81595856  0.51925103
	-1.02951524  0.41608413
	-1.05492927  0.95921666
	-0.76893194  0.95473779
	-0.58948677  0.93972756
	-0.70081698  0.75935441
	0.16206402  1.77362600
	0.74776836  2.47332822
	-0.19027542  1.69082899
	-0.38802933  1.03629228
	-0.37572628  0.68305215
	-0.10047543  0.66700822
	-0.08007609  0.67778113
	-0.16033452  0.26860769
	-0.20197775  0.53435013
	-0.13746704  0.50346976
	-0.03424450 -0.74617019
	-0.28407023 -0.22636267
	-0.36592709 -0.16986809
	-0.03730457 -0.56661146
	-0.20180620  0.21963110
	0.33083299 -0.63024769
	-0.03058469 -0.08645782
	0.15424482 -0.27749227
	-0.65572152 -0.48469733
	-0.67655403 -0.50406485
	0.57086580  0.31669074
	0.55031615  0.57092311
	0.51816903  1.09293282
	0.76383709  0.56262317
	0.70233864  0.26019267
	1.24375986 -0.59935712
	1.21835336 -0.58431394
	0.53507530  0.25399864
	1.23476274 -0.56383755
	-0.83678149 -0.93745708
	-0.52944999  0.40997990
	0.36339885 -0.60566717
	-0.43893783  1.01776727
	-0.26984133  1.18153428
	-0.24846397  0.28516881
	-0.28523162  0.66559679
	-1.12826886  0.50042576
	-0.99853836 -0.77477504
	-0.95763756 -0.69795535
	1.10469849  0.02134920
	1.68146555 -0.82850639
	1.55993910 -0.89724216
	-0.57359272  0.17500947
	-0.38296006  1.06774994
	-1.13132421 -0.28222939
	-1.02645608 -0.23194322
	-0.89628911  1.40964846
	-0.92313498  1.75583675
	-0.88421079  1.44398886
	2.18236139 -1.30185042
	1.30836386  2.82124745
	0.11585567 -1.34503114
	1.51741619  1.13432204
	-1.29631106 -1.15311939
	1.36130215 -1.78087570
	1.29147508  1.13207135
	0.69704084 -1.41528973
	1.08309738 -1.65376288
	-0.90870474 -1.92768552
	0.70283658 -1.42754090
	1.33686731 -1.40178561
	1.42309632  1.57538727
	1.42987734  0.09677811
	-1.09450340  0.01422195
	-1.35021143  1.22216848
	-1.67709115  0.71843204
	-1.52454728 -0.17441819
	1.83268008 -0.51037165
	-1.48965693 -0.62061127
	-1.28738562  0.08522303
	1.36315679 -1.74709218
	-1.28332860 -1.12133133
	-1.53926846 -0.19646101
	-0.78582231 -1.48561793
	0.74243655 -1.10784808
	1.08213181 -1.23095143
	1.81915630 -1.54354878
	2.22485807 -1.34245561
	1.67800454 -1.10383629
	-0.04288964 -0.82296436
	-0.07089813 -0.84745973
	2.11312897 -0.69535035
	-0.42123552 -1.55190881
	-1.52385179  1.14907913
	-1.47435618  0.66320897
	1.03447808 -1.63184768
	1.50661989  0.06830806]
s1 = m[1:n,1]
s2 = m[1:n,2]

cls = vec([1 2 2 2 3 2 2 2 2 2 2 2 2 2 2 1 2 2 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 3 3 2 3 2 2 3 2 1 2 2 2 2 2 1 3 3 2 2 2 2 2 1 2 3 1 3 1 2 3 1 3 3 4 3 3 1 3 2 2 2 2 3 2 2 3 2 2 2 3 3 3 3 3 2 2 3 2 1 2 3 3])
cluster_labels = fill(0, n)

# Choose a color palette with distinct colors
# color_palette = [:blue, :green, :red]
color_palette = distinguishable_colors(6)[3:end] 
# color_palette = distinguishable_colors(maximum(cls), [RGB(0.9,1,1), RGB(0,0,0)])

n = 105
color_labels = [color_palette[cls[i]] for i in 1:n];

# Scatter plot with colored points based on clusters
p = scatter(s1, s2, color=color_labels,ms=5,legend=false, xlabel="", ylabel="")
alpha_sp = 0.1
lg = false
offsets_x = [0,-0.14,0,-0.1]
offsets_y = [0,-0.24,0,0.2]
for k in unique(cls)
	# k=1
    cohesion = cohesion1(s1[cls .== k], s2[cls .== k],alpha_sp,lg)
	str = @sprintf "%.2E" cohesion
	annotate!(mean(s1[cls .== k])+offsets_x[k], mean(s2[cls .== k])+offsets_y[k], 
	text("$str", 10, :black))
	# text("$(round(cohesion,digits=6))", 10, :black))
end
display(p)



s1_cl1 = s1[cls .== 1]; s2_cl1 = s2[cls .== 1]
s1_cl2 = s1[cls .== 2]; s2_cl2 = s2[cls .== 2]
s1_cl3 = s1[cls .== 3]; s2_cl3 = s2[cls .== 3]
s1_cl4 = s1[cls .== 4]; s2_cl4 = s2[cls .== 4]



alphas = range(0, 1.0, length=1000)
cohesion1_vals_cl1 = [cohesion1(s1_cl1, s2_cl1, alpha, lg) for alpha in alphas]
cohesion1_vals_cl2 = [cohesion1(s1_cl1, s2_cl2, alpha, lg) for alpha in alphas]
cohesion1_vals_cl3 = [cohesion1(s1_cl3, s2_cl3, alpha, lg) for alpha in alphas]
cohesion1_vals_cl4 = [cohesion1(s1_cl4, s2_cl4, alpha, lg) for alpha in alphas]
begin
plot(alphas,  cohesion1_vals_cl1,  label="cl1", lw=2, color=color_palette[1])
plot!(alphas, cohesion1_vals_cl2,  label="cl2", lw=2, color=color_palette[2])
plot!(alphas, cohesion1_vals_cl3,  label="cl3", lw=2, color=color_palette[3])
plot!(alphas, cohesion1_vals_cl4,  label="cl4", lw=2, color=color_palette[4])
xlabel!("alpha values")
ylabel!("")
title!("Cohesion 1")
end


as = range(0, 8, length=1000)
cohesion2_vals_cl1 = [cohesion2(s1_cl1, s2_cl1, a, lg) for a in as]
cohesion2_vals_cl2 = [cohesion2(s1_cl1, s2_cl2, a, lg) for a in as]
cohesion2_vals_cl3 = [cohesion2(s1_cl3, s2_cl3, a, lg) for a in as]
cohesion2_vals_cl4 = [cohesion2(s1_cl4, s2_cl4, a, lg) for a in as]
begin
plot(as,  cohesion2_vals_cl1,  label="cl1", lw=2, color=color_palette[1])
plot!(as, cohesion2_vals_cl2,  label="cl2", lw=2, color=color_palette[2])
plot!(as, cohesion2_vals_cl3,  label="cl3", lw=2, color=color_palette[3])
plot!(as, cohesion2_vals_cl4,  label="cl4", lw=2, color=color_palette[4])
xlabel!("a values")
ylabel!("")
title!("Cohesion 2")
end


mu0 = @SVector [0.,0.]; k0 = 1; v0 = 5; M=1.
Psi = @SMatrix [1.0 0.; 0. 1.0]; S = @MMatrix zeros(2,2)
v0s = range(0, 6, length=1000)
lg = true
cohesion3_vals_cl1 = [cohesion3(s1_cl1, s2_cl1, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion3_vals_cl2 = [cohesion3(s1_cl1, s2_cl2, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion3_vals_cl3 = [cohesion3(s1_cl3, s2_cl3, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion3_vals_cl4 = [cohesion3(s1_cl4, s2_cl4, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
begin
plot(v0s,  cohesion3_vals_cl1,  label="cl1", lw=2, color=color_palette[1])
plot!(v0s, cohesion3_vals_cl2,  label="cl2", lw=2, color=color_palette[2])
plot!(v0s, cohesion3_vals_cl3,  label="cl3", lw=2, color=color_palette[3])
plot!(v0s, cohesion3_vals_cl4,  label="cl4", lw=2, color=color_palette[4])
xlabel!("v0 values")
ylabel!("")
title!("Cohesion 3")
end

mu0 = @SVector [0.,0.]; k0 = 1; v0 = 5; M=1.
Psi = @SMatrix [1.0 0.; 0. 1.0]; S = @MMatrix zeros(2,2)
v0s = range(0, 6, length=1000)
lg = true
cohesion4_vals_cl1 = [cohesion4(s1_cl1, s2_cl1, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion4_vals_cl2 = [cohesion4(s1_cl1, s2_cl2, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion4_vals_cl3 = [cohesion4(s1_cl3, s2_cl3, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion4_vals_cl4 = [cohesion4(s1_cl4, s2_cl4, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
begin
plot(v0s,  cohesion4_vals_cl1,  label="cl1", lw=2, color=color_palette[1], legend=:topright)
plot!(v0s, cohesion4_vals_cl2,  label="cl2", lw=2, color=color_palette[2])
plot!(v0s, cohesion4_vals_cl3,  label="cl3", lw=2, color=color_palette[3])
plot!(v0s, cohesion4_vals_cl4,  label="cl4", lw=2, color=color_palette[4])
xlabel!("v0 values")
ylabel!("")
title!("Cohesion 4")
end



phis = range(0, 3.0, length=1000)
cohesion5_vals_cl1 = [cohesion5(s1_cl1, s2_cl1, phi, lg) for phi in phis]
cohesion5_vals_cl2 = [cohesion5(s1_cl1, s2_cl2, phi, lg) for phi in phis]
cohesion5_vals_cl3 = [cohesion5(s1_cl3, s2_cl3, phi, lg) for phi in phis]
cohesion5_vals_cl4 = [cohesion5(s1_cl4, s2_cl4, phi, lg) for phi in phis]
begin
plot(phis,  cohesion5_vals_cl1,  label="cl1", lw=2, color=color_palette[1])
plot!(phis, cohesion5_vals_cl2,  label="cl2", lw=2, color=color_palette[2])
plot!(phis, cohesion5_vals_cl3,  label="cl3", lw=2, color=color_palette[3])
plot!(phis, cohesion5_vals_cl4,  label="cl4", lw=2, color=color_palette[4])
xlabel!("phi values")
ylabel!("")
title!("Cohesion 5")
end


phis = range(0, 3.0, length=1000)
cohesion6_vals_cl1 = [cohesion6(s1_cl1, s2_cl1, phi, lg) for phi in phis]
cohesion6_vals_cl2 = [cohesion6(s1_cl1, s2_cl2, phi, lg) for phi in phis]
cohesion6_vals_cl3 = [cohesion6(s1_cl3, s2_cl3, phi, lg) for phi in phis]
cohesion6_vals_cl4 = [cohesion6(s1_cl4, s2_cl4, phi, lg) for phi in phis]
begin
plot(phis,  cohesion6_vals_cl1,  label="cl1",  lw=2, color=color_palette[1])
plot!(phis, cohesion6_vals_cl2,  label="cl2", lw=2, color=color_palette[2])
plot!(phis, cohesion6_vals_cl3,  label="cl3", lw=2, color=color_palette[3])
plot!(phis, cohesion6_vals_cl4,  label="cl4", lw=2, color=color_palette[4])
xlabel!("phi values")
ylabel!("")
title!("Cohesion 6")
end











################################################
# SIMILARITIES
################################################

#### numerical 1234
X_jt_numerical = [1,1,1,2,2,2,2,2,3,4,2,2,1,4,4]
X_jt_numerical = rand(15)
X_jt_numerical = [290, 76, 65, 207, 212, 48, 181, 124, 151, 129, 139, 152, 141, 116, 149, 114, 115, 222, 153, 100, 137, 118, 182, 425, 206, 383, 202, 279, 369, 307, 1225, 262, 229, 214, 249, 191, 190, 50, 80, 83, 74, 85, 97, 64, 59, 108, 57, 39, 79, 70, 77, 154, 269, 221, 47, 180, 349, 147, 18, 37, 22, 90, 162, 43, 247, 1194, 88, 133, 273, 215, 73, 109, 113, 14, 13, 11, 12, 16, 122, 272, 112, 280, 203, 305, 8, 640, 210, 200, 95, 60, 380, 55, 25, 1601, 62, 344, 130, 157, 400, 34, 91, 131, 620, 30, 9, 4, 61, 6, 765, 259, 301]
X_jt_numerical = [
	repeat([.1],5)...,
	repeat([-.1],2)...,
	repeat([.2],3)...,
	repeat([-.2],1)...,
	repeat([.3],2)...]
# altitudes of agrimonia dataset
# X_jt_numerical = [1,1,1,2,2,2,2,2,2,2,2,2,2,2,2]

lg = true
alphas = range(0, 3.0, length=1000)
similarity1_vals = [similarity1(X_jt_numerical, alpha, lg) for alpha in alphas]
similarity2_vals = [similarity2(X_jt_numerical, alpha, lg) for alpha in alphas]
similarity3_vals = [similarity3(X_jt_numerical, alpha, lg) for alpha in alphas]
# similarity4_vals = [similarity4(X_jt_numerical, 0,1,2,2, lg) for alpha in alphas]

begin
plot(alphas, similarity1_vals,  label="similarity1", lw=2, color=:blue,legend=:right)
plot!(alphas, similarity2_vals, label="similarity2", lw=2, color=:green)
plot!(alphas, similarity3_vals, label="similarity3", lw=2, color=:red)
# plot!(alphas, similarity4_vals, label="similarity4", lw=2)
xlabel!("alpha")
ylabel!("Similarity Output")
title!("Effect of alpha on similarities (Numerical Data)")
end

#### categorical 123