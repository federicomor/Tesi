include("../utils.jl")

using Plots
using Random
using Colors
using Printf
using LaTeXStrings
using Distributions

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
color_palette = palette(:tab10)[[1,2,3,5]]
n = 105
color_labels = [color_palette[cls[i]] for i in 1:n];

save_plot = false
p = scatter(s1, s2, color=color_labels,ms=5,legend=false, xlabel="", ylabel="",
# aspect_ratio = 0.7,
markerstrokewidth=0,
# dpi=300, size=(500,400),
# title="Clusters test case")
title="Cohesions analysis test case")
# title="Cluster subdivision at test")
alpha_sp = 0.1
lg = false
offsets_x = [0,-0.14,0,-0.1]
offsets_y = [0,-0.24,0,0.2]
for k in unique(cls)
	# k=1
    cohesion = cohesion1(s1[cls .== k], s2[cls .== k],alpha_sp,lg)
	str = @sprintf "%.2E" cohesion
	annotate!(mean(s1[cls .== k])+offsets_x[k], mean(s2[cls .== k])+offsets_y[k], 
	# text("$str", 10, :black))
	# text("$str", 10, color_palette[k]))
	text("cluster $k", 10, color_palette[k]))
	# text("$(round(cohesion,digits=6))", 10, :black))
end
display(p) savefig(p,"cohesions_test.pdf")


s1_cl1 = s1[cls .== 1]; s2_cl1 = s2[cls .== 1]
s1_cl2 = s1[cls .== 2]; s2_cl2 = s2[cls .== 2]
s1_cl3 = s1[cls .== 3]; s2_cl3 = s2[cls .== 3]
s1_cl4 = s1[cls .== 4]; s2_cl4 = s2[cls .== 4]


lg = true
alphas = range(0, 1.0, length=1000)
cohesion1_vals_cl1 = [cohesion1(s1_cl1, s2_cl1, alpha, lg) for alpha in alphas]
cohesion1_vals_cl2 = [cohesion1(s1_cl1, s2_cl2, alpha, lg) for alpha in alphas]
cohesion1_vals_cl3 = [cohesion1(s1_cl3, s2_cl3, alpha, lg) for alpha in alphas]
cohesion1_vals_cl4 = [cohesion1(s1_cl4, s2_cl4, alpha, lg) for alpha in alphas]
begin
plot(alphas,  cohesion1_vals_cl1,  label="cluster 1", lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(alphas, cohesion1_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(alphas, cohesion1_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(alphas, cohesion1_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
# xlabel!("alpha values")
xlabel!(L"$\alpha$ values")
ylabel!("")
title!("Cohesion 1 (log=$lg)") savefig("cohesion1.pdf")
end

# scatter(s1_cl1,s2_cl1,color=:red)
# scatter!(s1_cl2,s2_cl2,color=:green)
# scatter!(s1_cl3,s2_cl3,color=:orange)


lg = false
as = range(0, 6, length=1000)
cohesion2_vals_cl1 = [cohesion2(s1_cl1, s2_cl1, a, lg) for a in as]
cohesion2_vals_cl2 = [cohesion2(s1_cl1, s2_cl2, a, lg) for a in as]
cohesion2_vals_cl3 = [cohesion2(s1_cl3, s2_cl3, a, lg) for a in as]
cohesion2_vals_cl4 = [cohesion2(s1_cl4, s2_cl4, a, lg) for a in as]
begin
plot(as,  cohesion2_vals_cl1,  label="cluster 1", lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(as, cohesion2_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(as, cohesion2_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(as, cohesion2_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
plot!(legend=:bottomright)
xlabel!(L"$a$ values")
ylabel!("")
title!("Cohesion 2 (log=$lg)")
# savefig("cohesion2.pdf")
end


begin
mu0 = @SVector [0.,0.]; k0 = 1; v0 = 5; M=1.
Psi = @SMatrix [1.0 0.; 0. 1.0]; S = @MMatrix zeros(2,2)
v0s = range(0, 6, length=1000)
lg = true
cohesion3_vals_cl1 = [cohesion3(s1_cl1, s2_cl1, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion3_vals_cl2 = [cohesion3(s1_cl1, s2_cl2, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion3_vals_cl3 = [cohesion3(s1_cl3, s2_cl3, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion3_vals_cl4 = [cohesion3(s1_cl4, s2_cl4, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
plot(v0s,  cohesion3_vals_cl1,  label="cluster 1", lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(v0s, cohesion3_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(v0s, cohesion3_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(v0s, cohesion3_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
plot!(legend=:bottomright)
xlabel!(L"$\nu_0$ values")
ylabel!("")
# title!("Cohesion 3 (log=$lg)")
savefig("cohesion3.pdf")
end


begin
mu0 = @SVector [0.,0.]; k0 = 1; v0 = 5; M=1.
Psi = @SMatrix [1.0 0.; 0. 1.0]; S = @MMatrix zeros(2,2)
v0s = range(0, 6, length=1000)
lg = true
cohesion4_vals_cl1 = [cohesion4(s1_cl1, s2_cl1, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion4_vals_cl2 = [cohesion4(s1_cl1, s2_cl2, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion4_vals_cl3 = [cohesion4(s1_cl3, s2_cl3, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
cohesion4_vals_cl4 = [cohesion4(s1_cl4, s2_cl4, mu0, k0, v0, Psi, lg, M, S) for v0 in v0s]
plot(v0s,  cohesion4_vals_cl1,  label="cluster 1", lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(v0s, cohesion4_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(v0s, cohesion4_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(v0s, cohesion4_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
# xlabel!("ν₀ values")
plot!(legend=:bottomright)
xlabel!(L"$\nu_0$ values")
ylabel!("")
title!("Cohesion 4 (log=$lg)")
# savefig("cohesion4.pdf")
end

begin
lg = true
phis = range(0, 3.0, length=1000)
cohesion5_vals_cl1 = [cohesion5(s1_cl1, s2_cl1, phi, lg) for phi in phis]
cohesion5_vals_cl2 = [cohesion5(s1_cl1, s2_cl2, phi, lg) for phi in phis]
cohesion5_vals_cl3 = [cohesion5(s1_cl3, s2_cl3, phi, lg) for phi in phis]
cohesion5_vals_cl4 = [cohesion5(s1_cl4, s2_cl4, phi, lg) for phi in phis]
plot(phis,  cohesion5_vals_cl1,  label="cluster 1", lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(phis, cohesion5_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(phis, cohesion5_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(phis, cohesion5_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
# xlabel!("ϕ values")
xlabel!(L"$\varphi$ values")
ylabel!("")
title!("Cohesion 5 (log=$lg)")
# title!("Cohesion 5 (log=$lg)" * L" wrt parameter $\phi$")
# savefig("cohesion5.pdf")
end


begin
lg = true
phis = range(0, 3.0, length=1000)
cohesion6_vals_cl1 = [cohesion6(s1_cl1, s2_cl1, phi, lg) for phi in phis]
cohesion6_vals_cl2 = [cohesion6(s1_cl1, s2_cl2, phi, lg) for phi in phis]
cohesion6_vals_cl3 = [cohesion6(s1_cl3, s2_cl3, phi, lg) for phi in phis]
cohesion6_vals_cl4 = [cohesion6(s1_cl4, s2_cl4, phi, lg) for phi in phis]
plot(phis,  cohesion6_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(phis, cohesion6_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(phis, cohesion6_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(phis, cohesion6_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
# xlabel!("phi values")
xlabel!(L"$\varphi$ values")
ylabel!("")
title!("Cohesion 6 (log=$lg)")
# savefig("cohesion6.pdf")
end











################################################
# SIMILARITIES
################################################
#### numerical - sims: 1234

# altitude
alts = vec([0.420135834,-0.472240779,-0.518110605,0.074027148,-0.589000336,-0.034392441,-0.201191808,-0.305441412,-0.313781381,0.136576910,-0.217871745,-0.297101444,0.094877069,0.069857163,0.807944363,0.053177227,0.374266008,0.749564584,0.491025565,4.319071039,0.303376277,0.103217037,0.374266008,0.249166483,0.249166483,0.007307401,0.003137417,-0.030222457,-0.580660368,-0.455560843,-0.434710922,-0.543130510,-0.338801286,-0.626530194,-0.459730827,-0.497260685,-0.468070795,-0.468070795,-0.167831935,0.332566167,0.132406926,0.666164901,-0.176171903,-0.714099862,-0.714099862,-0.497260685,-0.697419925,-0.413861001,-0.113622140,-0.609850257,0.240826515,4.189801530,-0.234551681,0.349246103,0.107387021,-0.484750732,-0.434710922,-0.317951365,-0.730779798,-0.743289751,-0.280421507,0.345076119,-0.322121349,-0.334631302,0.378435993,0.057347211,0.482685597,-0.755799703,1.879630296,0.086537100,0.044837258,-0.393011080,-0.626530194,-0.484750732,-0.538960526,-0.580660368,0.795434410,-0.559810447,-0.684909972,5.886985090,-0.530620558,-0.272081539,0.057347211,0.645314980,-0.247061634,-0.684909972,-0.305441412,-0.134472061,-0.647380115,-0.409691017,-0.242891650,1.796230613,-0.664060052,-0.697419925,-0.722439830,-0.751629719,-0.772479640,-0.534790542,-0.534790542,-0.764139672,2.400878318,0.290866325,0.466005660,-0.543130510,-0.589000336])

# wspeed10mmax
alts  = vec([-1.5831941135,0.6462411450,0.7692486349,0.0050318893,0.9707715439,-0.0289914590,-0.1729363940,0.0086556778,-0.0088593003,-0.0038262606,-0.0165095206,0.0020120655,-0.0716716355,-0.1815932222,-0.7527425656,-0.6007447671,-0.4934403610,-0.1511936625,-1.5781409416,-1.5945486510,-1.5277098840,-0.6389958687,0.0303984093,0.0040252813,-0.0265755999,-0.0285888158,0.0430816693,-0.0207372739,0.8403151553,0.6450332155,0.6355711009,0.8731305740,-0.0006051152,1.0140556852,0.4024407027,0.6653666958,0.6500662552,0.6703997355,-0.5012919029,-0.6832866179,-0.9705725231,-0.7485148122,-0.3925782458,0.9713755086,0.9206424686,-0.3366108445,0.8483680187,0.0255666912,-0.0469090803,1.0535147163,-0.6279231814,-0.8753474123,-0.0249650272,0.0330155899,-0.2540689937,0.4108962094,0.5135702190,0.1089138284,1.2808067884,1.1251852014,-0.0118791241,-0.8109245043,0.4104935662,0.3911666938,-1.4589988262,-1.6156471533,-1.4815065797,1.1004226461,-1.3345820852,0.0547583214,-1.8693727499,-0.1622663498,0.5665177964,-1.6715138938,0.6017490742,0.5208177961,0.6448318939,0.5973199993,0.5739666952,-1.9038591378,-0.0505328688,0.1238116258,-1.3019478559,-0.7946174557,0.0889829912,1.1181389458,0.5330984129,-0.0052355117,0.5876565631,-0.1634742793,0.0928081014,-0.1956857333,0.7946151549,0.6854988546,1.0796865226,1.0786799147,1.1755155982,0.7128785905,0.7064362997,1.1058583290,-0.3678156905,-1.2541138467,-0.7638152529,0.5512173558,0.0968345331])

# blh layer max
alts = vec([-0.502001858,0.088231371,0.156189698,0.115436612,0.204830478,0.028270289,-0.096983372,-0.008100745,0.076545898,-0.050716204,-0.066929797,-0.007772091,-0.077337171,-0.091688393,-0.163955738,-0.145551118,-0.046735840,-0.069157340,-0.502965909,-0.097786749,-0.373783008,-0.040418381,0.057301385,0.081074019,0.076253761,0.182335943,0.159403203,0.155824527,0.198403468,0.085894276,0.066759314,0.160498716,0.153377881,0.216515951,0.183102802,0.207642295,-0.001892838,-0.003682176,0.138405869,-0.009597947,-0.358445825,-0.062620779,0.084725729,0.030534349,0.031410759,0.149507068,0.029657938,-0.001965872,0.024874198,0.202091696,-0.032128999,-0.119623976,0.132051893,0.077714445,-0.104688481,-0.016171025,-0.029390216,-0.007333886,-0.014564272,-0.003499590,-0.001454633,-0.104213759,-0.064519668,-0.049182485,-0.238012422,-0.271535121,-0.242248405,0.047953007,-0.150283735,0.351410127,-0.390836495,-0.020772180,0.249125723,-0.424651332,0.285606309,0.272058464,0.386210926,0.294041759,0.183687076,-0.485233204,-0.107609849,-0.068682618,-0.251706335,-0.297608333,-0.102680040,-0.081719224,-0.092637837,-0.101146322,0.244122880,-0.022488484,-0.103592968,0.112807380,0.203041140,0.133695163,0.112113556,0.081585258,0.009390947,0.201945627,0.198366951,-0.074963560,0.324460506,-0.333431610,-0.219461734,0.273519148,-0.111845833])

# nh3sum
alts = vec([-0.42428648,0.52533376,0.43142692,0.28885668,0.21448279,0.30619978,-0.21348475,0.89531076,0.29568308,0.08895826,-0.02276479,0.88246903,0.07803603,-0.04186517,-0.21105158,-0.21553943,-0.25828212,-0.05048941,-0.42104496,-0.44637291,-0.37200577,-0.31247422,-0.16936327,0.13956820,0.11353328,0.64509979,0.12452309,0.43477929,0.02119448,0.52979457,0.42151851,0.33434344,0.64537014,0.21419892,0.77067839,0.84461972,-0.19142401,-0.10487075,0.27984044,-0.25698443,-0.33696813,-0.33126369,0.43342753,0.31670296,0.50694981,0.53303880,0.56777906,-0.38509352,0.38800836,0.40324270,-0.30562079,-0.38734826,0.39757883,-0.14918148,-0.17878505,-0.34105044,-0.33868486,0.28808618,-0.13051366,0.01647683,0.81988249,-0.34433522,-0.03825597,-0.17079614,-0.29991636,-0.30858114,-0.29169765,-0.20726665,-0.45565004,-0.32764098,-0.39937353,-0.15292586,0.24541108,-0.42481772,0.03219782,0.17702549,-0.41555681,0.03111641,-0.07395598,-0.42546657,0.16834719,-0.24241244,-0.40483870,-0.43834750,-0.28591211,0.07598136,-0.31132522,-0.14588318,0.16971247,-0.14254434,-0.30751325,-0.45118653,-0.07363156,0.20857560,-0.32884404,-0.17625726,-0.10638473,-0.09849044,-0.09950426,-0.03007782,-0.45641244,-0.42511106,-0.34862031,0.19877533,0.30221208])

# cls = vec([1 2 2 1 2 2 2 2 2 2 2 2 1 2 1 1 1 1 1 3 1 1 1 1 1 2 1 2 2 2 2 2 2 4 2 2 2 2 1 1 1 1 2 4 4 1 4 2 2 2 1 3 2 1 2 2 2 2 4 4 2 1 2 2 1 1 1 4 3 5 1 2 4 1 4 4 3 4 4 3 4 2 1 1 2 4 2 2 4 2 2 3 4 4 4 4 4 4 4 4 3 1 1 4 4])
cls = vec([1 2 2 1 2 2 2 2 2 2 2 2 1 2 1 1 1 1 1 3 1 1 1 1 1 2 1 2 2 2 2 2 2 4 2 2 2 2 1 1 1 1 2 4 4 1 4 2 2 2 1 3 2 1 2 2 2 2 4 4 2 1 2 2 1 1 1 4 3 2 1 2 4 1 4 4 3 4 4 3 4 2 1 1 2 4 2 2 4 2 2 3 4 4 4 4 4 4 4 4 3 1 1 4 4])
cls[80] = 5
maximum(cls)
n = 105
cluster_labels = fill(0, n)
color_palette = palette(:tab10)[[1,2,3,4,5]]
# color_palette = palette(:tab10)[[1,2,3,5]]
n = 105
color_labels = [color_palette[cls[i]] for i in 1:n];

scatter(1:105,alts,markerstrokewidth=0,color=color_labels,ms=4,legend=false, xlabel="", ylabel="Altitude",
# aspect_ratio = 0.7,
# title="Clusters test case")
title="Similarities analysis test case")
# savefig("similarity_test.pdf")


X_cl1 = alts[cls .== 1]
X_cl2 = alts[cls .== 2]
X_cl3 = alts[cls .== 3]
X_cl4 = alts[cls .== 4]
X_cl5 = alts[cls .== 5]

begin
p = scatter(1:length(X_cl1),X_cl1,markerstrokewidth=0,color=color_palette[1],ms=4,legend=false, xlabel="", ylabel="Altitude",title="Similarities analysis test case")
scatter!(length(X_cl1)+1:length(X_cl1)+length(X_cl2),X_cl2,
markerstrokewidth=0,color=color_palette[2],ms=4,legend=false)
scatter!(length(X_cl1)+length(X_cl2)+1:length(X_cl1)+length(X_cl2)+length(X_cl3),X_cl3,
markerstrokewidth=0,color=color_palette[3],ms=4,legend=false)
scatter!(length(X_cl1)+length(X_cl2)+length(X_cl3)+1:length(X_cl1)+length(X_cl2)+length(X_cl3)+length(X_cl4),X_cl4,
markerstrokewidth=0,color=color_palette[4],ms=4,legend=false)
scatter!(length(X_cl1)+length(X_cl2)+length(X_cl3)+length(X_cl4)+1:length(X_cl1)+length(X_cl2)+length(X_cl3)+length(X_cl4)+length(X_cl5),X_cl5,
markerstrokewidth=0,color=color_palette[5],ms=4,legend=false)
offsets_x = [18,55,75,92,103]
offsets_y = [1.08,0.5,3,-0.1,5.5]
for k in 1:5
	# k=1
	annotate!(offsets_x[k],offsets_y[k], 
	text("cluster $k", 10, color_palette[k]))
end
display(p)
end 
# savefig(p,"similarity_sorted_test.pdf")



begin
# lg = false
lg = true
cv_weight = 5
alphas = range(0, 2, length=1000)
similarity1_vals_cl1 = [similarity1(X_cl1, alpha, lg, cv_weight) for alpha in alphas]
similarity1_vals_cl2 = [similarity1(X_cl2, alpha, lg, cv_weight) for alpha in alphas]
similarity1_vals_cl3 = [similarity1(X_cl3, alpha, lg, cv_weight) for alpha in alphas]
similarity1_vals_cl4 = [similarity1(X_cl4, alpha, lg, cv_weight) for alpha in alphas]
similarity1_vals_cl5 = [similarity1(X_cl5, alpha, lg, cv_weight) for alpha in alphas]
plot(alphas,  similarity1_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(alphas, similarity1_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(alphas, similarity1_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(alphas, similarity1_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
plot!(alphas, similarity1_vals_cl5,  label="cluster 4", lw=2, color=color_palette[5])
xlabel!(L"$\alpha$ values")
ylabel!("")
title!("Similarity 1 (log=$lg)")
# savefig("similarity1_log_$lg.pdf")
end

lg = true
myX1_cat = ["a","a","b"]
myX2_cat = ["a","a","a"]
R=0
similarity2(myX1_cat,2,R,lg)
similarity2(myX2_cat,2,R,lg)

similarity3(myX1_cat,2,R,lg)
similarity3(myX2_cat,2,R,lg)

similarity1(myX1_cat,2,lg)
similarity1(myX2_cat,2,lg)


lg = true
myX1 = [0.0, 1.0,50,0, 20]
myX2 = [0.0, 0.0,0,0, 2.0]
R = max(maximum(myX2),maximum(myX1)) - min(minimum(myX2),minimum(myX1))
similarity2(myX1,2,R,lg)
similarity2(myX2,2,R,lg)

similarity3(myX1,2,R,lg)
similarity3(myX2,2,R,lg)

similarity1(myX1,1,true)
similarity1(myX2,1,true)
similarity4(myX1,1,5,2,2,true)
similarity4(myX2,1,5,2,2,true)






begin
lg = true
alphas = range(0, 3.0, length=1000)
R = maximum(alts)-minimum(alts)
similarity2_vals_cl1 = [similarity2(X_cl1, alpha, R, lg) for alpha in alphas]
similarity2_vals_cl2 = [similarity2(X_cl2, alpha, R, lg) for alpha in alphas]
similarity2_vals_cl3 = [similarity2(X_cl3, alpha, R, lg) for alpha in alphas]
similarity2_vals_cl4 = [similarity2(X_cl4, alpha, R, lg) for alpha in alphas]
similarity2_vals_cl5 = [similarity2(X_cl5, alpha, R, lg) for alpha in alphas]
plot(alphas,  similarity2_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(alphas, similarity2_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(alphas, similarity2_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(alphas, similarity2_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
plot!(alphas, similarity2_vals_cl5,  label="cluster 4", lw=2, color=color_palette[5])
xlabel!(L"$\alpha$ values")
ylabel!("")
title!("Similarity 2 (log=$lg)")
# savefig("similarity2.pdf")
end


begin
lg = true
alphas = range(0, 2.0, length=1000)
R = maximum(alts)-minimum(alts)
similarity3_vals_cl1 = [similarity3(X_cl1, alpha, R, lg) for alpha in alphas]
similarity3_vals_cl2 = [similarity3(X_cl2, alpha, R, lg) for alpha in alphas]
similarity3_vals_cl3 = [similarity3(X_cl3, alpha, R, lg) for alpha in alphas]
similarity3_vals_cl4 = [similarity3(X_cl4, alpha, R, lg) for alpha in alphas]
similarity3_vals_cl5 = [similarity3(X_cl5, alpha, R, lg) for alpha in alphas]
plot(alphas, similarity3_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(alphas, similarity3_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(alphas, similarity3_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(alphas, similarity3_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
plot!(alphas, similarity3_vals_cl5,  label="cluster 4", lw=2, color=color_palette[5])
xlabel!(L"$\alpha$ values")
ylabel!("")
title!("Similarity 3 (log=$lg)")
# savefig("similarity3.pdf")
end


begin
lg = true
mu = 0; lambda = 1; 

a_c = 1; b_c = 1.6
# a_c = 1; b_c = 1.5
# a_c = 1; b_c = 1
# a_c = 1.5; b_c = 2
# a_c = 1; b_c = .5
# a_c = 2; b_c = 2
# a_c = 2; b_c = 1.5

cv_weight = 1
println("$a_c $b_c")
println(quantile.(InverseGamma(a_c,b_c),[0.1 0.9]))
print(quantile.(Normal(0,sqrt(quantile(InverseGamma(a_c,b_c),0.9))),[0.1 0.9]))
lambdas = range(0, 2, length=1000)
# acs = range(0.1, 3.0, length=1000)
similarity4_vals_cl1 = [similarity4(X_cl1, mu, lambda, a_c, b_c, lg, cv_weight)
for lambda in lambdas]
# for b_c in acs]
# for b_c in acs]
similarity4_vals_cl2 = [similarity4(X_cl2, mu, lambda, a_c, b_c, lg, cv_weight)
for lambda in lambdas]
# for b_c in acs]
# for b_c in acs]
similarity4_vals_cl3 = [similarity4(X_cl3, mu, lambda, a_c, b_c, lg, cv_weight)
for lambda in lambdas]
# for b_c in acs]
# for b_c in acs]
similarity4_vals_cl4 = [similarity4(X_cl4, mu, lambda, a_c, b_c, lg, cv_weight)
for lambda in lambdas]
# for b_c in acs]
# for b_c in acs]
similarity4_vals_cl5 = [similarity4(X_cl5, mu, lambda, a_c, b_c, lg, cv_weight)
for lambda in lambdas]
# for b_c in acs]
# for b_c in acs]
plot(lambdas,  similarity4_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(lambdas, similarity4_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(lambdas, similarity4_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(lambdas, similarity4_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
plot!(lambdas, similarity4_vals_cl5,  label="cluster 4", lw=2, color=color_palette[5])
xlabel!(L"$\lambda_0$ values")
ylabel!("")
# title!("Similarity 4 (log=$lg)")
# annotate!(:bottomright,
# text("σ² ∼ invGamma($a_c, $b_c)",8,:gray,:right))
title!("Similarity 4 (log=$lg, σ²∼IG($a_c, $b_c))")
# savefig("similarity4_ac$(a_c)_bc$(b_c).pdf")
end



a=4.5
b=1.4
# nice for bhl max
a = 4
b = 2

a = 7.5
b = 2
quantile.(InverseGamma(a,b),[0.05 0.1 0.9 0.95 0.99])
print(quantile.(Normal(0,sqrt(quantile(InverseGamma(a,b),0.9))),[0.1 0.9]))
mean(InverseGamma(a,b))
print(quantile.(Normal(0,mean(InverseGamma(a,b))),[0.1 0.9]))
# ar r bl
begin
lg = true
mu = 0; lambda = 1; 
b_c = 2
cv_weight = 1
acs = range(0.01, 10.0, length=1000)
similarity4_vals_cl1 = [similarity4(X_cl1, mu, lambda, a_c, b_c, lg, cv_weight)
# for lambda in lambdas]
for a_c in acs]
# for b_c in acs]
similarity4_vals_cl2 = [similarity4(X_cl2, mu, lambda, a_c, b_c, lg, cv_weight)
# for lambda in lambdas]
for a_c in acs]
# for b_c in acs]
similarity4_vals_cl3 = [similarity4(X_cl3, mu, lambda, a_c, b_c, lg, cv_weight)
# for lambda in lambdas]
for a_c in acs]
# for b_c in acs]
similarity4_vals_cl4 = [similarity4(X_cl4, mu, lambda, a_c, b_c, lg, cv_weight)
# for lambda in lambdas]
for a_c in acs]
# for b_c in acs]
similarity4_vals_cl5 = [similarity4(X_cl5, mu, lambda, a_c, b_c, lg, cv_weight)
# for lambda in lambdas]
for a_c in acs]
# for b_c in acs]
plot(acs,  similarity4_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(acs, similarity4_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(acs, similarity4_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(acs, similarity4_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
plot!(acs, similarity4_vals_cl5,  label="cluster 4", lw=2, color=color_palette[5])
xlabel!(L"$a_0$ values")
ylabel!("")
# title!("Similarity 4 (log=$lg)")
# annotate!(:bottomright,
# text("σ² ∼ invGamma($a_c, $b_c)",8,:gray,:right))
title!("Similarity 4 (log=$lg, b₀=$b_c)")
# savefig("similarity4_lambda1_bc$(b_c).pdf")
end
	















#### categorical - sims: 123
#### precipitation type
Xj_cat_nums = vec([5 1 1 1 1 1 0 1 1 0 0 1 0 0 1 0 1 0 5 5 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0 0 0 5 0 1 5 1 1 1 1 1 1 1 0 1 1 0 1 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1])
Xj_cat = vec(string.([5 1 1 1 1 1 0 1 1 0 0 1 0 0 1 0 1 0 5 5 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0 0 0 5 0 1 5 1 1 1 1 1 1 1 0 1 1 0 1 0 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1]))


#### wind direction
Xj_cat = ["S"  "E"  "E"  "E"  "W"  "E"  "SE" "E"  "E"  "S"  "SE" "E"  "S"  "SE" "S"  "S"  "S"  "S"  "S"  "N"  "S"  "S"  "E"  "E"  "E"  "E"  "E"  "E"  "E"  "E"  "E"  "W"  "E"  "W"  "E"  "E"  "E"  "E"  "E"  "SW" "SW" "SW" "E"  "E"  "E"  "E"  "E"  "SE" "E"  "W"  "S"  "S"  "E"  "E"  "SE" "SE" "SE" "E"  "W"  "W"  "E"  "S"  "SE" "SE" "S"  "S"  "S"  "W"  "NW" "W"  "N"  "S"  "E"  "W"  "W"  "E"  "N"  "W"  "W"  "W"  "SE" "SE" "S"  "SW" "SE" "SE" "S"  "SE" "E"  "S"  "SE" "N"  "W"  "W"  "W"  "W"  "W"  "E"  "E"  "SE" "N"  "SW" "SE" "E"  "SE"]
cls = vec([1 2 2 1 2 2 2 2 2 2 2 2 1 2 1 1 1 1 1 3 1 1 1 1 1 2 1 2 2 2 2 2 2 4 2 2 2 2 1 1 1 1 2 4 4 1 4 2 2 2 1 3 2 1 2 2 2 2 4 4 2 1 2 2 1 1 1 4 3 2 1 2 4 1 4 4 3 4 4 3 4 2 1 1 2 4 2 2 4 2 2 3 4 4 4 4 4 4 4 4 3 1 1 4 4])
maximum(cls)
n = 105
cluster_labels = fill(0, n)
# color_palette = palette(:tab10)[[1,2,3,4,5]]
color_palette = palette(:tab10)[[1,2,3,5]]
n = 105
color_labels = [color_palette[cls[i]] for i in 1:n];
# Mapping categories to numerical values for visualization
category_dict = Dict("N" => 1, "S" => 2, "E" => 3, "W" => 4, "SE" => 5, "SW" => 6, "NW" => 7)
Xj_num = [category_dict[x] for x in Xj_cat]
# Define marker shape based on cluster
markers = [:circle, :square, :triangle, :star5]
marker_shape = [markers[c] for c in cls]
# Plot
scatter(1:length(Xj_cat), vec(Xj_num), marker=:o, c = color_labels, legend = false,
        xlabel = "Index", ylabel = "Category", title = "Categorical Data by Cluster",
		markerstrokewidth=0)
# Add category labels on the y-axis
yticks!(1:7, ["N", "S", "E", "W", "SE", "SW", "NW"])



using Plots.PlotMeasures
using Statistics
n = 105
s1 = m[1:n,1]
s2 = m[1:n,2]
s1m = mean(s1)
s2m = mean(s2)
R = zeros(105)

for i in 1:105
	R[i] = log((s1[i]-s1m)^2+(s2[i]-s2m)^2)
end
# Rscaled =  (R .- mean(R)) ./ sqrt(var(R))
# R = Rscaled
labels = String["E", "NE", "N", "NW", "W", "SW", "S", "SE"]
angles_degrees = [0, 45, 90, 135, 180, 225, 270, 315]
angles_radians = deg2rad.(angles_degrees)

θ = LinRange(0, 2pi, 106)
scatter(θ[1:105], R, proj=:polar, ms=4, m=:o,
markercolor=color_labels, fa=1,
lims=(0,maximum(R)),
markerstrokewidth=0,
xaxis=false,yaxis=false,bottom_margin = 10px,legend=false)
# size=(600,600))
z = 1.1*exp.(im*2π*(0:7)/8)
annotate!(real.(z), imag.(z), text.(labels,12))


# Mapping compass directions to angles in degrees (0 = East, counterclockwise)
direction_to_angle = Dict(
    "E"  => 0,
    "NE" => 45,
    "N"  => 90,
    "NW" => 135,
    "W"  => 180,
    "SW" => 225,
    "S"  => 270,
    "SE" => 315
)
# Convert the random directions to angles using the dictionary
angles_degrees = [direction_to_angle[dir] for dir in Xj_cat]
angles_radians = deg2rad.(angles_degrees)
R = LinRange(0,1,105)
scatter(angles_radians, R, proj=:polar, ms=4, m=:o,
markercolor=color_labels, fa=1,
markerstrokewidth=0,
xaxis=false,yaxis=false,bottom_margin = 10px,legend=false)
# size=(600,600))
z = 1.1*exp.(im*2π*(0:7)/8)
annotate!(real.(z), imag.(z), text.(labels,12))


X_cl2 = Xj_cat[cls .== 2]
X_cl1 = Xj_cat[cls .== 1]
X_cl3 = Xj_cat[cls .== 3]
X_cl4 = Xj_cat[cls .== 4]

begin
R = LinRange(0,1,length(X_cl1))
angles_degrees = [direction_to_angle[dir] for dir in X_cl1]
angles_radians = deg2rad.(angles_degrees)
p1 = scatter(angles_radians, R, proj=:polar, ms=4, m=:o,title="cluster 1",
markercolor=color_palette[1],markerstrokewidth=0,
xaxis=false,yaxis=false,
bottom_margin = 10px,
left_margin = -20px,
legend=false,
size=(350,350),dpi=300)
z = 1.1*exp.(im*2π*(0:7)/8)
annotate!(real.(z), imag.(z), text.(labels,9)) savefig("Xcat_cluster1.pdf")
end


R = LinRange(0,1,length(X_cl2))
angles_degrees = [direction_to_angle[dir] for dir in X_cl2]
angles_radians = deg2rad.(angles_degrees)
p2 = scatter(angles_radians, R, proj=:polar, ms=4, m=:o,title="cluster 2",
markercolor=color_palette[2],markerstrokewidth=0,
xaxis=false,yaxis=false,
bottom_margin = 10px,
left_margin = -20px,
legend=false,
size=(350,350),dpi=300)
z = 1.1*exp.(im*2π*(0:7)/8)
annotate!(real.(z), imag.(z), text.(labels,9)) savefig("Xcat_cluster2.pdf")

R = LinRange(0,1,length(X_cl3))
angles_degrees = [direction_to_angle[dir] for dir in X_cl3]
angles_radians = deg2rad.(angles_degrees)
p3 = scatter(angles_radians, R, proj=:polar, ms=4, m=:o,title="cluster 3",
markercolor=color_palette[3],markerstrokewidth=0,
xaxis=false,yaxis=false,
bottom_margin = 10px,
left_margin = -20px,
legend=false,
size=(350,350),dpi=300)
z = 1.1*exp.(im*2π*(0:7)/8)
annotate!(real.(z), imag.(z), text.(labels,9)) savefig("Xcat_cluster3.pdf")

R = LinRange(0,1,length(X_cl4))
angles_degrees = [direction_to_angle[dir] for dir in X_cl4]
angles_radians = deg2rad.(angles_degrees)
p4 = scatter(angles_radians, R, proj=:polar, ms=4, m=:o,title="cluster 4",
markercolor=color_palette[4],markerstrokewidth=0,
xaxis=false,yaxis=false,
bottom_margin = 10px,
left_margin = -20px,
legend=false,
size=(350,350),dpi=300)
z = 1.1*exp.(im*2π*(0:7)/8)
annotate!(real.(z), imag.(z), text.(labels,9)) savefig("Xcat_cluster4.pdf")

l = @layout [a b; c d]
plot(p1, p2, p3, p4, layout = l,dpi=300, size=(600,500))
savefig("all_Xcat.pdf")


begin
lg = true
alphas = range(0, 2.0, length=1000)
similarity1_vals_cl1 = [similarity1(X_cl1, alpha, lg) for alpha in alphas]
similarity1_vals_cl2 = [similarity1(X_cl2, alpha, lg) for alpha in alphas]
similarity1_vals_cl3 = [similarity1(X_cl3, alpha, lg) for alpha in alphas]
similarity1_vals_cl4 = [similarity1(X_cl4, alpha, lg) for alpha in alphas]
plot(alphas,  similarity1_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(alphas, similarity1_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(alphas, similarity1_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(alphas, similarity1_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
# xlabel!("phi values")
xlabel!(L"$\alpha$ values")
ylabel!("")
title!("Similarity 1 (log=$lg)")
# savefig("similarity1_cat_lg$lg.pdf")
end


begin
lg = true
alphas = range(0, 1.0, length=1000)
similarity2_vals_cl1 = [similarity2(X_cl1, alpha,0, lg) for alpha in alphas]
similarity2_vals_cl2 = [similarity2(X_cl2, alpha,0, lg) for alpha in alphas]
similarity2_vals_cl3 = [similarity2(X_cl3, alpha,0, lg) for alpha in alphas]
similarity2_vals_cl4 = [similarity2(X_cl4, alpha,0, lg) for alpha in alphas]
plot(alphas,  similarity2_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(alphas, similarity2_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(alphas, similarity2_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(alphas, similarity2_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
# xlabel!("phi values")
xlabel!(L"$\alpha$ values")
ylabel!("")
title!("Similarity 2 (log=$lg)")
# savefig("similarity2_cat_lg$lg.pdf")
end


begin
lg = true
alphas = range(0, 3.0, length=1000)
similarity3_vals_cl1 = [similarity3(X_cl1, alpha,0, lg) for alpha in alphas]
similarity3_vals_cl2 = [similarity3(X_cl2, alpha,0, lg) for alpha in alphas]
similarity3_vals_cl3 = [similarity3(X_cl3, alpha,0, lg) for alpha in alphas]
similarity3_vals_cl4 = [similarity3(X_cl4, alpha,0, lg) for alpha in alphas]
plot(alphas,  similarity3_vals_cl1,  label="cluster 1",  lw=2, color=color_palette[1],
dpi=300, size=(760,400)
)
plot!(alphas, similarity3_vals_cl2,  label="cluster 2", lw=2, color=color_palette[2])
plot!(alphas, similarity3_vals_cl3,  label="cluster 3", lw=2, color=color_palette[3])
plot!(alphas, similarity3_vals_cl4,  label="cluster 4", lw=2, color=color_palette[4])
# xlabel!("phi values")
xlabel!(L"$\alpha$ values")
ylabel!("")
title!("Similarity 3 (log=$lg)")
# savefig("similarity3_cat_lg$lg.pdf")
end
