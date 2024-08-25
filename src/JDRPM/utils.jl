using SpecialFunctions
using LinearAlgebra
using Statistics
using StaticArrays

logit(x::Real) = log(x / (one(x) - x))
const logpi = log(π)

# function section(title::String)
# 	total_width = length(title) + 4
# 	debug("┌" * "─" ^ (total_width - 2) * "┐")
# 	debug("│" * " " * title * " " * "│")
# 	debug("└" * "─" ^ (total_width - 2) * "┘")
# end

##################
##   RELABEL    ##
##################

# dont overwrite Si - ignore corollary variables
function relabel(Si::AbstractVector{Int}, n::Int)
	Sirelab = zeros(Int,n)
	shuffle = n
	loc = 1
	lab = 1 # new label index

	while shuffle > 0
		for j in 1:n
			if Si[j] == Si[loc]
				Sirelab[j] = lab
				shuffle -= 1
			end
		end
		lab += 1
		loc = findfirst(x -> Sirelab[x] == 0, 1:n) # find the next unprocessed label
	end
	return Sirelab
end
relabel(Si::AbstractVector{Int}) = relabel(Si,length(Si))

# overwrite Si - ignore the corollary variables
function relabel!(Si::AbstractVector{Int}, n::Int)
	Sirelab = zeros(Int,n)
	shuffle = n
	loc = 1
	lab = 1 # new label index

	while shuffle > 0
		for j in 1:n
			if Si[j] == Si[loc]
				Sirelab[j] = lab
				shuffle -= 1
			end
		end
		lab += 1
		loc = findfirst(x -> Sirelab[x] == 0, 1:n) # find the next unprocessed label
	end
	for j in 1:n
		Si[j] = Sirelab[j]
	end
end
relabel!(Si::AbstractVector{Int}) = relabel!(Si,length(Si))

# dont overwrite Si - deal with the corollary variables
function relabel_full(Si::AbstractVector{Int}, n::Int)
	Sirelab = zeros(Int,n)
	nhrelab = zeros(Int,n)
	oldLab = zeros(Int,n)
	shuffle = n
	loc = 1
	lab = 1 # new label index

	while shuffle > 0
		for j in 1:n
			if Si[j] == Si[loc]
				Sirelab[j] = lab
				shuffle -= 1
				nhrelab[lab] += 1
			end
		end
		oldLab[lab] = Si[loc]
		lab += 1
		loc = findfirst(x -> Sirelab[x] == 0, 1:n) # find the next unprocessed label
	end
	return Sirelab, nhrelab, oldLab
end
relabel_full(Si::AbstractVector{Int}) = relabel_full(Si,length(Si))


# dont overwrite Si - deal with the corollary variables
function relabel_full!(Si::AbstractVector{Int}, n::Int, Sirelab::Vector{Int}, nhrelab::Vector{Int}, oldLab::Vector{Int})
	fill!(Sirelab, 0)
	fill!(nhrelab, 0)
	fill!(oldLab, 0)
	shuffle = n
	loc = 1
	lab = 1 # new label index

	while shuffle > 0
		for j in 1:n
			if Si[j] == Si[loc]
				Sirelab[j] = lab
				shuffle -= 1
				nhrelab[lab] += 1
			end
		end
		oldLab[lab] = Si[loc]
		lab += 1
		loc = findfirst(x -> Sirelab[x] == 0, 1:n) # find the next unprocessed label
	end
end


# n = 10
# Si = rand((1:5),n)
# Sirelab = zeros(Int, n)
# nhrelab = zeros(Int, n)
# oldLab = zeros(Int, n)
# @timev relabel!(Si, n, Sirelab, nhrelab, oldLab)
# println("""
# 		 Si = $Si
# 	Sirelab = $Sirelab
# 	nhrelab = $nhrelab
# 	 oldLab = $oldLab
# 	""")




########################
##   COMPATIBILITY    ##
########################

function compatibility(rho1::Vector{Int}, rho2::Vector{Int})
	n = length(rho1)
	scr1 = relabel(rho1, n)
	scr2 = relabel(rho2, n)
	# check
	for i in 1:n
		if scr1[i] != scr2[i]
			return 0
		end
	end
	return 1
end

# rho1 = [1, 2, 1, 3, 2]
# rho2 = [1, 1, 2, 3, 2]
# rho3 = [3, 1, 3, 2, 1]
# println(rho1, "\n", rho2, " -> ", compatibility(rho1,rho2))
# println(rho1, "\n", rho2, " -> ", compatibility(rho1,rho3))



#####################################
##   COHESION FUNCTIONS (space)    ##
#####################################

# paper 3 section 3.1
function cohesion1(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, alpha::Real; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	if sdim==1 
		return lg ? log(M) : M
	end
	# out = log(M) + lgamma(sdim)
	out = 0.0
	# compute the centroids
	cent1 = mean(s1)
	cent2 = mean(s2)
	# compute the sum of the distances (the D_h in the paper)
	sum_dist = sum(sqrt((s1[i] - cent1)^2 + (s2[i] - cent2)^2) for i in 1:sdim)

	# decide what to return
	if sum_dist >= 1
		out -= lgamma(alpha*sum_dist) # minus since we are at the denominator
	elseif sum_dist != 0
		out -= log(sum_dist) # minus since we are at the denominator
	else
		out = log(1)
	end
	return lg ? out : exp(out)
end

# paper 3 section 3.1
function cohesion2(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, a::Real; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	# out = log(M) + lgamma(sdim)
	out = 1.0
	for i in 1:sdim
		for j in 1:sdim
			dist = sqrt((s1[i] - s1[j])^2 + (s2[i] - s2[j])^2)
			if dist > a
				return lg ? log(0.0) : 0.0 # in this case the rest vanishes
			end
		end
	end
	return lg ? out : exp(out)
end

function G2a(a::Real, lg::Bool)
	out = logpi + lgamma(a) + lgamma(a - 0.5)
	return lg ? out : exp(out)
end

# paper 3 section 3.1
function cohesion3_4(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}; Cohesion::Int, lg::Bool, M::Real=1.0)
	sdim = length(s1)
	sp = [s1 s2]
	sbar = @SVector [mean(s1), mean(s2)]
	# S = sum( (sp[i,:] - sbar)*(sp[i,:] - sbar)' for i in 1:sdim) # sum is slow
	# FIXED: sum is slow because i didnt initialize S
	# S = zeros(2,2); S = sum( (sp[i,:] - sbar)*(sp[i,:] - sbar)' for i in 1:sdim) # this is fast now
	S = zeros(2,2)
	for i in 1:sdim
		vtemp1 = sp[i,:] - sbar
		S += (vtemp1)*(vtemp1)'
	end
	
	# compute updated parameters
	# for cohesion 3
	kn = k0+sdim
	vn = v0+sdim
	vtemp2 = sbar-mu_0
	Psi_n = Psi + S + (k0*sdim)/(k0+sdim)*vtemp2*vtemp2'
	
	if Cohesion == 3
		out = -sdim * logpi + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * logdet(Psi) - 0.5 * vn * logdet(Psi_n) + log(k0) - log(kn)
		return lg ? out : exp(out)
	end

	# for cohesion 4
	knn = kn+sdim
	vnn = vn+sdim
	mu_n = (k0*mu_0 + sdim*sbar)/(k0+sdim)
	Psi_nn = Psi_n + S + (kn*sdim)/(kn+sdim)*(sbar-mu_n)*(sbar-mu_n)'
	
	if Cohesion == 4
		out = -sdim * logpi + G2a(0.5 * vnn, true) - G2a(0.5 * vn, true) + 0.5 * vn * logdet(Psi_n) - 0.5 * vnn * logdet(Psi_nn) + log(kn) - log(knn)
		return lg ? out : exp(out)
	end

	return NaN
end


##################
# below there is the one "translated" from C, but the above is more readable/efficient/easy
# but are the same, test them on these if you want
# (s1, s2, Psi, Psi_vec) = ([0.06541394973925674, 0.1875903839556078, 0.7065742551867602, 0.8223492385591462], [0.6711654699192571, 0.6199278925430733, 0.36880242735326396, 0.9723482028752322], [2.0 1.0; 1.0 3.0], [2.0, 1.0, 1.0, 3.0])
##################
function cohesion3_4_C(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}; Cohesion::Int, lg::Bool, M::Real=1.0)
		dim = length(s1)
		# Compute sample means
		sbar1 = mean(s1)
		sbar2 = mean(s2)
		# Compute deviations from the sample mean
		Vs1, Vs2, Vs3, Vs4 = 0.0, 0.0, 0.0, 0.0
		for i in 1:dim
			s_sbar1 = s1[i] - sbar1
			s_sbar2 = s2[i] - sbar2
			Vs1 += s_sbar1 * s_sbar1
			Vs2 += s_sbar1 * s_sbar2
			Vs3 += s_sbar2 * s_sbar1
			Vs4 += s_sbar2 * s_sbar2
		end
		# Updated parameters
		kn = k0 + dim
		vn = v0 + dim
		knn = kn + dim
		vnn = vn + dim
		mun1 = k0 / (k0 + dim) * mu_0[1] + dim / (k0 + dim) * sbar1
		mun2 = k0 / (k0 + dim) * mu_0[2] + dim / (k0 + dim) * sbar2
		sbar_mu_01 = sbar1 - mu_0[1]
		sbar_mu_02 = sbar2 - mu_0[2]
		sbar_mun1 = sbar1 - mun1
		sbar_mun2 = sbar2 - mun2
		Vsbarmu_01 = sbar_mu_01^2
		Vsbarmu_02 = sbar_mu_01 * sbar_mu_02
		Vsbarmu_03 = copy(Vsbarmu_02)
		Vsbarmu_04 = sbar_mu_02^2
		Vsbarmun1 = sbar_mun1^2
		Vsbarmun2 = sbar_mun1 * sbar_mun2
		Vsbarmun3 = copy(Vsbarmun2)
		Vsbarmun4 = sbar_mun2^2
		Ln1 = Psi[1] + Vs1 + k0 * dim / (k0 + dim) * Vsbarmu_01
		Ln2 = Psi[2] + Vs2 + k0 * dim / (k0 + dim) * Vsbarmu_02
		Ln3 = Psi[3] + Vs3 + k0 * dim / (k0 + dim) * Vsbarmu_03
		Ln4 = Psi[4] + Vs4 + k0 * dim / (k0 + dim) * Vsbarmu_04
		Lnn1 = Ln1 + Vs1 + kn * dim / (kn + dim) * Vsbarmun1
		Lnn2 = Ln2 + Vs2 + kn * dim / (kn + dim) * Vsbarmun2
		Lnn3 = Ln3 + Vs3 + kn * dim / (kn + dim) * Vsbarmun3
		Lnn4 = Ln4 + Vs4 + kn * dim / (kn + dim) * Vsbarmun4
		dPsi = Psi[1] * Psi[4] - Psi[2] * Psi[3]
		dLn = Ln1 * Ln4 - Ln2 * Ln3
		dLnn = Lnn1 * Lnn4 - Lnn2 * Lnn3
		if Cohesion == 3
			out = -dim * log(π) + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * log(dPsi) - 0.5 * vn * log(dLn) + log(k0) - log(kn)
		elseif Cohesion == 4
			out = -dim * log(π) + G2a(0.5 * vnn, true) - G2a(0.5 * vn, true) + 0.5 * vn * log(dLn) - 0.5 * vnn * log(dLnn) + log(kn) - log(knn)
		end
		return lg ? out : exp(out)
end

# paper 6 pag 4, cluster variance/entropy similarity function
function cohesion5(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, phi::Real; lg::Bool, M::Real=1.0)
	dim = length(s1)
	# compute the centroids
	cent1 = mean(s1)
	cent2 = mean(s2)
	# compute the sum of the distances
	sum_dist = sum(sqrt((s1[i] - cent1)^2 + (s2[i] - cent2)^2) for i in 1:dim)
		
	out = -phi*sum_dist
	return lg ? out : exp(out)
end

# non trovata su nessun paper
function cohesion6(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, phi::Real; lg::Bool, M::Real=1.0)
	dim = length(s1)
	if dim==1
		return lg ? 0.0 : 1.0
	end
	# compute the centroids
	cent1 = mean(s1)
	cent2 = mean(s2)
	# compute the sum of the distances
	sum_dist = sum(sqrt((s1[i] - cent1)^2 + (s2[i] - cent2)^2) for i in 1:dim)
	
	out = -phi*log(sum_dist)
	return lg ? out : exp(out)
end


# s1 = [1,2,3].*1.0
# s2 = [4,5,6].*1.0
# alpha = 0.1
# a = 10
# mu_0 = [1.,2.]
# k0 = 0.5
# v0 = 0.5
# Psi = [1. 2.; 2 4]
# phi=0.5

# x1 = cohesion1(s1,s2,alpha,false)
# x2 = cohesion2(s1,s2,a,false)
# x3 = cohesion3_4(s1, s2, mu_0, k0, v0, Psi,3,false)
# x4 = cohesion3_4(s1, s2, mu_0, k0, v0, Psi,4,false)
# x5 = cohesion5(s1,s2,phi,false)
# x6 = cohesion6(s1,s2,phi,false)
# println()
# y1 = cohesion1(s1,s2,alpha,true)
# y2 = cohesion2(s1,s2,a,true)
# y3 = cohesion3_4(s1, s2, mu_0, k0, v0, Psi,3,true)
# y4 = cohesion3_4(s1, s2, mu_0, k0, v0, Psi,4,true)
# y5 = cohesion5(s1,s2,phi,true)
# y6 = cohesion6(s1,s2,phi,true)

# @show y1 - log(x1)
# @show y2 - log(x2)
# @show y3 - log(x3)
# @show y4 - log(x4)
# @show y5 - log(x5)
# @show y6 - log(x6)

function spatial_cohesion(idx::Real, s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, sp_params; lg::Bool, M::Real)
	idx==1.0 && return cohesion1(s1,s2,sp_params[1],lg=lg,M=M) 
	idx==2.0 && return cohesion2(s1,s2,sp_params[1],lg=lg,M=M) 
	idx==3.0 && return cohesion3_4_C(s1,s2,sp_params[1],sp_params[2],sp_params[3],sp_params[4],lg=lg,M=M,Cohesion=3) 
	idx==4.0 && return cohesion3_4_C(s1,s2,sp_params[1],sp_params[2],sp_params[3],sp_params[4],lg=lg,M=M,Cohesion=4) 
	idx==5.0 && return cohesion5(s1,s2,sp_params[1],lg=lg,M=M) 
	idx==6.0 && return cohesion6(s1,s2,sp_params[1],lg=lg,M=M) 
end

############################################
##   SIMILARITY FUNCTIONS (covariates)    ##
############################################


# farle singole, concentrate su una singola covariata per volta
# e come sopra con la parte spaziale, in modo cioè sceglibile/flessibile
# normalizzare le covariate