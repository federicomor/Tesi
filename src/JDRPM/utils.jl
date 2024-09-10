using SpecialFunctions
using Statistics

logit(x::Real) = log(x / (one(x) - x))
const logpi = log(π)
const log2pi = log(2*π)

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

# dont overwrite Si - consider corollary variables
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


# dont overwrite Si - consider corollary variables
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
	sbar = [mean(s1), mean(s2)]
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
end

# below there is the one "translated" from C, which is way more efficient
function cohesion3(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	# Compute sample means
	sbar1 = mean(s1)
	sbar2 = mean(s2)
	# Compute deviations from the sample mean
	S1, S2, S3, S4 = 0.0, 0.0, 0.0, 0.0
	@inbounds for i in 1:sdim
		s_sbar1 = s1[i] - sbar1
		s_sbar2 = s2[i] - sbar2

		S1 += s_sbar1 * s_sbar1
		S4 += s_sbar2 * s_sbar2
		S2 += s_sbar1 * s_sbar2
	end
	S3 = copy(S2) # to avoid repeating computations

	# Updated parameters for cohesion 3
	kn = k0 + sdim
	vn = v0 + sdim

	auxvec1_1 = sbar1 - mu_0[1]
	auxvec1_2 = sbar2 - mu_0[2]

	auxmat1_1 = auxvec1_1^2
	auxmat1_2 = auxvec1_1 * auxvec1_2
	auxmat1_3 = copy(auxmat1_2)
	auxmat1_4 = auxvec1_2^2

	auxconst1 = k0 * sdim
	auxconst2 = k0 + sdim
	Psi_n_1 = Psi[1] + S1 + auxconst1 / (auxconst2) * auxmat1_1
	Psi_n_2 = Psi[2] + S2 + auxconst1 / (auxconst2) * auxmat1_2
	Psi_n_3 = Psi[3] + S3 + auxconst1 / (auxconst2) * auxmat1_3
	Psi_n_4 = Psi[4] + S4 + auxconst1 / (auxconst2) * auxmat1_4

	detPsi_n = Psi_n_1 * Psi_n_4 - Psi_n_2 * Psi_n_3
	detPsi = Psi[1] * Psi[4] - Psi[2] * Psi[3]

	out = -sdim * logpi + G2a(0.5 * vn, true) - G2a(0.5 * v0, true) + 0.5 * v0 * log(detPsi) - 0.5 * vn * log(detPsi_n) + log(k0) - log(kn)
	return lg ? out : exp(out)

end

function cohesion4(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, mu_0::AbstractVector{Float64}, k0::Real, v0::Real, Psi::Matrix{Float64}; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	# Compute sample means
	sbar1 = mean(s1)
	sbar2 = mean(s2)
	# Compute deviations from the sample mean
	S1, S2, S3, S4 = 0.0, 0.0, 0.0, 0.0
	@inbounds for i in 1:sdim
		s_sbar1 = s1[i] - sbar1
		s_sbar2 = s2[i] - sbar2

		S1 += s_sbar1 * s_sbar1
		S4 += s_sbar2 * s_sbar2
		S2 += s_sbar1 * s_sbar2
	end
	S3 = copy(S2) # to avoid repeating computations

	# Updated parameters for cohesion 3
	kn = k0 + sdim
	vn = v0 + sdim

	auxvec1_1 = sbar1 - mu_0[1]
	auxvec1_2 = sbar2 - mu_0[2]

	auxmat1_1 = auxvec1_1^2
	auxmat1_2 = auxvec1_1 * auxvec1_2
	auxmat1_3 = copy(auxmat1_2)
	auxmat1_4 = auxvec1_2^2

	auxconst1 = k0 * sdim
	auxconst2 = k0 + sdim
	Psi_n_1 = Psi[1] + S1 + auxconst1 / (auxconst2) * auxmat1_1
	Psi_n_2 = Psi[2] + S2 + auxconst1 / (auxconst2) * auxmat1_2
	Psi_n_3 = Psi[3] + S3 + auxconst1 / (auxconst2) * auxmat1_3
	Psi_n_4 = Psi[4] + S4 + auxconst1 / (auxconst2) * auxmat1_4

	detPsi_n = Psi_n_1 * Psi_n_4 - Psi_n_2 * Psi_n_3
	
	# Updated parameters for cohesion 4
	knn = kn + sdim
	vnn = vn + sdim

	mu_n_1 = (k0 * mu_0[1] + sdim * sbar1) / (auxconst2)
	mu_n_2 = (k0 * mu_0[2] + sdim * sbar2) / (auxconst2)
	sbar_mun1 = sbar1 - mu_n_1
	sbar_mun2 = sbar2 - mu_n_2

	auxmat2_1 = sbar_mun1^2
	auxmat2_2 = sbar_mun1 * sbar_mun2
	auxmat2_3 = copy(auxmat2_2)
	auxmat2_4 = sbar_mun2^2

	auxconst3 = kn * sdim
	auxconst4 = kn + sdim
	Psi_nn_1 = Psi_n_1 + S1 + auxconst3 / (auxconst4) * auxmat2_1
	Psi_nn_2 = Psi_n_2 + S2 + auxconst3 / (auxconst4) * auxmat2_2
	Psi_nn_3 = Psi_n_3 + S3 + auxconst3 / (auxconst4) * auxmat2_3
	Psi_nn_4 = Psi_n_4 + S4 + auxconst3 / (auxconst4) * auxmat2_4
	detPsi_nn = Psi_nn_1 * Psi_nn_4 - Psi_nn_2 * Psi_nn_3
	
	out = -sdim * logpi + G2a(0.5 * vnn, true) - G2a(0.5 * vn, true) + 0.5 * vn * log(detPsi_n) - 0.5 * vnn * log(detPsi_nn) + log(kn) - log(knn)
	return lg ? out : exp(out)
end


# paper 6 pag 4, cluster variance/entropy similarity function
function cohesion5(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, phi::Real; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	# compute the centroids
	cent1 = mean(s1)
	cent2 = mean(s2)
	# compute the sum of the distances
	sum_dist = sum(sqrt((s1[i] - cent1)^2 + (s2[i] - cent2)^2) for i in 1:sdim)
		
	out = -phi*sum_dist
	return lg ? out : exp(out)
end

# non trovata su nessun paper
function cohesion6(s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, phi::Real; lg::Bool, M::Real=1.0)
	sdim = length(s1)
	if sdim==1
		return lg ? 0.0 : 1.0
	end
	# compute the centroids
	cent1 = mean(s1)
	cent2 = mean(s2)
	# compute the sum of the distances
	sum_dist = sum(sqrt((s1[i] - cent1)^2 + (s2[i] - cent2)^2) for i in 1:sdim)
	
	out = -phi*log(sum_dist)
	return lg ? out : exp(out)
end

# n = 20
# s1 = rand(n)
# s2 = rand(n)
# alpha = 0.13
# a = 10
# mu_0 = [1.,2.]
# k0 = 0.5
# v0 = 0.5
# Psi = [1. 2.; 2 4]
# phi=0.5

# cohesion1(s1,s2,alpha,lg=false)
# cohesion2(s1,s2,a,lg=false)
# cohesion3(s1, s2, mu_0, k0, v0, Psi,lg=false)
# cohesion4(s1, s2, mu_0, k0, v0, Psi,lg=false)
# cohesion5(s1,s2,phi,lg=false)
# cohesion6(s1,s2,phi,lg=false)

function spatial_cohesion(idx::Real, s1::AbstractVector{Float64}, s2::AbstractVector{Float64}, sp_params::Vector; lg::Bool, M::Real)
	idx==1.0 && return cohesion1(s1,s2,sp_params[1],lg=lg,M=M) 
	idx==2.0 && return cohesion2(s1,s2,sp_params[1],lg=lg,M=M) 
	idx==3.0 && return cohesion3(s1,s2,sp_params[1],sp_params[2],sp_params[3],sp_params[4],lg=lg,M=M) 
	idx==4.0 && return cohesion3(s1,s2,sp_params[1],sp_params[2],sp_params[3],sp_params[4],lg=lg,M=M) 
	idx==5.0 && return cohesion5(s1,s2,sp_params[1],lg=lg,M=M) 
	idx==6.0 && return cohesion6(s1,s2,sp_params[1],lg=lg,M=M) 
end

############################################
##   SIMILARITY FUNCTIONS (covariates)    ##
############################################

# farle singole, concentrate su una singola covariata per volta
# e come sopra con la parte spaziale, in modo cioè sceglibile/flessibile
# normalizzare le covariate

# eg X_jt would be X_cl[indexes,covariate_index,t] 

##### Cluster variance/entropy similarity function - paper 6 pag 4
function similarity1(X_jt::Union{AbstractVector{<:Real}, AbstractVector{String}}, alpha::Real; lg::Bool)
	if isa(first(X_jt),Real) # numerical case
		xbar_j = mean(X_jt)
		card_Sjt = length(X_jt)
		Hx = 0.0
		@inbounds for i in eachindex(X_jt)
			Hx += (X_jt[i] - xbar_j)^2 
		end
		Hx /= card_Sjt
		return lg ? -alpha * Hx : exp(-alpha * Hx) 

	else # categorical case
		# X_jt = ["Old", "Young", "Middle", "Young", "Old", "Old"]
		unique_keys = unique(X_jt)
		# @show unique_keys
		counts = Dict(key => 0.0 for key in unique_keys)
		for item in X_jt
			counts[item] += 1
		end
		# @show counts
		total = length(X_jt)
		for key in keys(counts)
			counts[key] /= total
		end
		Hx = 0.0
		for r in keys(counts)
			Hx -= counts[r] * log(counts[r])
		end
		return lg ? -alpha * Hx : exp(-alpha * Hx) 
	end
end
# similarity1([repeat(["a"],1000)...,"b"], 5, lg=false)


##### Gower dissimilarity
function gower_d(x1::Real, x2::Real, R::Real) # numerical case
	return abs(x1-x2)/R
end
function gower_d(x1::String, x2::String, R::Real) # categorical case
	# R is fictitious here, just to allow the same function call later
	return Int(!(x1 == x2))
	# add one when values are different => higher values for dissimilar data
	# so eg fully dissimilar values will give around zero as value returned
	# which i.e. would be -Inf with the lg=true of our algorithm
	# which seems right, i.e. is all in a "penalization" perspective now
	# where we just penalize less the vectors which are similar,
	# rahter than givine them a positive high value
end

##### Total Gower dissimilarity - paper 6 pag 4
function similarity2(X_jt::Union{AbstractVector{<:Real}, AbstractVector{String}}, alpha::Real; lg::Bool)
	H = 0.0
	if isa(first(X_jt),Real) 
		R = (maximum(X_jt) - minimum(X_jt))
		if R==0.0 R=eps() end
	else R=0.0 end
	n_j = length(X_jt)
	for l in 1:(n_j-1)
		for k in (l+1):n_j
			H += gower_d(X_jt[l], X_jt[k], R)
		end
	end 
	out = -alpha * H
	return lg ? out : exp(out)
end

##### Average Gower dissimilarity - paper 6 pag 4
function similarity3(X_jt::Union{AbstractVector{<:Real}, AbstractVector{String}}, alpha::Real; lg::Bool)
	H = 0.0
	if isa(first(X_jt),Real) 
		R = (maximum(X_jt) - minimum(X_jt))
		if R==0.0 R=eps() end
	else R=0.0 end
	n_j = length(X_jt)
	for l in 1:(n_j-1)
		for k in l:n_j
			H += gower_d(X_jt[l], X_jt[k], R)
		end
	end 
	out = -2*alpha / (n_j*(n_j-1)) * H
	return lg ? out : exp(out)
end

# # 0 => completely dissimilar
# # 1 => completely similar

# X_jt = [1,1,1]
# similarity2(X_jt,1,lg=false)
# similarity3(X_jt,1,lg=false)
# XX_jt = [1,1,2]
# similarity2(XX_jt,1,lg=false)
# similarity3(XX_jt,1,lg=false)

# X_jt = [repeat([1],5)...,2]
# similarity2(X_jt,1,lg=false)
# similarity3(X_jt,1,lg=false)
# similarity3(X_jt,2,lg=false)
# X_jt = collect(1:201)
# similarity2(X_jt,1,lg=false) # 0 ie completely dissimilar
# similarity3(X_jt,1,lg=false) # strangely high, better to increase alpha
# similarity3(X_jt,2,lg=false)
# similarity3(X_jt,10,lg=false)

# XX_jt = [1,1,1,1,1,1,7,8,9,10,11]
# similarity2(XX_jt,1,lg=false)
# similarity3(XX_jt,3,lg=false)
# XX_jt = [1,2,3,4,5,6,7,8,9,10,11]
# similarity2(XX_jt,1,lg=false)
# similarity3(XX_jt,3,lg=false)


function similarity4(X_jt::AbstractVector{<:Real}, mu_c::Real, lambda_c::Real, a_c::Real, b_c::Real; lg::Bool)
	n = length(X_jt)
	nm = n/2
	xbar = mean(X_jt)
	aux2 = 0.0
	@inbounds @fastmath @simd for i in eachindex(X_jt)
		aux2 += X_jt[i]^2
	end
	# @show aux2
	aux1 = b_c + 0.5 * (aux2 - (n*xbar + lambda_c*mu_c)^2/(n+lambda_c) + lambda_c*mu_c^2 )
	out = -nm*log2pi + 0.5*log(lambda_c/(lambda_c+n)) + lgamma(a_c+nm) - lgamma(a_c) + a_c*log(b_c) + (-a_c-nm)*log(aux1)
	return lg ? out : exp(out)
end


function similarity5(X_jt::AbstractVector{<:Real}, mu_c::Real, lambda_c::Real, a_c::Real, b_c::Real; lg::Bool)
	return lg ? 0.0 : 1.0
end


# 1,2,3 work on both covariate types
# 4 and 5 only on numerical

# numerical covariates specialization
function covariate_similarity(idx::Real, X_jt::AbstractVector{<:Real}, cv_params::Vector; lg::Bool)
	# println("numerical")
	idx==1.0 && return similarity1(X_jt,cv_params[1],lg=lg) 
	idx==2.0 && return similarity2(X_jt,cv_params[1],lg=lg) 
	idx==3.0 && return similarity3(X_jt,cv_params[1],lg=lg) 
	idx==4.0 && return similarity4(X_jt,cv_params[1],cv_params[2],cv_params[3],cv_params[4],lg=lg) 
	idx==5.0 && return similarity5(X_jt,cv_params[1],cv_params[2],cv_params[3],cv_params[4],lg=lg) 
end

# categorical covariates specialization
function covariate_similarity(idx::Real, X_jt::AbstractVector{String}, cv_params::Vector; lg::Bool)
	# println("categorical")
	idx==1.0 && return similarity1(X_jt,cv_params[1],lg=lg) 
	idx==2.0 && return similarity2(X_jt,cv_params[1],lg=lg) 
	idx==3.0 && return similarity3(X_jt,cv_params[1],lg=lg) 
end

# mu_c = 0.0
# lambda_c = 1.0
# a_c = 2.0
# b_c = 2.0
# cv_params = [mu_c,lambda_c,a_c,b_c]
# X_jt = [1,1,2,2,2,1,2,1,1,1] ./ 10
# covariate_similarity(4,X_jt,cv_params,lg=true) # -7
# X_jt = [1,1,2,2,2,1,2,1,1,100] ./ 10
# covariate_similarity(4,X_jt,cv_params,lg=true) # -29, correctly more penalized

# covariate_similarity(1,[repeat(["a"],100000)...,"b"], [2], lg=true) # -0.0002
# covariate_similarity(1,[repeat(["a"],10)...,"b"], [2], lg=true) # -0.609, correctly more penalized
