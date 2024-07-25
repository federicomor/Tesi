##################
##   RELABEL    ##
##################

# dont overwrite Si - ignore corollary variables
function relabel(Si::Vector{Int}, n::Int)
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

# overwrite Si - ignore the corollary variables
function relabel!(Si::Vector{Int}, n::Int)
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

# dont overwrite Si - deal with the corollary variables
function relabel_full(Si::Vector{Int}, n::Int)
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

# dont overwrite Si - deal with the corollary variables
function relabel_full!(Si::Vector{Int}, n::Int, Sirelab::Vector{Int}, nhrelab::Vector{Int}, oldLab::Vector{Int})
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




############################################
##   SIMILARITY FUNCTIONS (covariates)    ##
############################################






