
begin
lg_weights = [-1000.0, -2000, -1500, -100, -200]
ph_tmp = copy(lg_weights)
sort!(ph_tmp)
max_ph = ph_tmp[end]
sum_ph = 0.0
# exponentiate...
for k in 1:length(lg_weights)
	 # for numeric purposes we subract max_ph
	# lg_weights[k] = exp(lg_weights[k]) # without
	lg_weights[k] = exp(lg_weights[k]-max_ph) # with
	global sum_ph += lg_weights[k]
end
@show lg_weights
# ... and normalize
for k in 1:length(lg_weights)
	lg_weights[k] /= sum_ph
end
@show lg_weights
end



begin 
	ph = rand((0:100),5) * 1.0
	ph ./= sum(ph)
	uu = rand()
	# uu = 1
	begin cph = cumsum(ph)
		@show ph
		@show cph
			@assert isapprox(cph[end], 1, atol=1e-6)
		println(findfirst(x -> x==1, uu .<= cph))
		end
		@show uu
		begin cph = cumsum(ph)
		iaux = 0
		for k in 1:length(ph)
			if uu < cph[k]
				global iaux = k
				break
			end
		end
		println(iaux)
	end
end


begin
A_star = [0.5695141296842557 -0.5958323106341032; -0.5958323106341032 5.607107634272117]
invA_star = inv(A_star)
# invA_star = [1.9755086272843259 0.20992496432168561; 0.20992496432168564 0.20065248429953192]

using LinearAlgebra
println(eigvals(A_star)) # [0.5000000000009021, 5.67662176395547] all >0
println(eigvals(invA_star)) # [0.1761611115874664, 1.9999999999963913] still all >0

println(isposdef(A_star)) # true
isposdef(invA_star) # false

end