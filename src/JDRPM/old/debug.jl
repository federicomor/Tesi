macro showd(exs...)
	args = []
	for ex in exs
		push!(args, sprint(Base.show_unquoted, ex), " = ", esc(ex), '\n')
	end
	return :(string($(args...)))
end

function tostr(obj)
    io = IOBuffer()
    show(io, "text/plain", obj)
    String(take!(io))
end

function debug(str)
	print(log_file,str * (str[end]=='\n' ? "" : "\n"))
end

function printlgln(str)
	print(log_file,str * (str[end]=='\n' ? "" : "\n"))
end

function pretty_log(str)
	if str=="Si_iter" println(log_file,"Si_iter\n",tostr(Si_iter)); return; end
	if str=="gamma_iter" println(log_file,"gamma_iter\n",tostr(gamma_iter)); return; end
	if str=="nh" println(log_file,"nh\n",tostr(nh)); return; end
	if str=="nclus_iter" println(log_file,"nclus_iter\n",tostr(nclus_iter)); return; end
	if str=="muh_iter" println(log_file,"muh_iter\n",tostr(muh_iter)); return; end
	if str=="sig2h_iter" println(log_file,"sig2h_iter\n",tostr(sig2h_iter)); return; end
	if str=="alpha_iter" println(log_file,"alpha_iter\n",tostr(alpha_iter)); return; end
	if str=="theta_iter" println(log_file,"theta_iter\n",tostr(theta_iter)); return; end
	if str=="tau2_iter" println(log_file,"tau2_iter\n",tostr(tau2_iter)); return; end
	if str=="phi0_iter" println(log_file,"phi0_iter\n",tostr(phi0_iter)); return; end
	if str=="phi1_iter" println(log_file,"phi1_iter\n",tostr(phi1_iter)); return; end
	if str=="lambda2_iter" println(log_file,"lambda2_iter\n",tostr(lambda2_iter)); return; end
	if str=="sp_coords" println(log_file,"sp_coords\n",tostr(sp_coords)); return; end
	if str=="sp1" println(log_file,"sp1\n",tostr(sp1)); return; end
	if str=="sp2" println(log_file,"sp2\n",tostr(sp2)); return; end
	if str=="beta_iter" println(log_file,"beta_iter\n",tostr(beta_iter)); return; end
end

# function section(title::String)
# 	total_width = length(title) + 4
# 	debug("┌" * "─" ^ (total_width - 2) * "┐")
# 	debug("│" * " " * title * " " * "│")
# 	debug("└" * "─" ^ (total_width - 2) * "┘")
# end