nh_tmp = rand((-2:2),10)
using BenchmarkTools

function method1(nh_tmp)
	nclus_temp = count(x->(x>0),nh_tmp)
	return nclus_temp
end
function method2(nh_tmp)
	nclus_temp = sum(nh_tmp .> 0)
	return nclus_temp
end

method2(nh_tmp)
method2(nh_tmp)

@btime method2(nh_tmp)
@btime method1(nh_tmp)

using Random
Random.seed!(1)
function f()
lg_weights = rand(10)
max_ph = maximum(lg_weights)
sum_ph = 0.0
# exponentiate...
@simd for k in eachindex(lg_weights)
	 # for numerical purposes we subract max_ph
	lg_weights[k] = exp(lg_weights[k] - max_ph)
	sum_ph += lg_weights[k]
end
# sum_ph
# lg_weights
end

@btime f()

#######################################
using Random
Random.seed!(1)
sp1 = rand(10)
indexes = [1,2,3,4,6,8]

function f_test(indexes,j)
	sp1_red = sp1[indexes]
	# sp1_red = @view sp1[indexes]
	result = 0.
	for k in 1:10
		aux_idxs = randsubseq(1:length(sp1_red),0.5)
		s1o = @view sp1_red[aux_idxs]
		# copy!(s1o,sp1_red[aux_idxs])
		s1n = copy(s1o); push!(s1n, sp1[j])
		# println(typeof(s1n))
		# resize!(s1n,length(s1o)+1)
		# copy!(s1n,s1o)
		# s1n[end] = sp1[j]
		# copy!(s1n,s1o); push!(s1n, sp1[j])
		result += (sum(s1o)-sum(s1n))*k
	end
	return result
end
function f_test_prealloc(indexes,j)
	sp1_red = sp1[indexes]
	s1n = zeros(length(sp1))
	# sp1_red = @view sp1[indexes]
	result = 0.
	for k in 1:10
		aux_idxs = randsubseq(1:length(sp1_red),0.5)
		s1o = @view sp1_red[aux_idxs]
		copy!(s1n,s1o); push!(s1n, sp1[j])
		result += (sum(s1o)-sum(s1n))*k
	end
	return result
end
f_test(indexes,1)
f_test_prealloc(indexes,1)



using BenchmarkTools
nh_tmp = rand(100)
@btime sum(nh_tmp .> 0)
@btime count(x->(x>0), nh_tmp)


n = 50; rho_tmp = rand((1:5),n); k = 1
include("../utils.jl")
@btime findall(j -> rho_tmp[j]==k, 1:n) 
@btime findall(rho_tmp .== k)
@btime findall_faster(j -> rho_tmp[j]==k, 1:n)

# v = rand(20)
# @btime begin
# 	extra = extrema(v)
# 	R = extra[2] - extra[1]
# end
# @btime begin
# 	R = maximum(v) - minimum(v)
# end


function Gdist(a::Real,b::Real,R::Real)
	a==b && return 1
	return 1-abs(a-b)/R
end
function Gdist(a::String,b::String,R::Real)
	a==b && return 1
	return 0
end
function gower_sim(X_jt::Vector{<:Real})
	S = 0.
	isa(first(X_jt),Real) ? R = maximum(X_jt) - minimum(X_jt) : R = 0.
	cdim = length(v)
	for i in 1:(cdim-1)
		for j in (i+1):cdim
			S += Gdist(X_jt[i],X_jt[j],R)
		end
	end
	return S / (2*(cvdim*(cvdim-1))) # scale to [0,1]

end
v = [1,1000,2,1,0,-1,15,1,1,2] + rand(10)/10
v = rand(10)
v = repeat([1],10)
gower_sim(v)



######################################
using Plots
using PDFmerger: append_pdf!

# Create a PDF to store all plots
pdf_filename = "trace_plots.pdf"
# pdf_backend = pdf(pdf_filename)

# Number of iterations and variables (replace with your settings)
n_iters = 1000
n_vars = 5

# To simulate the trace values (you'll use your own variables in the loop)
variable_traces = zeros(n_iters, n_vars)

# Initialize a plot for each variable
plots = [plot(title="Variable $i", xlabel="Iteration", ylabel="Value") for i in 1:n_vars]

# Iteratively update the plots and add points
for iter in 1:n_iters
    # Simulate updates (you would use your actual variable updates)
    for i in 1:n_vars
        variable_traces[iter, i] = randn()*iter/10 # Replace with the actual draw of variable `i`

        # Update the plot for this variable by adding the new point
        scatter!(plots[i], [iter], [variable_traces[iter, i]], seriestype=:line, color=:blue, label="",markersize=2)

        # Display the current plot (optional, slows down for large n_iters)
        # display(plots[i])
    end
end

p = plot(plots[1:5]...,layout=5)
savefig(p,"trace_tmp.pdf")
append_pdf!("all_trace_plots.pdf", "trace_tmp.pdf", cleanup=true)

for i in 1:n_vars
	savefig(plots[i], "trace_tmp.pdf")
    append_pdf!("all_trace_plots.pdf", "trace_tmp.pdf", cleanup=true)
end


# After the loop, save all plots to the single PDF

################################

function fun!(v::Vector{Float64},aux_array::Vector{Float64})
	aux_array .= 0
	aux_array = rand(10)
	v = aux_array
end
v = zeros(10)
w = zeros(10)
aux_array = zeros(10)
fun!(v,aux_array);
@show v w aux_array
fun!(w,aux_array)
@show v w aux_array