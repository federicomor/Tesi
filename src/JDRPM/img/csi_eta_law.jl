using Plots

function f(x, mu, b)
    term1 = (1/(2*b)) * exp(-1/b * abs(log((x + 1) / (1 - x)) - mu))
    # term2 = (1/(x+1)) + (1/(1-x))
    return term1
end

mu = 0.0
for b in [0.5, 1, 2]
    x_vals = range(-0.999, 0.999, length=500)
    y_vals = [f(x,mu,b) for x in x_vals]

    plot(x_vals, y_vals, xlabel="x", ylabel="f(x)", title="ξᵢ ∼ Laplace(μ=$mu, b=$b)", lw=2, legend=false)
    ylims!(0,1.2)
    savefig("img/csi_eta_law_plot_b_$(replace(string(b),"."=>"_")).svg")
end