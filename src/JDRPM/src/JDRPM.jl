module JDRPM

export my_example, MCMC_fit, close_log_file

function my_example(iters::Number)
	result = 0.0
	for i in 1:iters
		print("iteration $i of $iters\r")
		if i%20 == 0
			println()
		end
		# do my heavy work
		result += rand() # for example
		sleep(0.02)
	end
	return result
end

include("../MCMC_fit.jl")

end # module JDRPM
