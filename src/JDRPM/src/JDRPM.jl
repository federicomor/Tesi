module JDRPM

export my_example, MCMC_fit, close_log_file

function my_example(iters::Number)
	result = 0.0

	# for i in 1:iters
	# 	print("iteration $i of $iters\r")
	# 	if i%20 == 0
	# 		println()
	# 	end
	# 	# do my heavy work
	# 	result += rand() # for example
	# 	sleep(0.02)
	# end
	
	p = Progress(round(Int64(iters));
		showspeed=true,
		dt=0.5,
        # barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        # barglyphs=BarGlyphs("[=> ]"),
        color=:yellow,
        barlen=0,
        # enabled=false,
        )
	for i in 1:iters
		result += rand() # for example
		sleep(0.05)
		next!(p;
			# showvalues = [("iteration",i), (:result,result)]
			)
	end

	return result
end

include("../MCMC_fit.jl")

end # module JDRPM
