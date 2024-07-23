module JDRPM

function my_example(iters::Number)
	result = 0.0
	for i in 1:iters
		print("iteration $i of $iters\r")
		# do my heavy work
		result += rand() # for example
		sleep(0.02)
	end
	return result
end
export my_example

end # module JDRPM
