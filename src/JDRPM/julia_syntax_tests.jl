x = 1
println("(before the loop) x=$x")
for i in 1:3
	global x = 5
	y = x+i
	println("($i) y=$y x=$x")
end
println("(after the loop) x=$x")