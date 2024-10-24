title = "\t"^12 # space for title sections in log prints

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

macro showd(exs...)
	args = []
	for ex in exs
		push!(args, sprint(Base.show_unquoted, ex), " = ", esc(ex), '\n')
	end
	return :(string($(args...)))
end
