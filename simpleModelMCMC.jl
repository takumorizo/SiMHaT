


macro ran(s::Integer, e::Integer)
    (s+1):e
end

macro ran(len::Integer)
    1:len
end

println(@ran(10))

function main()
    inputfile = ARGS[1]
    open(inputfile) do f
        for line::String in eachline(f)
            line = strip(line)
            lineCol::Array{String,1} = split(line, '\t')

            outputStr::String = join(lineCol, '\t')
            println(outputStr)
        end
    end
end

main()
