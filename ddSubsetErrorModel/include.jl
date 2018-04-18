let
    global Include
    existingFiles = Set{AbstractString}([])
    function Include(path::AbstractString, guard::Bool = true)
        if guard
            if (path ∉ existingFiles)
                print("include the julia file of: "); println(path)
                include(path)
                push!(existingFiles, path)
            else
                print("avoid re-inclusion of the julia file of: "); println(path)
                println(existingFiles)
            end
        else
            include(path)
            push!(existingFiles, path)
        end
    end
end
