let
    global Include
    existing_files = Set{AbstractString}([])
    function Include(path::AbstractString, guard::Bool = true)
        if guard
            if (path ∉ existing_files)
                print("include the julia file of: "); println(path)
                include(path)
                push!(existing_files, path)
            else
                print("avoid re-inclusion of the julia file of: "); println(path)
                println(existing_files)
            end
        else
            include(path)
            push!(existing_files, path)
        end
    end
end
