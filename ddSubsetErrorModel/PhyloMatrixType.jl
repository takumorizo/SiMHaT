Include("SortingCache2DType.jl")


module PhyloMatrixType
    using ..SortingCache2DType

    mutable struct PhyloMatrix{I}
        cache::SortingCache2D{I}
        Ly::Array{Bool, 1}  # left most min 1 consistent hold true or not.
        Lmin::Array{I, 1} # left most min 1, index is a private version
    end
    export PhyloMatrix

    function init(R::I, C::I;
                  descend::Bool = true, linkedValue::I = (I)(1)) where {I <: Integer}
        cache::SortingCache2D{I} = SortingCache2DType.init(R, C,
                                                       descend = descend, linkedValue = linkedValue)
        Ly::Array{Bool, 1}  = ones(Bool, C+(I)(2))
        Lmin::Array{I, 1}   = ones(I, C+(I)(2))

        return PhyloMatrix(cache, Ly, Lmin)
    end

    function add!(phylo::PhyloMatrix{I}, col::I, v::AbstractArray{I,1};
                  updateLy::Bool = true)::Nothing where {I <: Integer}
        privateCol::I = SortingCache2DType.at(phylo.cache, col)
        @assert privateCol ∉ phylo.cache.sortedCols && ((I)(1) <= col <= phylo.cache.C)
        addAt::I = SortingCache2DType._binSearch!(phylo.cache.matrix, phylo.cache.sortedCols,
                                               v, privateCol, phylo.cache.descend)
        LyUpdates::Array{I,1} =  SortingCache2DType._addLinkedList!(phylo.cache.matrix,
                                                                 phylo.cache.LLLeft, phylo.cache.LLRight,
                                                                 phylo.cache.sortedCols, addAt)
        if updateLy
            _updateL!(phylo.cache, phylo.Ly, phylo.Lmin, LyUpdates)
        end
        return nothing
        # println("add!")
        # println(phylo.cache.sortedCols)
        # println(phylo.cache.matrix[:, phylo.cache.sortedCols])
        # println(phylo.cache.LLLeft)
        # println(LyUpdates)
        # println(phylo.Ly)
        # println(phylo.Lmin)
    end

    function rm!(phylo::PhyloMatrix{I}, col::I; updateLy::Bool = true)::Nothing where {I <: Integer}
        privateCol::I = SortingCache2DType.at(phylo.cache, col)
        @assert privateCol ∈ phylo.cache.sortedCols && ((I)(1) <= col <= phylo.cache.C)
        # rmAt::I = findfirst(phylo.cache.sortedCols, privateCol)
        rmAt::I = something(findfirst(isequal(privateCol), phylo.cache.sortedCols), 0)
        LyUpdates::Array{I,1} = SortingCache2DType._rmLinkedList!(phylo.cache.matrix,
                                                               phylo.cache.LLLeft,
                                                               phylo.cache.LLRight,
                                                               phylo.cache.sortedCols, rmAt)
        deleteat!(phylo.cache.sortedCols, something(findfirst(isequal(privateCol), phylo.cache.sortedCols), 0))
        # deleteat!(phylo.cache.sortedCols, findfirst(phylo.cache.sortedCols, privateCol))
        if updateLy
            _updateL!(phylo.cache, phylo.Ly, phylo.Lmin, LyUpdates)
        end
        return nothing
    end

    function edit!(phylo::PhyloMatrix{I}, row::I, col::I, value::I)::Nothing where {I <: Integer}
        privateCol::I = SortingCache2DType.at(phylo.cache, col)
        isIn::Bool = (privateCol ∈ phylo.cache.sortedCols)
        if isIn && value == phylo.cache.matrix[row, privateCol]
            # println("ignore edit")
            return nothing
        else
            # println("start edit")
            v::Array{I,1} = deepcopy(phylo.cache.matrix[:, privateCol])
            v[row] = value
            (isIn) && rm!(phylo, col, updateLy = true)
            add!(phylo, col, v, updateLy = true)
            return nothing
        end

        # v::Array{I,1} = deepcopy(phylo.cache.matrix[:, privateCol])
        # v[row] = value
        # if privateCol ∈ phylo.cache.sortedCols
        #     rm!(phylo, col, updateLy = false)
        # end
        # add!(phylo, col, v, updateLy = false)
        # updateLy::Array{I,1} = [ privateCol ]
        # (value == phylo.cache.linkedValue) && push!(updateLy, phylo.cache.LLRight[row, privateCol])
        # _updateL!(phylo.cache, phylo.Ly, phylo.Lmin, updateLy)
        # return nothing
    end

    function isTree(phylo::PhyloMatrix{I})::Bool where {I <: Integer}
        return mapreduce(x->x, &, phylo.Ly)
    end

    function _updateL!(cache::SortingCache2D{I}, Ly::Array{Bool,1}, Lmin::Array{I,1}, LyUpdates::Array{I,1};
                     linkedValue::I = (I)(1))::Nothing where {I <: Integer}
        for c in LyUpdates
            Lmin[c] = cache.C + (I)(2)
            Ly[c] = true
            for r in (I)(1):cache.R
                if cache.matrix[r, c] == linkedValue
                    if Lmin[c] != (cache.C + (I)(2)) && Lmin[c] != cache.LLLeft[r,c]
                        Ly[c] = false
                        break
                    end
                    Lmin[c] = min(Lmin[c],cache.LLLeft[r,c])
                end
            end
        end
        return nothing
    end

end
