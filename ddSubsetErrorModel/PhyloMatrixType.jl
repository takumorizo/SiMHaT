@Include "SortingCache2DType.jl"


module PhyloMatrixType
    using ..SortingCache2DType

    mutable struct PhyloMatrix{I}
        cache::SortingCache2D{I}
        Ly::Array{Bool, 1}  # left most min 1 consistent hold true or not.
        Lmin::Array{I, 1} # left most min 1, index is a private version
    end
    export PhyloMatrix

    function init(R::I, C::I;
                  descend::Bool = true, linked_value::I = (I)(1)) where {I <: Integer}
        cache::SortingCache2D{I} = SortingCache2DType.init(R, C,
                                                           descend = descend, linked_value = linked_value)
        Ly::Array{Bool, 1}  = ones(Bool, C+(I)(2))
        Lmin::Array{I, 1}   = ones(I, C+(I)(2))

        return PhyloMatrix(cache, Ly, Lmin)
    end

    function add!(phylo::PhyloMatrix{I}, col::I, v::AbstractArray{I,1};
                  update_Ly::Bool = true)::Nothing where {I <: Integer}
        privateCol::I = SortingCache2DType.at(phylo.cache, col)
        @assert privateCol ∉ phylo.cache.sorted_cols && ((I)(1) <= col <= phylo.cache.C)
        addAt::I = SortingCache2DType._binsearch!(phylo.cache.matrix, phylo.cache.sorted_cols,
                                                  v, privateCol, phylo.cache.descend)
        Ly_update_list::Array{I,1} =  SortingCache2DType._addlinkedlist!(phylo.cache.matrix,
                                                                         phylo.cache.llleft, phylo.cache.llright,
                                                                         phylo.cache.sorted_cols, addAt)
        if update_Ly
            _updatel!(phylo.cache, phylo.Ly, phylo.Lmin, Ly_update_list)
        end
        return nothing
        # println("add!")
        # println(phylo.cache.sorted_cols)
        # println(phylo.cache.matrix[:, phylo.cache.sorted_cols])
        # println(phylo.cache.llleft)
        # println(Ly_update_list)
        # println(phylo.Ly)
        # println(phylo.Lmin)
    end

    function rm!(phylo::PhyloMatrix{I}, col::I; update_Ly::Bool = true)::Nothing where {I <: Integer}
        privateCol::I = SortingCache2DType.at(phylo.cache, col)
        @assert privateCol ∈ phylo.cache.sorted_cols && ((I)(1) <= col <= phylo.cache.C)
        # rmAt::I = findfirst(phylo.cache.sorted_cols, privateCol)
        rmAt::I = something(findfirst(isequal(privateCol), phylo.cache.sorted_cols), 0)
        Ly_update_list::Array{I,1} = SortingCache2DType._rmlinkedlist!(phylo.cache.matrix,
                                                                       phylo.cache.llleft,
                                                                       phylo.cache.llright,
                                                                       phylo.cache.sorted_cols, rmAt)
        deleteat!(phylo.cache.sorted_cols, something(findfirst(isequal(privateCol), phylo.cache.sorted_cols), 0))
        # deleteat!(phylo.cache.sorted_cols, findfirst(phylo.cache.sorted_cols, privateCol))
        if update_Ly
            _updatel!(phylo.cache, phylo.Ly, phylo.Lmin, Ly_update_list)
        end
        return nothing
    end

    function edit!(phylo::PhyloMatrix{I}, row::I, col::I, value::I)::Nothing where {I <: Integer}
        privateCol::I = SortingCache2DType.at(phylo.cache, col)
        isIn::Bool = (privateCol ∈ phylo.cache.sorted_cols)
        if isIn && value == phylo.cache.matrix[row, privateCol]
            # println("ignore edit")
            return nothing
        else
            # println("start edit")
            v::Array{I,1} = deepcopy(phylo.cache.matrix[:, privateCol])
            v[row] = value
            (isIn) && rm!(phylo, col, update_Ly = true)
            add!(phylo, col, v, update_Ly = true)
            return nothing
        end

        # v::Array{I,1} = deepcopy(phylo.cache.matrix[:, privateCol])
        # v[row] = value
        # if privateCol ∈ phylo.cache.sorted_cols
        #     rm!(phylo, col, update_Ly = false)
        # end
        # add!(phylo, col, v, update_Ly = false)
        # update_Ly::Array{I,1} = [ privateCol ]
        # (value == phylo.cache.linked_value) && push!(update_Ly, phylo.cache.llright[row, privateCol])
        # _updatel!(phylo.cache, phylo.Ly, phylo.Lmin, update_Ly)
        # return nothing
    end

    function istree(phylo::PhyloMatrix{I})::Bool where {I <: Integer}
        return mapreduce(x->x, &, phylo.Ly)
    end

    function _updatel!(cache::SortingCache2D{I}, Ly::Array{Bool,1}, Lmin::Array{I,1}, Ly_update_list::Array{I,1};
                       linked_value::I = (I)(1))::Nothing where {I <: Integer}
        for c in Ly_update_list
            Lmin[c] = cache.C + (I)(2)
            Ly[c] = true
            for r in (I)(1):cache.R
                if cache.matrix[r, c] == linked_value
                    if Lmin[c] != (cache.C + (I)(2)) && Lmin[c] != cache.llleft[r,c]
                        Ly[c] = false
                        break
                    end
                    Lmin[c] = min(Lmin[c],cache.llleft[r,c])
                end
            end
        end
        return nothing
    end

end
