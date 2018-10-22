
module SortingCache2DType
    mutable struct SortingCache2D{I}
        matrix::Array{I,2}  # R * (1+C+1+bufferSize)
        llleft::Array{I,2}  # R * (1+C+1+bufferSize)
        llright::Array{I,2} # R * (1+C+1+bufferSize)

        sorted_cols::Array{I,1} # col private indexing
        R::I # row size public indexing
        C::I # col size public indexing
        descend::Bool
        linked_value::I
    end
    export SortingCache2D

    @inline function at(cache::SortingCache2D{I}, i::I)::I where {I <: Integer}
        return (I)(1)+i
    end

    function init(R::I, C::I;
                  descend::Bool = true, linked_value::I = (I)(1)) where {I <: Integer}
        matrix::Array{I,2}  = zeros(R, C + (I)(2))
        llleft::Array{I,2}  = zeros(R, C + (I)(2))
        llright::Array{I,2} = zeros(R, C + (I)(2))

        for r in (I)(1):R
            llleft[r, C + (I)(2)] = (I)(1)
            llright[r, (I)(1)]    = C + (I)(2)
        end

        sorted_cols::Array{I, 1} = []
        bufferUsed::I = (I)(0)

        return SortingCache2D(
                matrix, llleft, llright,
                sorted_cols, R, C,
                descend, linked_value)
    end

    function add!(cache::SortingCache2D{I}, col::I, v::AbstractArray{I,1})::Nothing where {I <: Integer}
        privateCol::I = at(cache, col)
        @assert privateCol ∉ cache.sorted_cols && ((I)(1) <= col <= cache.C)
        addAt::I = _binsearch!(cache.matrix, cache.sorted_cols,
                                               v, privateCol, cache.descend)
        _addlinkedlist!(cache.matrix, cache.llleft, cache.llright,
                                        cache.sorted_cols, addAt)
        return nothing
    end

    function rm!(cache::SortingCache2D{I}, col::I)::Nothing where {I <: Integer}
        privateCol::I = at(cache, col)
        @assert privateCol ∈ cache.sorted_cols && ((I)(1) <= col <= cache.C)
        # rmAt::I = findfirst(cache.sorted_cols, privateCol)
        rmAt::I = something(findfirst(isequal(privateCol), cache.sorted_cols), 0)
        _rmlinkedlist!(cache.matrix, cache.llleft, cache.llright,
                                       cache.sorted_cols, rmAt)
        deleteat!(cache.sorted_cols, rmAt)
        return nothing
    end

    function edit!(cache::SortingCache2D{I}, row::I, col::I, value::I)::Nothing where {I <: Integer}
        @assert (I)(1) <= row <= cache.R
        privateCol::I = at(cache, col)
        if cache.matrix[row, privateCol] == value
            return
        end
        v::Array{I,1} = deepcopy(cache.matrix[:, privateCol])
        v[row] = value
        if privateCol ∈ cache.sorted_cols
            rm!(cache, col)
        end
        add!(cache, col, v)
        return nothing
    end

    function _eq(a::AbstractArray{I,1}, b::AbstractArray{I,1})::Bool where{I <: Integer}
        return mapreduce(x->x, &, a.==b)
    end

    # a > b
    function _gt(a::AbstractArray{I,1}, b::AbstractArray{I,1})::Bool where{I <: Integer}
        isGtEq::Bool = true
        isGt::Bool   = false
        for i in 1:length(a)
            isGtEq &= (a[i] >= b[i])
            isGt    = (a[i] >  b[i])
            if !(isGtEq)
                return false
            end
            if isGt
                return true
            end
        end
        return false
    end

    # a < b
    function _lt(a::AbstractArray{I,1}, b::AbstractArray{I,1})::Bool where{I <: Integer}
        isLtEq::Bool = true
        isLt::Bool   = false
        for i in 1:length(a)
            isLtEq &= (a[i] <= b[i])
            isLt    = (a[i] <  b[i])
            if !(isLtEq)
                return false
            end
            if isLt
                return true
            end
        end
        return false
    end

    # matrix: private indexing
    # sorted_cols: private indexing
    # v : added vector
    # col : private indexing
    function _binsearch_ascend!(matrix::AbstractArray{I, 2},
                                sorted_cols::AbstractArray{I,1},
                                v::AbstractArray{I, 1}, col::I)::I where {I <: Integer}
        R::I, C::I = size(matrix)
        matrix[:,col] = deepcopy(v)
        if length(sorted_cols) == 0
            push!(sorted_cols, col)
            return (I)(1)
        else
            left::I  = (I)(1)
            right::I = length(sorted_cols)
            if !(_gt(matrix[:,sorted_cols[right]], v))
                push!(sorted_cols, col)
                return length(sorted_cols)
            elseif _lt(v, matrix[:,sorted_cols[left]])
                insert!(sorted_cols, 1, col)
                return (I)(1)
            else
                while right - left > 1
                    m::I = div(left+right, (I)(2))
                    if _eq(matrix[:,sorted_cols[m]], v)
                        left  = m
                        right = left+(I)(1)
                        break
                    end
                    trimRight::Bool = _gt(matrix[:, sorted_cols[m]], v)
                    if trimRight
                        right = m
                    else
                        left = m
                    end
                end
                insert!(sorted_cols, left+1, col)
                return left+(I)(1)
            end
        end
    end

    # matrix: private indexing
    # sorted_cols: private indexing
    # v : added vector
    # col : private indexing
    function _binsearch_descend!(matrix::AbstractArray{I, 2},
                                 sorted_cols::AbstractArray{I,1},
                                 v::AbstractArray{I, 1}, col::I)::I where {I <: Integer}
        R::I, C::I = size(matrix)
        matrix[:,col] = deepcopy(v)
        if length(sorted_cols) == 0
            push!(sorted_cols, col)
            return (I)(1)
        else
            left::I  = (I)(1)
            right::I = length(sorted_cols)
            # [left, right ], not left >= v > right
            if !(_lt(matrix[:,sorted_cols[right]], v))
                push!(sorted_cols, col)
                return length(sorted_cols)
            elseif _gt(v, matrix[:,sorted_cols[left]])
                insert!(sorted_cols, 1, col)
                return (I)(1)
            else
                while right - left > 1
                    m::I = div(left+right, (I)(2))
                    if _eq(matrix[:,sorted_cols[m]], v)
                        left  = m
                        right = left+(I)(1)
                        break
                    end
                    trimRight::Bool = _lt(matrix[:, sorted_cols[m]], v)
                    if trimRight # m < v <= left
                        right = m
                    else # m >= v > right
                        left = m
                    end
                end
                insert!(sorted_cols, left+1, col)
                return left+(I)(1)
            end
        end
    end

    function _binsearch!(matrix::AbstractArray{I, 2},
                         sorted_cols::AbstractArray{I,1},
                         v::AbstractArray{I, 1}, col::I,
                         descend::Bool)::I where {I <: Integer}
        if descend
            return _binsearch_descend!(matrix, sorted_cols, v, col)
        else
            return _binsearch_ascend!(matrix, sorted_cols, v, col)
        end
    end

    # matrix: private indexing
    # llleft: private indexing
    # llright: private indexing
    # sorted_cols: private indexing
    # col_order: order in sorted_cols
    # link_val : value which is linked bi-directionally
    # interval : window size for neighborhood search
    # Ly_update_list : return private index array where Ly must be updated
    function _rmlinkedlist!(matrix::AbstractArray{I, 2},
                            llleft::AbstractArray{I, 2},
                            llright::AbstractArray{I, 2},
                            sorted_cols::AbstractArray{I, 1},
                            col_order::I;
                            link_val::I = (I)(1),
                            interval::I = (I)(5))::Array{I,1} where {I <: Integer}
        @assert size(matrix) == size(llleft) == size(llright)
        R::I, C::I = size(matrix)
        col::I = sorted_cols[col_order]
        Ly_update_list::Array{I,1} = [col]
        for r in (I)(1):R; if matrix[r,col] == link_val;
                leftMost::I  = llleft[r,col]
                rightMost::I = llright[r,col]
                # println((leftMost, rightMost))
                llright[r,leftMost] = rightMost
                llleft[r,rightMost] = leftMost
                llright[r,col] = llleft[r,col] = (I)(0)
                push!(Ly_update_list, rightMost)
        end;end
        return Ly_update_list
    end

    # matrix: private indexing
    # llleft: private indexing
    # llright: private indexing
    # sorted_cols: private indexing
    # col_order: order in sorted_cols
    # link_val : value which is linked bi-directionally
    # interval : window size for neighborhood search
    # Ly_update_list : return private index array where Ly must be updated
    function _addlinkedlist!(matrix::AbstractArray{I, 2},
                             llleft::AbstractArray{I, 2},
                             llright::AbstractArray{I, 2},
                             sorted_cols::AbstractArray{I, 1},
                             col_order::I;
                             link_val::I = (I)(1),
                             interval::I = (I)(5))::Array{I,1} where {I <: Integer}
        @assert size(matrix) == size(llleft) == size(llright)
        R::I, C::I  = size(matrix)
        listSize::I = length(sorted_cols)
        start::I, fin::I = (I)(1), C
        colAt::I    = sorted_cols[col_order]

        v::Array{I,1} = matrix[:, sorted_cols[col_order]]
        colToOrder::Dict{I, I} = Dict{I, I}()
        for i in (I)(1):length(sorted_cols)
            colToOrder[sorted_cols[i]] = i
        end
        colToOrder[start] = (I)(0); colToOrder[fin] = fin
        Ly_update_list::Array{I,1} = [colAt]
        for r in (I)(1):R
            if v[r] == link_val
                before::Array{I,1} = sorted_cols[ max((I)(1), col_order - interval):(col_order - (I)(1)) ]
                before = before[length(before):-1:1]
                after::Array{I,1}  = sorted_cols[ min(listSize + (I)(1), col_order + (I)(1)):min(listSize, col_order + interval) ]
                leftMost::I = start
                rightMost::I = fin
                if link_val ∈ matrix[r, before] # checkLeft
                    for i in before; if matrix[r, i] == link_val
                        leftMost = i; break
                    end;end
                    rightMost = llright[r,leftMost]
                elseif link_val ∈ matrix[r, after] # checkRight
                    for i in after; if matrix[r, i] == link_val
                        rightMost = i; break
                    end;end
                    leftMost = llleft[r,rightMost]
                else # parseFromLinkedList
                    if col_order < div(length(sorted_cols), (I)(2))
                        rightMost = start
                        while colToOrder[rightMost] <= col_order
                            rightMost = llright[r,rightMost]
                        end
                        leftMost = llleft[r,rightMost]
                    else
                        leftMost = fin
                        while colToOrder[leftMost] >= col_order
                            leftMost = llleft[r,leftMost]
                        end
                        rightMost = llright[r,leftMost]
                    end
                end
                llright[r,leftMost] = colAt;  llright[r,colAt]    = rightMost
                llleft[r,rightMost] = colAt;  llleft[r,colAt]     = leftMost
                push!(Ly_update_list, rightMost)
            else
                # llright[r,colAt]    = 0
                # llleft[r,colAt]     = 0
            end
        end
        return Ly_update_list
    end

end
