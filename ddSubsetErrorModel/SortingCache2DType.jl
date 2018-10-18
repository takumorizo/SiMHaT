Include("phylogenyTree.jl")

module SortingCache2DType
    mutable struct SortingCache2D{I}
        matrix::Array{I,2}  # R * (1+C+1+bufferSize)
        LLLeft::Array{I,2}  # R * (1+C+1+bufferSize)
        LLRight::Array{I,2} # R * (1+C+1+bufferSize)

        sortedCols::Array{I,1} # col private indexing
        R::I # row size public indexing
        C::I # col size public indexing
        descend::Bool
        linkedValue::I
    end
    export SortingCache2D

    @inline function at(cache::SortingCache2D{I}, i::I)::I where {I <: Integer}
        return (I)(1)+i
    end

    function init(R::I, C::I;
                  descend::Bool = true, linkedValue::I = (I)(1)) where {I <: Integer}
        matrix::Array{I,2}  = zeros(R, C + (I)(2))
        LLLeft::Array{I,2}  = zeros(R, C + (I)(2))
        LLRight::Array{I,2} = zeros(R, C + (I)(2))

        for r in (I)(1):R
            LLLeft[r, C + (I)(2)] = (I)(1)
            LLRight[r, (I)(1)]    = C + (I)(2)
        end

        sortedCols::Array{I, 1} = []
        bufferUsed::I = (I)(0)

        return SortingCache2D(
                matrix, LLLeft, LLRight,
                sortedCols, R, C,
                descend, linkedValue)
    end

    function add!(cache::SortingCache2D{I}, col::I, v::AbstractArray{I,1})::Nothing where {I <: Integer}
        privateCol::I = at(cache, col)
        @assert privateCol ∉ cache.sortedCols && ((I)(1) <= col <= cache.C)
        addAt::I = _binSearch!(cache.matrix, cache.sortedCols,
                                               v, privateCol, cache.descend)
        _addLinkedList!(cache.matrix, cache.LLLeft, cache.LLRight,
                                        cache.sortedCols, addAt)
        return nothing
    end

    function rm!(cache::SortingCache2D{I}, col::I)::Nothing where {I <: Integer}
        privateCol::I = at(cache, col)
        @assert privateCol ∈ cache.sortedCols && ((I)(1) <= col <= cache.C)
        # rmAt::I = findfirst(cache.sortedCols, privateCol)
        rmAt::I = something(findfirst(isequal(privateCol), cache.sortedCols), 0)
        _rmLinkedList!(cache.matrix, cache.LLLeft, cache.LLRight,
                                       cache.sortedCols, rmAt)
        deleteat!(cache.sortedCols, rmAt)
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
        if privateCol ∈ cache.sortedCols
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
    # sortedCols: private indexing
    # v : added vector
    # col : private indexing
    function _binSearchAscend!(matrix::AbstractArray{I, 2},
                              sortedCols::AbstractArray{I,1},
                              v::AbstractArray{I, 1}, col::I)::I where {I <: Integer}
        R::I, C::I = size(matrix)
        matrix[:,col] = deepcopy(v)
        if length(sortedCols) == 0
            push!(sortedCols, col)
            return (I)(1)
        else
            left::I  = (I)(1)
            right::I = length(sortedCols)
            if !(_gt(matrix[:,sortedCols[right]], v))
                push!(sortedCols, col)
                return length(sortedCols)
            elseif _lt(v, matrix[:,sortedCols[left]])
                insert!(sortedCols, 1, col)
                return (I)(1)
            else
                while right - left > 1
                    m::I = div(left+right, (I)(2))
                    if _eq(matrix[:,sortedCols[m]], v)
                        left  = m
                        right = left+(I)(1)
                        break
                    end
                    trimRight::Bool = _gt(matrix[:, sortedCols[m]], v)
                    if trimRight
                        right = m
                    else
                        left = m
                    end
                end
                insert!(sortedCols, left+1, col)
                return left+(I)(1)
            end
        end
    end

    # matrix: private indexing
    # sortedCols: private indexing
    # v : added vector
    # col : private indexing
    function _binSearchDescend!(matrix::AbstractArray{I, 2},
                               sortedCols::AbstractArray{I,1},
                               v::AbstractArray{I, 1}, col::I)::I where {I <: Integer}
        R::I, C::I = size(matrix)
        matrix[:,col] = deepcopy(v)
        if length(sortedCols) == 0
            push!(sortedCols, col)
            return (I)(1)
        else
            left::I  = (I)(1)
            right::I = length(sortedCols)
            # [left, right ], not left >= v > right
            if !(_lt(matrix[:,sortedCols[right]], v))
                push!(sortedCols, col)
                return length(sortedCols)
            elseif _gt(v, matrix[:,sortedCols[left]])
                insert!(sortedCols, 1, col)
                return (I)(1)
            else
                while right - left > 1
                    m::I = div(left+right, (I)(2))
                    if _eq(matrix[:,sortedCols[m]], v)
                        left  = m
                        right = left+(I)(1)
                        break
                    end
                    trimRight::Bool = _lt(matrix[:, sortedCols[m]], v)
                    if trimRight # m < v <= left
                        right = m
                    else # m >= v > right
                        left = m
                    end
                end
                insert!(sortedCols, left+1, col)
                return left+(I)(1)
            end
        end
    end

    function _binSearch!(matrix::AbstractArray{I, 2},
                        sortedCols::AbstractArray{I,1},
                        v::AbstractArray{I, 1}, col::I,
                        descend::Bool)::I where {I <: Integer}
        if descend
            return _binSearchDescend!(matrix, sortedCols, v, col)
        else
            return _binSearchAscend!(matrix, sortedCols, v, col)
        end
    end

    # matrix: private indexing
    # LLLeft: private indexing
    # LLRight: private indexing
    # sortedCols: private indexing
    # colOrder: order in sortedCols
    # linkVal : value which is linked bi-directionally
    # interval : window size for neighborhood search
    # LyUpdates : return private index array where Ly must be updated
    function _rmLinkedList!(matrix::AbstractArray{I, 2},
                          LLLeft::AbstractArray{I, 2},
                          LLRight::AbstractArray{I, 2},
                          sortedCols::AbstractArray{I, 1},
                          colOrder::I;
                          linkVal::I = (I)(1),
                          interval::I = (I)(5))::Array{I,1} where {I <: Integer}
        @assert size(matrix) == size(LLLeft) == size(LLRight)
        R::I, C::I = size(matrix)
        col::I = sortedCols[colOrder]
        LyUpdates::Array{I,1} = [col]
        for r in (I)(1):R; if matrix[r,col] == linkVal;
                leftMost::I  = LLLeft[r,col]
                rightMost::I = LLRight[r,col]
                # println((leftMost, rightMost))
                LLRight[r,leftMost] = rightMost
                LLLeft[r,rightMost] = leftMost
                LLRight[r,col] = LLLeft[r,col] = (I)(0)
                push!(LyUpdates, rightMost)
        end;end
        return LyUpdates
    end

    # matrix: private indexing
    # LLLeft: private indexing
    # LLRight: private indexing
    # sortedCols: private indexing
    # colOrder: order in sortedCols
    # linkVal : value which is linked bi-directionally
    # interval : window size for neighborhood search
    # LyUpdates : return private index array where Ly must be updated
    function _addLinkedList!(matrix::AbstractArray{I, 2},
                            LLLeft::AbstractArray{I, 2},
                            LLRight::AbstractArray{I, 2},
                            sortedCols::AbstractArray{I, 1},
                            colOrder::I;
                            linkVal::I = (I)(1),
                            interval::I = (I)(5))::Array{I,1} where {I <: Integer}
        @assert size(matrix) == size(LLLeft) == size(LLRight)
        R::I, C::I  = size(matrix)
        listSize::I = length(sortedCols)
        start::I, fin::I = (I)(1), C
        colAt::I    = sortedCols[colOrder]

        v::Array{I,1} = matrix[:, sortedCols[colOrder]]
        colToOrder::Dict{I, I} = Dict{I, I}()
        for i in (I)(1):length(sortedCols)
            colToOrder[sortedCols[i]] = i
        end
        colToOrder[start] = (I)(0); colToOrder[fin] = fin
        LyUpdates::Array{I,1} = [colAt]
        for r in (I)(1):R
            if v[r] == linkVal
                before::Array{I,1} = sortedCols[ max((I)(1), colOrder - interval):(colOrder - (I)(1)) ]
                before = before[length(before):-1:1]
                after::Array{I,1}  = sortedCols[ min(listSize + (I)(1), colOrder + (I)(1)):min(listSize, colOrder + interval) ]
                leftMost::I = start
                rightMost::I = fin
                if linkVal ∈ matrix[r, before] # checkLeft
                    for i in before; if matrix[r, i] == linkVal
                        leftMost = i; break
                    end;end
                    rightMost = LLRight[r,leftMost]
                elseif linkVal ∈ matrix[r, after] # checkRight
                    for i in after; if matrix[r, i] == linkVal
                        rightMost = i; break
                    end;end
                    leftMost = LLLeft[r,rightMost]
                else # parseFromLinkedList
                    if colOrder < div(length(sortedCols), (I)(2))
                        rightMost = start
                        while colToOrder[rightMost] <= colOrder
                            rightMost = LLRight[r,rightMost]
                        end
                        leftMost = LLLeft[r,rightMost]
                    else
                        leftMost = fin
                        while colToOrder[leftMost] >= colOrder
                            leftMost = LLLeft[r,leftMost]
                        end
                        rightMost = LLRight[r,leftMost]
                    end
                end
                LLRight[r,leftMost] = colAt;  LLRight[r,colAt]    = rightMost
                LLLeft[r,rightMost] = colAt;  LLLeft[r,colAt]     = leftMost
                push!(LyUpdates, rightMost)
            else
                # LLRight[r,colAt]    = 0
                # LLLeft[r,colAt]     = 0
            end
        end
        return LyUpdates
    end

end
