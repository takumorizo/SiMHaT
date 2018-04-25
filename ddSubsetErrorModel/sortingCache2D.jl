Include("phylogenyTree.jl")

module __sortingCache2D
    # a == b
    function eq(a::AbstractArray{I,1}, b::AbstractArray{I,1})::Bool where{I <: Integer}
        return mapreduce(x->x, &, a.==b)
    end

    # a > b
    function gt(a::AbstractArray{I,1}, b::AbstractArray{I,1})::Bool where{I <: Integer}
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
    function lt(a::AbstractArray{I,1}, b::AbstractArray{I,1})::Bool where{I <: Integer}
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
    function binSearchAscend!(matrix::AbstractArray{I, 2},
                              sortedCols::AbstractArray{I,1},
                              v::AbstractArray{I, 1}, col::I)::I where {I <: Integer}
        R::I, C::I = size(matrix)
        matrix[:,col] = deepcopy(v)
        if length(sortedCols) == 0
            push!(sortedCols, col)
            return 1
        else
            left::I  = 1
            right::I = length(sortedCols)
            if !(__sortingCache2D.gt(matrix[:,sortedCols[right]], v))
                push!(sortedCols, col)
                return length(sortedCols)
            elseif __sortingCache2D.lt(v, matrix[:,sortedCols[left]])
                insert!(sortedCols, 1, col)
                return 1
            else
                while right - left > 1
                    m::I = div(left+right, 2)
                    if eq(matrix[:,sortedCols[m]], v)
                        left  = m
                        right = left+1
                        break
                    end
                    trimRight::Bool = gt(matrix[:, sortedCols[m]], v)
                    if trimRight
                        right = m
                    else
                        left = m
                    end
                end
                insert!(sortedCols, left+1, col)
                return left+1
            end
        end
    end

    # matrix: private indexing
    # sortedCols: private indexing
    # v : added vector
    # col : private indexing
    function binSearchDescend!(matrix::AbstractArray{I, 2},
                               sortedCols::AbstractArray{I,1},
                               v::AbstractArray{I, 1}, col::I)::I where {I <: Integer}
        R::I, C::I = size(matrix)
        matrix[:,col] = deepcopy(v)
        if length(sortedCols) == 0
            push!(sortedCols, col)
            return 1
        else
            left::I  = 1
            right::I = length(sortedCols)
            # [left, right ], not left >= v > right
            if !(__sortingCache2D.lt(matrix[:,sortedCols[right]], v))
                push!(sortedCols, col)
                return length(sortedCols)
            elseif __sortingCache2D.gt(v, matrix[:,sortedCols[left]])
                insert!(sortedCols, 1, col)
                return 1
            else
                while right - left > 1
                    m::I = div(left+right, 2)
                    if eq(matrix[:,sortedCols[m]], v)
                        left  = m
                        right = left+1
                        break
                    end
                    trimRight::Bool = lt(matrix[:, sortedCols[m]], v)
                    if trimRight # m < v <= left
                        right = m
                    else # m >= v > right
                        left = m
                    end
                end
                insert!(sortedCols, left+1, col)
                return left+1
            end
        end
    end

    function binSearch!(matrix::AbstractArray{I, 2},
                        sortedCols::AbstractArray{I,1},
                        v::AbstractArray{I, 1}, col::I,
                        descend::Bool)::I where {I <: Integer}
        if descend
            return binSearchDescend!(matrix, sortedCols, v, col)
        else
            return binSearchAscend!(matrix, sortedCols, v, col)
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
    function rmLinkedList!(matrix::AbstractArray{I, 2},
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
        for r in 1:R; if matrix[r,col] == linkVal;
                leftMost::I  = LLLeft[r,col]
                rightMost::I = LLRight[r,col]
                # println((leftMost, rightMost))
                LLRight[r,leftMost] = rightMost
                LLLeft[r,rightMost] = leftMost
                LLRight[r,col] = LLLeft[r,col] = 0
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
    function addLinkedList!(matrix::AbstractArray{I, 2},
                            LLLeft::AbstractArray{I, 2},
                            LLRight::AbstractArray{I, 2},
                            sortedCols::AbstractArray{I, 1},
                            colOrder::I;
                            linkVal::I = (I)(1),
                            interval::I = (I)(5))::Array{I,1} where {I <: Integer}
        @assert size(matrix) == size(LLLeft) == size(LLRight)
        R::I, C::I  = size(matrix)
        listSize::I = length(sortedCols)
        start::I, fin::I = 1, C
        colAt::I    = sortedCols[colOrder]

        v::Array{I,1} = matrix[:, sortedCols[colOrder]]
        colToOrder::Dict{I, I} = Dict{I, I}()
        for i in 1:length(sortedCols)
            colToOrder[sortedCols[i]] = i
        end
        colToOrder[start] = 0; colToOrder[fin] = fin
        LyUpdates::Array{I,1} = [colAt]
        for r in 1:R
            if v[r] == linkVal
                before::Array{I,1} = sortedCols[ max(1, colOrder - interval):(colOrder - 1) ]
                before = before[length(before):-1:1]
                after::Array{I,1}  = sortedCols[ min(listSize + 1, colOrder + 1):min(listSize, colOrder + interval) ]
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
                    if colOrder < div(length(sortedCols), 2)
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


module sortingCache2D
    using __sortingCache2D

    type SortingCache2D{I}
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
        return 1+i
    end

    function init(R::I, C::I;
                  descend::Bool = true, linkedValue::I = (I)(1)) where {I <: Integer}
        matrix::Array{I,2}  = zeros(R, C + 2)
        LLLeft::Array{I,2}  = zeros(R, C + 2)
        LLRight::Array{I,2} = zeros(R, C + 2)

        for r in 1:R
            LLLeft[r, C + 2] = 1
            LLRight[r, 1]    = C + 2
        end

        sortedCols::Array{I, 1} = []
        bufferUsed::I = 0

        return SortingCache2D(
                matrix, LLLeft, LLRight,
                sortedCols, R, C,
                descend, linkedValue)
    end

    function add!(cache::SortingCache2D{I}, col::I, v::AbstractArray{I,1})::Void where {I <: Integer}
        privateCol::I = at(cache, col)
        @assert privateCol ∉ cache.sortedCols && (1 <= col <= cache.C)
        addAt::I = __sortingCache2D.binSearch!(cache.matrix, cache.sortedCols,
                                               v, privateCol, cache.descend)
        __sortingCache2D.addLinkedList!(cache.matrix, cache.LLLeft, cache.LLRight,
                                        cache.sortedCols, addAt)
        return nothing
    end

    function rm!(cache::SortingCache2D{I}, col::I)::Void where {I <: Integer}
        privateCol::I = at(cache, col)
        @assert privateCol ∈ cache.sortedCols && (1 <= col <= cache.C)
        rmAt::I = findfirst(cache.sortedCols, privateCol)
        __sortingCache2D.rmLinkedList!(cache.matrix, cache.LLLeft, cache.LLRight,
                                       cache.sortedCols, rmAt)
        deleteat!(cache.sortedCols, rmAt)
        return nothing
    end

    function edit!(cache::SortingCache2D{I}, row::I, col::I, value::I)::Void where {I <: Integer}
        @assert 1 <= row <= cache.R
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

end
