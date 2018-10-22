module PhylogenyTree
    function _isphylogenic(sorted_mat::AbstractArray{I, 2},
                           buffer2d::AbstractArray{I, 2};
                           offset::I = (I)(0), see_driver::Bool = false)::Bool where {I <: Integer}
        X, Y = size(sorted_mat)
        if (X,Y) != size(buffer2d)
            error("inconstent buffer2d size @__phylogenyTree.isphylogenic")
        end
        fill!(buffer2d, (I)(0))
        for x in (I)(1):X
            latest = (I)(0)
            for y in (I)(1):Y
                if sorted_mat[x,y] == ((I)(1) + offset)
                    buffer2d[x,y] = latest
                    latest = y
                end
            end
        end
        zeroCount::I = 0
        for y in (I)(1):Y
            Ly = maximum(buffer2d[:,y])
            zeroCount += (I)( (Ly == 0) && ( sum(sorted_mat[:,y]) > X * offset ) )
            for x in (I)(1):X
                if sorted_mat[x,y] == ((I)(1) + offset) && buffer2d[x,y] != Ly
                    return false
                end
            end
        end
        # if see_driver; return zeroCount == 1;
        # else; return true; end
        return (!see_driver) || (zeroCount==1)
    end

    function _radixsort(matrix::AbstractArray{I, 2};
                        ascend_digit::Bool = true,
                        base::I = (I)(2),
                        descend::Bool = false,
                        offset::I = (I)(0))::Array{I,1} where {I <: Integer}
        X::I, Y::I = size(matrix)
        bucket::Array{I, 2}     = zeros(base, Y)
        bucketSize::Array{I, 1} = zeros(base)
        ans::Array{I, 1}        = collect((I)(1):Y)
        digitOrder = (I)(1):X
        if !ascend_digit; digitOrder = X:(I)(-1):(I)(1); end

        for x in digitOrder
            # create bucket @ x pposition
            for y in ans
                bucketSize[ (matrix[x,y]-offset) + (I)(1) ] += (I)(1)
                idx = bucketSize[ (matrix[x,y]-offset) + (I)(1) ]
                bucket[(matrix[x,y]-offset)+(I)(1), idx] = y
            end
            # bucket to ans list conversion
            updateIdx = (I)(1)
            for b in (I)(1):base
                for elem in (I)(1):bucketSize[b]
                    ans[updateIdx] = bucket[b,elem]
                    updateIdx += (I)(1)
                end
            end
            for b in (I)(1):base
                bucketSize[b] = 0
                for y in (I)(1):Y; bucket[b,y] = (I)(0); end
            end
        end
        if descend
            return reverse(ans)
        else
            return ans
        end
    end

    function isphylogenic(binMat::AbstractArray{I, 2};
                          offset::I = (I)(0),
                          see_driver::Bool = false,
                          min_x::I = (I)(1), min_y::I = (I)(1))::Bool where {I <: Integer}
            # if min_x == 1 || min_y == 1
                # error("refactoring not done")
            # end
            X,Y = size(binMat)
            if X < min_x || Y < min_y
                return false
            end
            orderY::Array{I, 1} = _radixsort(binMat,ascend_digit = true, base = (I)(2), descend = true, offset = offset)
            buffer2d = fill((I)(0), X, Y)
            return _isphylogenic( view((binMat), (I)(1):X, orderY), buffer2d, offset = offset, see_driver = see_driver)
    end
end
