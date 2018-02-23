module __phylogenyTree
    # assume that sorted Mat is coloum-wise-sorted in an descending order.
    function isPhylogenic(sortedMat::AbstractArray{I, 2},
                          buffer2D::AbstractArray{I, 2};
                          offset::I = 0, checkDriver::Bool = false)::Bool where {I <: Integer}
        X, Y = size(sortedMat)
        if (X,Y) != size(buffer2D)
            error("inconstent buffer2D size @__phylogenyTree.isPhylogenic")
        end
        fill!(buffer2D, (I)(0))
        for x in 1:X
            latest = 0
            for y in 1:Y
                if sortedMat[x,y] == (1 + offset)
                    buffer2D[x,y] = latest
                    latest = y
                end
            end
        end

        zeroCount::I = 0
        for y in 1:Y
            Ly = maximum(buffer2D[:,y])
            zeroCount += (I)( (Ly == 0) && ( sum(sortedMat[:,y]) > X * offset ) )
            for x in 1:X
                if buffer2D[x,y] > 0 && buffer2D[x,y] != Ly
                    return false
                end
            end
        end
        # if checkDriver; return zeroCount == 1;
        # else; return true; end
        return (!checkDriver) || (zeroCount==1)
    end

    function radixSort(matrix::AbstractArray{I, 2};
                       ascendDigit::Bool = true,
                       base::I = 2,
                       descend::Bool = false,
                       offset::I = 0)::Array{I,1} where {I <: Integer}
        X,Y = size(matrix)
        bucket::Array{I, 2}     = zeros(base, Y)
        bucketSize::Array{I, 1} = zeros(base)
        ans::Array{I, 1}        = collect(1:Y)
        digitOrder = 1:X
        if !ascendDigit; digitOrder = X:-1:1; end

        for x in digitOrder
            # create bucket @ x pposition
            for y in ans
                bucketSize[ (matrix[x,y]-offset) + 1 ] += 1
                idx = bucketSize[ (matrix[x,y]-offset) + 1 ]
                bucket[(matrix[x,y]-offset)+1, idx] = y
            end
            # bucket to ans list conversion
            updateIdx = 1
            for b in 1:base
                for elem in 1:bucketSize[b]
                    ans[updateIdx] = bucket[b,elem]
                    updateIdx += 1
                end
            end
            for b in 1:base
                bucketSize[b] = 0
                for y in 1:Y; bucket[b,y] = 0; end
            end
        end
        if descend
            return reverse(ans)
        else
            return ans
        end
    end
end


module phylogenyTree
    using __phylogenyTree
    function isPhylogenic(binMat::AbstractArray{I, 2};
                          offset::I = 0,
                          checkDriver::Bool = false,
                          minX::I = 1, minY::I = 1)::Bool where {I <: Integer}
            if minX == 1 || minY == 1
                error("refactoring not done")
            end
            X,Y = size(binMat)
            if X < minX || Y < minY
                return false
            end
            orderY::Array{I, 1} = __phylogenyTree.radixSort(binMat,
                                                            ascendDigit = true, base = (I)(2), descend = true, offset = offset)
            buffer2D = fill((I)(0), X, Y)
            return __phylogenyTree.isPhylogenic( view((binMat), 1:X, orderY), buffer2D, offset = offset, checkDriver = checkDriver)
    end
end
