Include("phyloMatrix.jl")

module __buffPhyloMatrix
    function binaryArrayToInt(v::AbstractArray{I,1}; ascend::Bool = true)::I where {I <: Integer}
        ans::I = 0
        for i in 1:length(v)
            digit::I = (ascend) * (i-1) + (!ascend) * (length(v) - i)
            ans += (v[i]) * (2) ^ (digit)
        end
        return ans
    end

    function intToBinaryArray(num::I, vectorSize::I; ascend::Bool = true)::Array{I, 1} where {I <: Integer}
        v::Array{I, 1} = digits(num, 2, vectorSize)
        if !ascend
            v = v[length(v):-1:1]
        end
        return v
    end

    function summaryBits(B::Dict{Tuple{I,I}, I}, usageS::Dict{I,Array{I,1}},
                         usageV::Dict{I,Array{I,1}},
                         vectorSize::I;
                         blockType::I = 1)::Set{I} where {I <: Integer}
        ans::Set{I} = Set{I}()
        for m in keys(usageV)
            v::Array{I,1} = zeros(I, vectorSize)
            for c in keys(usageS)
                if B[(c,m)] == blockType
                    v[usageS[c]] .= 1
                end
            end
            push!(ans, binaryArrayToInt(v))
        end
        return ans
    end
end

module buffPhyloMatrix
    using phyloMatrix
    using __phyloMatrix
    using __buffPhyloMatrix

    type BuffPhyloMatrix{I}
        phylo::PhyloMatrix{I}
        colSize::I
        bufferSize::I
        bufferUsed::I
        addedSetInBuffer::Set{I}
    end
    export BuffPhyloMatrix

    function init(R::I, C::I;
                  descend::Bool = true, linkedValue::I = 1,
                  bufferSize::I = 10) where {I <: Integer}
        phylo::PhyloMatrix{I} = phyloMatrix.init(R, C+bufferSize, descend = descend, linkedValue = linkedValue)
        return BuffPhyloMatrix(phylo, C, bufferSize, 0, Set{I}())
    end

    function novelIndex(buffPhylo::BuffPhyloMatrix{I})::I where {I <: Integer}
        return buffPhylo.colSize + buffPhylo.bufferUsed + 1
    end

    function add!(buffPhylo::BuffPhyloMatrix{I}, col::I, v::AbstractArray{I,1})::Void where {I <: Integer}
        @assert 1 <= col <= buffPhylo.colSize
        phyloMatrix.add!(buffPhylo.phylo, col, v)
        return nothing
    end

    function rm!(buffPhylo::BuffPhyloMatrix{I}, col::I)::Void where {I <: Integer}
        @assert 1 <= col <= buffPhylo.colSize
        phyloMatrix.rm!(buffPhylo.phylo, col)
        return nothing
    end

    function edit!(buffPhylo::BuffPhyloMatrix{I}, row::I, col::I, value::I)::Void where {I <: Integer}
        @assert 1 <= col <= buffPhylo.colSize
        phyloMatrix.edit!(buffPhylo.phylo, row, col, value)
        return nothing
    end

    function isTree(buffPhylo::BuffPhyloMatrix{I})::Bool where {I <: Integer}
        return phyloMatrix.isTree(buffPhylo.phylo)
    end

    function push!(buffPhylo::BuffPhyloMatrix{I}, v::AbstractArray{I,1})::Void where {I <: Integer}
        @assert buffPhylo.bufferUsed < buffPhylo.bufferSize
        nextIndex::I = novelIndex(buffPhylo)
        Base.push!(buffPhylo.addedSetInBuffer, __buffPhyloMatrix.binaryArrayToInt(v))
        buffPhylo.bufferUsed += 1
        phyloMatrix.add!(buffPhylo.phylo, nextIndex, v)
        return nothing
    end

    function pop!(buffPhylo::BuffPhyloMatrix{I})::Void where {I <: Integer}
        if buffPhylo.bufferUsed > 0
            rmIndex::I = novelIndex(buffPhylo) - 1
            Base.pop!(buffPhylo.addedSetInBuffer, __buffPhyloMatrix.binaryArrayToInt(view(buffPhylo.phylo.cache.matrix, : , rmIndex + 1 )))
            buffPhylo.bufferUsed -= 1
            phyloMatrix.rm!(buffPhylo.phylo, rmIndex)
        end
        return nothing
    end

    function clearBuffer!(buffPhylo::BuffPhyloMatrix{I})::Void where {I <: Integer}
        while buffPhylo.bufferUsed > 0
            pop!(buffPhylo)
        end
        return nothing
    end

    function update!(buffPhylo::BuffPhyloMatrix{I},
                     B::Dict{Tuple{I,I}, I},
                     usageS::Dict{I,Array{I,1}},
                     usageV::Dict{I,Array{I,1}};
                     blockType::I = 1)::Void where {I <: Integer}
        R::I = buffPhylo.phylo.cache.R
        setNext::Set{I} = __buffPhyloMatrix.summaryBits(B, usageS, usageV, R,
                                                        blockType = blockType)
        prevDiff::Set{I} = setdiff(buffPhylo.addedSetInBuffer, setNext)
        nextDiff::Set{I} = setdiff(setNext, buffPhylo.addedSetInBuffer)

        if length(prevDiff) > 0
            clearBuffer!(buffPhylo)
            for num in setNext
                push!(buffPhylo, __buffPhyloMatrix.intToBinaryArray(num, R))
            end
        elseif length(nextDiff) > 0
            for num in nextDiff
                push!(buffPhylo, __buffPhyloMatrix.intToBinaryArray(num, R))
            end
        elseif length(prevDiff) == 0 && length(nextDiff) == 0
            # do nothing
        end
        return nothing
    end

end
