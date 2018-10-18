Include("phyloMatrix.jl")

module BuffPhyloMatrixType
    using ..PhyloMatrixType

    mutable struct BuffPhyloMatrix{I}
        phylo::PhyloMatrix{I}
        colSize::I
        bufferSize::I
        bufferUsed::I
        addedSetInBuffer::Set{I}
    end
    export BuffPhyloMatrix

    function init(R::I, C::I;
                  descend::Bool = true, linkedValue::I = (I)(1),
                  bufferSize::I = (I)(10)) where {I <: Integer}
        phylo::PhyloMatrix{I} = PhyloMatrixType.init(R, C+bufferSize, descend = descend, linkedValue = linkedValue)
        return BuffPhyloMatrix(phylo, C, bufferSize, (I)(0), Set{I}())
    end

    function novelIndex(buffPhylo::BuffPhyloMatrix{I})::I where {I <: Integer}
        return buffPhylo.colSize + buffPhylo.bufferUsed + (I)(1)
    end

    function add!(buffPhylo::BuffPhyloMatrix{I}, col::I, v::AbstractArray{I,1})::Nothing where {I <: Integer}
        @assert (I)(1) <= col <= buffPhylo.colSize
        PhyloMatrixType.add!(buffPhylo.phylo, col, v)
        return nothing
    end

    function rm!(buffPhylo::BuffPhyloMatrix{I}, col::I)::Nothing where {I <: Integer}
        @assert (I)(1) <= col <= buffPhylo.colSize
        PhyloMatrixType.rm!(buffPhylo.phylo, col)
        return nothing
    end

    function edit!(buffPhylo::BuffPhyloMatrix{I}, row::I, col::I, value::I)::Nothing where {I <: Integer}
        @assert (I)(1) <= col <= buffPhylo.colSize
        PhyloMatrixType.edit!(buffPhylo.phylo, row, col, value)
        return nothing
    end

    function isTree(buffPhylo::BuffPhyloMatrix{I})::Bool where {I <: Integer}
        return PhyloMatrixType.isTree(buffPhylo.phylo)
    end

    function push!(buffPhylo::BuffPhyloMatrix{I}, v::AbstractArray{I,1})::Nothing where {I <: Integer}
        @assert buffPhylo.bufferUsed < buffPhylo.bufferSize
        nextIndex::I = novelIndex(buffPhylo)
        Base.push!(buffPhylo.addedSetInBuffer, _binaryArrayToInt(v))
        buffPhylo.bufferUsed += (I)(1)
        PhyloMatrixType.add!(buffPhylo.phylo, nextIndex, v)
        return nothing
    end

    function pop!(buffPhylo::BuffPhyloMatrix{I})::Nothing where {I <: Integer}
        if buffPhylo.bufferUsed > 0
            rmIndex::I = novelIndex(buffPhylo) - (I)(1)
            Base.pop!(buffPhylo.addedSetInBuffer, _binaryArrayToInt(view(buffPhylo.phylo.cache.matrix, : , rmIndex + 1 )))
            buffPhylo.bufferUsed -= (I)(1)
            PhyloMatrixType.rm!(buffPhylo.phylo, rmIndex)
        end
        return nothing
    end

    function clearBuffer!(buffPhylo::BuffPhyloMatrix{I})::Nothing where {I <: Integer}
        while buffPhylo.bufferUsed > 0
            pop!(buffPhylo)
        end
        return nothing
    end

    function update!(buffPhylo::BuffPhyloMatrix{I},
                     B::Dict{Tuple{I,I}, I},
                     usageS::Dict{I,Array{I,1}},
                     usageV::Dict{I,Array{I,1}};
                     blocktype::I = (I)(1))::Nothing where {I <: Integer}
        R::I = buffPhylo.phylo.cache.R
        setNext::Set{I} = _summaryBits(B, usageS, usageV, R,
                                                        blocktype = blocktype)
        prevDiff::Set{I} = setdiff(buffPhylo.addedSetInBuffer, setNext)
        nextDiff::Set{I} = setdiff(setNext, buffPhylo.addedSetInBuffer)

        if length(prevDiff) > 0
            clearBuffer!(buffPhylo)
            for num in setNext
                push!(buffPhylo, _intToBinaryArray(num, R))
            end
        elseif length(nextDiff) > 0
            for num in nextDiff
                push!(buffPhylo, _intToBinaryArray(num, R))
            end
        elseif length(prevDiff) == 0 && length(nextDiff) == 0
            # do nothing
        end
        return nothing
    end

    function _binaryArrayToInt(v::AbstractArray{I,1}; ascend::Bool = true)::I where {I <: Integer}
        ans::I = (I)(0)
        for i in (I)(1):length(v)
            digit::I = (ascend) * (i-1) + (!ascend) * (length(v) - i)
            ans += (v[i]) * (I)(2) ^ (digit)
        end
        return ans
    end

    function _intToBinaryArray(num::I, vectorSize::I; ascend::Bool = true)::Array{I, 1} where {I <: Integer}
        v::Array{I, 1} = digits(typeof(num), num, base=(I)(2), pad=(Int)(vectorSize))
        # The reason for 'pad' to be typed as ::Int is the inner implementations of intfuncs.jl.
        # This should be fixed in the later versions of julia...
        if !ascend
            v = v[length(v):-1:1]
        end
        return v
    end

    function _summaryBits(B::Dict{Tuple{I,I}, I}, usageS::Dict{I,Array{I,1}},
                         usageV::Dict{I,Array{I,1}},
                         vectorSize::I;
                         blocktype::I = (I)(1))::Set{I} where {I <: Integer}
        ans::Set{I} = Set{I}()
        for m in keys(usageV)
            v::Array{I,1} = zeros(I, vectorSize)
            for c in keys(usageS)
                if B[(c,m)] == blocktype
                    v[usageS[c]] .= (I)(1)
                end
            end
            Base.push!(ans, _binaryArrayToInt(v))
        end
        return ans
    end

end
