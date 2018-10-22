@Include "PhyloMatrixType.jl"

module BuffPhyloMatrixType
    using ..PhyloMatrixType

    mutable struct BuffPhyloMatrix{I}
        phylo::PhyloMatrix{I}
        col_size::I
        buffer_size::I
        buffer_used::I
        added_set_in_buffer::Set{I}
    end
    export BuffPhyloMatrix

    function init(R::I, C::I;
                  descend::Bool = true, linked_value::I = (I)(1),
                  buffer_size::I = (I)(10)) where {I <: Integer}
        phylo::PhyloMatrix{I} = PhyloMatrixType.init(R, C+buffer_size, descend = descend, linked_value = linked_value)
        return BuffPhyloMatrix(phylo, C, buffer_size, (I)(0), Set{I}())
    end

    function novel_index(buff_phylo::BuffPhyloMatrix{I})::I where {I <: Integer}
        return buff_phylo.col_size + buff_phylo.buffer_used + (I)(1)
    end

    function add!(buff_phylo::BuffPhyloMatrix{I}, col::I, v::AbstractArray{I,1})::Nothing where {I <: Integer}
        @assert (I)(1) <= col <= buff_phylo.col_size
        PhyloMatrixType.add!(buff_phylo.phylo, col, v)
        return nothing
    end

    function rm!(buff_phylo::BuffPhyloMatrix{I}, col::I)::Nothing where {I <: Integer}
        @assert (I)(1) <= col <= buff_phylo.col_size
        PhyloMatrixType.rm!(buff_phylo.phylo, col)
        return nothing
    end

    function edit!(buff_phylo::BuffPhyloMatrix{I}, row::I, col::I, value::I)::Nothing where {I <: Integer}
        @assert (I)(1) <= col <= buff_phylo.col_size
        PhyloMatrixType.edit!(buff_phylo.phylo, row, col, value)
        return nothing
    end

    function istree(buff_phylo::BuffPhyloMatrix{I})::Bool where {I <: Integer}
        return PhyloMatrixType.istree(buff_phylo.phylo)
    end

    function push!(buff_phylo::BuffPhyloMatrix{I}, v::AbstractArray{I,1})::Nothing where {I <: Integer}
        @assert buff_phylo.buffer_used < buff_phylo.buffer_size
        nextIndex::I = novel_index(buff_phylo)
        Base.push!(buff_phylo.added_set_in_buffer, _binaryarray_to_int(v))
        buff_phylo.buffer_used += (I)(1)
        PhyloMatrixType.add!(buff_phylo.phylo, nextIndex, v)
        return nothing
    end

    function pop!(buff_phylo::BuffPhyloMatrix{I})::Nothing where {I <: Integer}
        if buff_phylo.buffer_used > 0
            rmIndex::I = novel_index(buff_phylo) - (I)(1)
            Base.pop!(buff_phylo.added_set_in_buffer, _binaryarray_to_int(view(buff_phylo.phylo.cache.matrix, : , rmIndex + 1 )))
            buff_phylo.buffer_used -= (I)(1)
            PhyloMatrixType.rm!(buff_phylo.phylo, rmIndex)
        end
        return nothing
    end

    function clearbuffer!(buff_phylo::BuffPhyloMatrix{I})::Nothing where {I <: Integer}
        while buff_phylo.buffer_used > 0
            pop!(buff_phylo)
        end
        return nothing
    end

    function update!(buff_phylo::BuffPhyloMatrix{I},
                     B::Dict{Tuple{I,I}, I},
                     usage_s::Dict{I,Array{I,1}},
                     usage_v::Dict{I,Array{I,1}};
                     blocktype::I = (I)(1))::Nothing where {I <: Integer}
        R::I = buff_phylo.phylo.cache.R
        setNext::Set{I} = _summarybits(B, usage_s, usage_v, R,
                                                        blocktype = blocktype)
        prevDiff::Set{I} = setdiff(buff_phylo.added_set_in_buffer, setNext)
        nextDiff::Set{I} = setdiff(setNext, buff_phylo.added_set_in_buffer)

        if length(prevDiff) > 0
            clearbuffer!(buff_phylo)
            for num in setNext
                push!(buff_phylo, _int_to_binaryarray(num, R))
            end
        elseif length(nextDiff) > 0
            for num in nextDiff
                push!(buff_phylo, _int_to_binaryarray(num, R))
            end
        elseif length(prevDiff) == 0 && length(nextDiff) == 0
            # do nothing
        end
        return nothing
    end

    function _binaryarray_to_int(v::AbstractArray{I,1}; ascend::Bool = true)::I where {I <: Integer}
        ans::I = (I)(0)
        for i in (I)(1):length(v)
            digit::I = (ascend) * (i-1) + (!ascend) * (length(v) - i)
            ans += (v[i]) * (I)(2) ^ (digit)
        end
        return ans
    end

    function _int_to_binaryarray(num::I, vector_size::I; ascend::Bool = true)::Array{I, 1} where {I <: Integer}
        v::Array{I, 1} = digits(typeof(num), num, base=(I)(2), pad=(Int)(vector_size))
        # The reason for 'pad' to be typed as ::Int is the inner implementations of intfuncs.jl.
        # This should be fixed in the later versions of julia...
        if !ascend
            v = v[length(v):-1:1]
        end
        return v
    end

    function _summarybits(B::Dict{Tuple{I,I}, I}, usage_s::Dict{I,Array{I,1}},
                          usage_v::Dict{I,Array{I,1}},
                          vector_size::I;
                          blocktype::I = (I)(1))::Set{I} where {I <: Integer}
        ans::Set{I} = Set{I}()
        for m in keys(usage_v)
            v::Array{I,1} = zeros(I, vector_size)
            for c in keys(usage_s)
                if B[(c,m)] == blocktype
                    v[usage_s[c]] .= (I)(1)
                end
            end
            Base.push!(ans, _binaryarray_to_int(v))
        end
        return ans
    end

end
