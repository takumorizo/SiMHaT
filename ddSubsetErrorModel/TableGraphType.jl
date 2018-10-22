
module TableGraphType
    export TableGraph
    mutable struct TableGraph{I, R}
        E::Array{I, 2} # edge[to, from] = 1
        V::I
        W::Array{R, 2}
    end

    function init(vertex_size::I, distance::Array{R, 2})::TableGraph{I, R} where {I <: Integer, R <: Real}
        @assert size(distance)[1] == vertex_size && size(distance)[2] == vertex_size
        g::TableGraph{I, R} = TableGraph{I,R}(zeros(I, vertex_size, vertex_size), vertex_size, distance)
        for n in (I)(1):vertex_size
            g.E[n, n] = (I)(1)
        end
        return g
    end

    function dfs!(g::TableGraph{I}, from::I, S::Set{I}, both = true)::Nothing where {I <: Integer}
        push!(S, from)
        for to in (I)(1):g.V
            ( g.E[to, from] == (I)(1) && to ∉ S ) && (dfs!(g, (I)(to), S, both))
            ( both && g.E[from, to] == (I)(1) && to ∉ S ) && (dfs!(g, (I)(to), S, both))
        end
        return nothing
    end

    function rmedge!(g::TableGraph{I}, from::I, to::I)::Nothing where {I <: Integer}
        g.E[to, from] = (I)(0)
        return nothing
    end

    function addedge!(g::TableGraph{I}, from::I, to::I)::Nothing where {I <: Integer}
        g.E[to, from] = (I)(1)
        return nothing
    end

    function getlink(g::TableGraph{I}, i::I) where {I <: Integer}
        # return findfirst(g.E[:, i], (I)(1))
        return something(findfirst(isequal((I)(1)), g.E[:, i]), 0)
    end


    function diffgroup!(g::TableGraph{I, R},
                        s_pp::Set{I}, s_pn::Set{I}, s_np::Set{I}, s_nn::Set{I},
                        from::I, to_prev::I, to_next::I; both = true)::Nothing where {I <: Integer, R <: Real}
        dfs!(g, to_prev, s_pp, both)
        dfs!(g, to_next, s_pn, both)

        if from ∈ s_pp && from ∈ s_pn
            while length(s_pn) > 0
                pop!(s_pn)
            end
        end

        rmedge!(g, from, to_prev)
        addedge!(g, from, to_next)

        dfs!(g, to_prev, s_np, both)
        dfs!(g, to_next, s_nn, both)

        if from ∈ s_np && from ∈ s_nn
            while length(s_nn) > 0
                pop!(s_nn)
            end
        end
        rmedge!(g, from, to_next)
        addedge!(g, from, to_prev)
        return nothing
    end
end
