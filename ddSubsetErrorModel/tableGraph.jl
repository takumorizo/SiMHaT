
module tableGraph
    export TableGraph
    mutable struct TableGraph{I, R}
        E::Array{I, 2} # edge[to, from] = 1
        V::I
        W::Array{R, 2}
    end

    function init(vertexSize::I, distance::Array{R, 2})::TableGraph{I, R} where {I <: Integer, R <: Real}
        @assert size(distance)[1] == vertexSize && size(distance)[2] == vertexSize
        g::TableGraph{I, R} = TableGraph{I,R}(zeros(I, vertexSize, vertexSize), vertexSize, distance)
        for n in (I)(1):vertexSize
            g.E[n, n] = (I)(1)
        end
        return g
    end

    function DFS!(g::TableGraph{I}, from::I, S::Set{I}, both = true)::Nothing where {I <: Integer}
        push!(S, from)
        for to in (I)(1):g.V
            ( g.E[to, from] == (I)(1) && to ∉ S ) && (DFS!(g, (I)(to), S, both))
            ( both && g.E[from, to] == (I)(1) && to ∉ S ) && (DFS!(g, (I)(to), S, both))
        end
        return nothing
    end

    function rmEdge!(g::TableGraph{I}, from::I, to::I)::Nothing where {I <: Integer}
        g.E[to, from] = (I)(0)
        return nothing
    end

    function addEdge!(g::TableGraph{I}, from::I, to::I)::Nothing where {I <: Integer}
        g.E[to, from] = (I)(1)
        return nothing
    end

    function getLink(g::TableGraph{I}, i::I) where {I <: Integer}
        # return findfirst(g.E[:, i], (I)(1))
        return something(findfirst(isequal((I)(1)), g.E[:, i]), 0)
    end


    function diffGroup!(g::TableGraph{I, R},
                        S_pp::Set{I}, S_pn::Set{I}, S_np::Set{I}, S_nn::Set{I},
                        from::I, to_prev::I, to_next::I; both = true)::Nothing where {I <: Integer, R <: Real}
        DFS!(g, to_prev, S_pp, both)
        DFS!(g, to_next, S_pn, both)

        if from ∈ S_pp && from ∈ S_pn
            while length(S_pn) > 0
                pop!(S_pn)
            end
        end

        rmEdge!(g, from, to_prev)
        addEdge!(g, from, to_next)

        DFS!(g, to_prev, S_np, both)
        DFS!(g, to_next, S_nn, both)

        if from ∈ S_np && from ∈ S_nn
            while length(S_nn) > 0
                pop!(S_nn)
            end
        end
        rmEdge!(g, from, to_next)
        addEdge!(g, from, to_prev)
        return nothing
    end
end
