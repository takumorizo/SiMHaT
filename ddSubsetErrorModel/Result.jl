@Include "BuffPhyloMatrixType.jl"
@Include "SamplerType.jl"
@Include "InputParser.jl"

module Result
    using Plots
    using ..SamplerType
    using JLD2, FileIO
    using ..InputParser

    function writejld(smax, lnprobs, jldpath)
        # jldopen(jldpath, "w") do file
        #     addrequire(file, sampler)
        #     write(file, "smax", smax)
        #     write(file, "lnprobs", lnprobs)
        # end
        # save(jldpath, "smax", smax, "lnprobs", lnprobs)
        jldopen(jldpath, "w") do file
            file["smax"]    = smax
            file["lnprobs"] = lnprobs
        end
    end

    function readjld(jldpath)
        # smax = jldopen(jldpath, "r") do file
        #     read(file, "smax")
        # end
        # lnprobs = jldopen(jldpath, "r") do file
        #     read(file, "lnprobs")
        # end
        # return smax, lnprobs
        smax, lnprobs = load(jldpath, "smax", "lnprobs")
        println(typeof(smax))
        return smax, lnprobs
    end

    #= TODO:
    smax.usage_s -> smax.usage_s^(1)
    =#
    function viewmap(smax, lnprobs, output_dir::String)
        pyplot()
        iters = collect(1:length(lnprobs))
        plot(iters, lnprobs, size=(2000,2000))
        savefig(output_dir * "/lnprobs.png")

        acf = _acf(lnprobs)
        bar(acf)
        savefig(output_dir * "/acf.png")

        println("isValid state")
        println(SamplerType.isvalid(smax, debug = true))
        mat, rbreak, sort_r, cbreak, sort_c = _sort_biclustermatrix(smax.Z, smax.usage_s, smax.usage_v, smax.B)

        # println("MAP, B")
        # println(smax.B)
        # println("isError valid matrix")
        # println(BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(buffPhyloMatrix.BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix{Int64})
        # isTree(smax.treeCache)

        _view_matrix_heatmap(mat, rbreak, sort_r, cbreak, sort_c, "MAP.Z", file_path = output_dir * "/MAP_Z.png")
        _view_matrix_heatmap(smax.ln_p_data[sort_r,sort_c,1], rbreak, sort_r, cbreak, sort_c, "err", file_path = output_dir * "/err.png")
        _view_matrix_heatmap(smax.ln_p_data[sort_r,sort_c,2], rbreak, sort_r, cbreak, sort_c, "mat", file_path = output_dir * "/mat.png")
        _view_matrix_heatmap(smax.ln_p_data[sort_r,sort_c,3], rbreak, sort_r, cbreak, sort_c, "pat", file_path = output_dir * "/pat.png")

        _view_matrix_heatmap(smax.ln_p_data[:,:,1], rbreak, sort_r, cbreak, sort_c, "err", file_path = output_dir * "/err.nonSort.png")
        _view_matrix_heatmap(smax.ln_p_data[:,:,2], rbreak, sort_r, cbreak, sort_c, "mat", file_path = output_dir * "/mat.nonSort.png")
        _view_matrix_heatmap(smax.ln_p_data[:,:,3], rbreak, sort_r, cbreak, sort_c, "pat", file_path = output_dir * "/pat.nonSort.png")
        _view_matrix_heatmap(smax.H[sort_r,sort_c], rbreak, sort_r, cbreak, sort_c, "MAP.H", file_path = output_dir * "/MAP_H.png", clim = [1,2])
        _view_f_heatmap(smax.f, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.f", file_path = output_dir * "/MAP_f.png")
        _view_u_heatmap(smax.u, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.u", file_path = output_dir * "/MAP_u.png")
        _view_a_heatmap(smax.a, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.a", file_path = output_dir * "/MAP_a.png")
        _view_blocktype_heatmap(smax.B, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.B", file_path = output_dir * "/MAP_B.png")
    end

    function viewdata(sampled, output_dir::String)
        pyplot()

        cutoff = 200
        lnprobs = []; iters   = []; for (s,c) in sampled; push!(lnprobs, s.lnProb); push!(iters, c); end; plot(iters[cutoff:length(iters)], lnprobs[cutoff:length(iters)], size=(2000,2000))
        savefig(output_dir * "/lnprobs.png")

        smax = _get_mapstate(sampled)
        println("isValid state")
        println(SamplerType.isvalid(smax, debug = true))
        mat, rbreak, sort_r, cbreak, sort_c = _sort_biclustermatrix(smax.Z, smax.usage_s, smax.usage_v)

        # println("MAP, B")
        # println(smax.B)
        # println("isError valid matrix")
        # println(BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(buffPhyloMatrix.BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix{Int64})
        # isTree(smax.treeCache)

        _view_matrix_heatmap(mat, rbreak, sort_r, cbreak, sort_c, "MAP.Z", file_path = output_dir * "/MAP_Z.png")
        _view_matrix_heatmap(smax.ln_p_data[sort_r,sort_c,1], rbreak, sort_r, cbreak, sort_c, "err", file_path = output_dir * "/err.png")
        _view_matrix_heatmap(smax.ln_p_data[sort_r,sort_c,2], rbreak, sort_r, cbreak, sort_c, "mat", file_path = output_dir * "/mat.png")
        _view_matrix_heatmap(smax.ln_p_data[sort_r,sort_c,3], rbreak, sort_r, cbreak, sort_c, "pat", file_path = output_dir * "/pat.png")

        _view_matrix_heatmap(smax.H, rbreak, sort_r, cbreak, sort_c, "MAP.H", file_path = output_dir * "/MAP_H.png", clim = [1,2])
        _view_f_heatmap(smax.f, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.f", file_path = output_dir * "/MAP_f.png")
        _view_u_heatmap(smax.u, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.u", file_path = output_dir * "/MAP_u.png")
        _view_a_heatmap(smax.a, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.a", file_path = output_dir * "/MAP_a.png")
        _view_blocktype_heatmap(smax.B, smax.usage_s, rbreak, sort_r, smax.usage_v, cbreak, sort_c, "MAP.B", file_path = output_dir * "/MAP_B.png")
    end

    function evaluatefvalue(summary_path::String, score_path::String, thres::REAL = (REAL)(0.0), rm_error::Bool = true,
                            INT::Type{<:Integer} = Int32)::Tuple{REAL,REAL,REAL} where {REAL <: Real}
        if endswith(score_path, "jld2")
            sampled::Sampler{INT, REAL}, lnprobs::Array{REAL, 1} = readjld(score_path)
            return _evaluatefvalue(summary_path, sampled, thres, rm_error)
        else
            score::Array{REAL, 2}       = InputParser.parse_input_summary(score_path)
            return _evaluatefvalue(summary_path, score, thres, INT)
        end
    end

    function savefvalue(output_path::String, prescision::REAL, recall::REAL, fvalue::REAL, tag::String) where {REAL <: Real}
        @assert endswith(output_path, "jld2")
        if isfile(output_path)
            println(load(output_path))
            precisions::Dict{String, REAL} = load(output_path, "precisions")
            recalls::Dict{String, REAL}    = load(output_path, "recalls")
            fvalues::Dict{String, REAL}    = load(output_path, "fvalues")
            precisions[tag] = prescision
            recalls[tag]    = recall
            fvalues[tag]    = fvalue
            # precisions::Dict{String, REAL} = jldopen(output_path, "r") do file
            #     read(file, "precisions")
            # end
            # recalls::Dict{String, REAL} = jldopen(output_path, "r") do file
            #     read(file, "recalls")
            # end
            # fvalues::Dict{String, REAL} = jldopen(output_path, "r") do file
            #     read(file, "fvalues")
            # end
            # precisions[tag] = prescision
            # recalls[tag] = recall
            # fvalues[tag] = fvalue
        else
            precisions = Dict{String,REAL}(tag=>prescision)
            recalls    = Dict{String,REAL}(tag=>recall)
            fvalues    = Dict{String,REAL}(tag=>fvalue)
        end
        println(precisions)
        println(recalls)
        println(fvalues)
        jldopen(output_path, "w") do file
            file["precisions"] = precisions
            file["recalls"]    = recalls
            file["fvalues"]    = fvalues
        end
    end


    function _sort_biclustermatrix(matrix, usage_r, usage_c, B)
        cKeys = []
        for c in keys(usage_c)
            num = 0
            for r in keys(usage_r)
                num += (B[r,c] == 1)
                num -= (B[r,c] == 4)
            end
            append!(cKeys, tuple([c, num]))
        end
        sort!(cKeys, by = x -> x[2], rev = true)
        sort_r = []
        sort_c = []
        r_breaks = []
        c_breaks = []
        for x in cKeys
            append!(sort_c, usage_c[x[1]])
            push!(c_breaks, length(usage_c[x[1]]))
        end
        for (key, array) in usage_r
            append!(sort_r, array)
            push!(r_breaks, length(array))
        end
        for i in 1:(length(r_breaks)-1)
            r_breaks[i+1] += r_breaks[i]
        end
        for i in 1:(length(c_breaks)-1)
            c_breaks[i+1] += c_breaks[i]
        end
        return (matrix[sort_r, sort_c], r_breaks, sort_r, c_breaks, sort_c)
    end


    function _sort_biclustermatrix(matrix, usage_r, usage_c)
        sort_r = []
        sort_c = []
        r_breaks = []
        c_breaks = []
        for (key, array) in usage_r
            append!(sort_r, array)
            push!(r_breaks, length(array))
        end
        for (key, array) in usage_c
            append!(sort_c, array)
            push!(c_breaks, length(array))
        end
        for i in 1:(length(r_breaks)-1)
            r_breaks[i+1] += r_breaks[i]
        end
        for i in 1:(length(c_breaks)-1)
            c_breaks[i+1] += c_breaks[i]
        end
        return (matrix[sort_r, sort_c], r_breaks, sort_r, c_breaks, sort_c)
    end

    function _view_matrix_heatmap(matrix, r_breaks, sort_r, c_breaks, sort_c, title;
                                  file_path = "", clim = nothing)
        pyplot()
        (clim == nothing) && heatmap(yflip=true, xticks = (1:size(matrix,2),sort_c),
                             yticks = (1:size(matrix,1),sort_r), matrix, title = title, size=(2000,2000), c=:Blues)
        (clim != nothing) && heatmap(yflip=true, xticks = (1:size(matrix,2),sort_c),
                             yticks = (1:size(matrix,1),sort_r), matrix, title = title, clim = clim, size=(2000,2000), c=:Blues)
        hline!(r_breaks .+ 0.5)
        vline!(c_breaks .+ 0.5)
        if file_path != ""
            savefig(file_path)
        end
    end

    function _view_f_heatmap(f, usage_r, r_breaks, sort_r, usage_c, c_breaks, sort_c, title;
                             file_path = "")
        pyplot()
        matrix = zeros(Float64, length(sort_r), length(sort_c))
        X, Y = size(matrix)
        R = length(sort_r)
        C = length(sort_c)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = f[y]
            end
        end
        matrix = matrix[sort_r, sort_c]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sort_c),
                yticks = (1:size(matrix,1),sort_r), clim = (0.7, 1.0), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(r_breaks .+ 0.5)
        vline!(c_breaks .+ 0.5)
        if file_path != ""
            savefig(file_path)
        end
    end

    function _view_u_heatmap(u, usage_r, r_breaks, sort_r, usage_c, c_breaks, sort_c, title;
                             file_path = "")
        pyplot()
        matrix = zeros(Float64, length(sort_r), length(sort_c))
        X, Y = size(matrix)
        R = length(sort_r)
        C = length(sort_c)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = (u[y] == x)
            end
        end
        matrix = matrix[sort_r, sort_c]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sort_c),
                yticks = (1:size(matrix,1),sort_r), clim = (0,1), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(r_breaks .+ 0.5)
        vline!(c_breaks .+ 0.5)
        if file_path != ""
            savefig(file_path)
        end
    end

    function _view_a_heatmap(a, usage_r, r_breaks, sort_r, usage_c, c_breaks, sort_c, title;
                             file_path = "")
        pyplot()
        matrix = zeros(Float64, length(sort_r), length(sort_c))
        X, Y = size(matrix)
        R = length(sort_r)
        C = length(sort_c)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = a[x]
            end
        end
        matrix = matrix[sort_r, sort_c]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sort_c),
                yticks = (1:size(matrix,1),sort_r), clim = (1,2), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(r_breaks .+ 0.5)
        vline!(c_breaks .+ 0.5)
        if file_path != ""
            savefig(file_path)
        end
    end

    function _view_blocktype_heatmap(B, usage_r, r_breaks, sort_r, usage_c, c_breaks, sort_c, title;
                                     file_path = "")
        pyplot()
        matrix = zeros(Float64, length(sort_r), length(sort_c))
        R = length(sort_r)
        C = length(sort_c)
        for (r,c) in Iterators.product(keys(usage_r),keys(usage_c))
            for (i,j) in Iterators.product(usage_r[r],usage_c[c])
                matrix[i,j] = B[r,c]
            end
        end
        matrix = matrix[sort_r, sort_c]
        print(matrix)
        heatmap(yflip=true, xticks = (1:size(matrix,2),sort_c),
                yticks = (1:size(matrix,1),sort_r), clim = (1,4), matrix, title = title, size=(2000,2000), c=:pu_or)
        hline!(r_breaks .+ 0.5)
        vline!(c_breaks .+ 0.5)
        if file_path != ""
            savefig(file_path)
        end
    end

    function _get_mapstate(sampled)
        maxLik     = -Inf
        maxCountAt = 0
        for (s, count) in sampled
            if s.lnProb > maxLik
                maxCountAt = count
                maxLik = s.lnProb
            end
        end
        smax = nothing
        for (s, count) in sampled;
            if count == maxCountAt; smax = s; end
        end
        return smax
    end

    function _acf(s::Array{R, 1})::Array{R, 1} where { R <: Real }
        S = length(s)
        ans::Array{R, 1} = zeros(R, S)
        s_avg::R = sum(s) / S
        s_se::R  = sum( (s .- s_avg) .* (s .- s_avg) ) / (S - 1.0)

        for t in 1:(S-1)
            corr::R = 0.0
            for i in 1:(S-t)
                corr += (s[i] - s_avg) * (s[i+t] - s_avg)
            end
            corr = corr / (S-t)
            ans[t] = corr / s_se
        end
        return ans
    end

    function _evaluatefvalue(summary_path::String, score::Array{REAL, 2}, thres::REAL = (REAL)(0.0),
                             INT::Type{<:Integer} = Int32)::Tuple{REAL,REAL,REAL} where {REAL <: Real}
        summaryTemp::Array{REAL, 2} = InputParser.parse_input_summary(summary_path)
        summary::Array{REAL, 2}     = zeros(size(score))
        for (x,y) in Iterators.product(1:size(summaryTemp)[1],1:size(summaryTemp)[2])
            summary[x,y] = summaryTemp[x,y]
        end
        TP::INT = sum( (score .>  thres) .* (summary .>  0.0) )
        TN::INT = sum( (score .<= thres) .* (summary .<= 0.0) )
        FP::INT = sum( (score .>  thres) .* (summary .<= 0.0) )
        FN::INT = sum( (score .<= thres) .* (summary .>  0.0) )
        precision::REAL = TP / (TP + FP)
        recall::REAL = TP / (TP + FN)
        print("Precision: "); println(precision);
        print("Recall: ");    println(recall);
        fvalue::REAL = 2.0 * precision * recall / (precision + recall)
        print("Fvalue: ");    println(fvalue);
        return (precision, recall, fvalue)
    end

    function _evaluatefvalue(summary_path::String, sampled::Sampler{INT, REAL},
                             thres::REAL = (REAL)(0.0),
                             rm_error::Bool = true)::Tuple{REAL,REAL,REAL} where {INT <: Integer, REAL <: Real}
        score::Array{REAL, 2}       = convert.(REAL, sampled.Z) .- 1.0
        S::INT, M::INT = size(score)
        if rm_error
            for (i,j) in Iterators.product(1:S, 1:M)
                (sampled.er[j] == 2) && (score[i,j] = 0.0)
            end
        end
        return _evaluatefvalue(summary_path, score, thres, INT)
    end

end
