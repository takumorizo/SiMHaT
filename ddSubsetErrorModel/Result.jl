Include("config.jl")
Include("BuffPhyloMatrixType.jl")
Include("SamplerType.jl")
Include("InputParser.jl")

module Result
    using Plots
    using ..SamplerType
    using JLD2, FileIO
    using ..config
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
    function viewmap(smax, lnprobs, outputDir::String)
        pyplot()
        iters = collect(1:length(lnprobs))
        plot(iters, lnprobs, size=(2000,2000))
        savefig(outputDir * "/lnprobs.png")

        acf = _acf(lnprobs)
        bar(acf)
        savefig(outputDir * "/acf.png")

        println("isValid state")
        println(SamplerType.isvalid(smax, debug = true))
        mat, rbreak, sortR, cbreak, sortC = _sort_biclustermatrix(smax.Z, smax.usage_s, smax.usage_v, smax.B)

        # println("MAP, B")
        # println(smax.B)
        # println("isError valid matrix")
        # println(BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(buffPhyloMatrix.BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix{Int64})
        # isTree(smax.treeCache)

        _view_matrix_heatmap(mat, rbreak, sortR, cbreak, sortC, "MAP.Z", filePath = outputDir * "/MAP_Z.png")
        _view_matrix_heatmap(smax.ln_p_data[sortR,sortC,1], rbreak, sortR, cbreak, sortC, "err", filePath = outputDir * "/err.png")
        _view_matrix_heatmap(smax.ln_p_data[sortR,sortC,2], rbreak, sortR, cbreak, sortC, "mat", filePath = outputDir * "/mat.png")
        _view_matrix_heatmap(smax.ln_p_data[sortR,sortC,3], rbreak, sortR, cbreak, sortC, "pat", filePath = outputDir * "/pat.png")

        _view_matrix_heatmap(smax.ln_p_data[:,:,1], rbreak, sortR, cbreak, sortC, "err", filePath = outputDir * "/err.nonSort.png")
        _view_matrix_heatmap(smax.ln_p_data[:,:,2], rbreak, sortR, cbreak, sortC, "mat", filePath = outputDir * "/mat.nonSort.png")
        _view_matrix_heatmap(smax.ln_p_data[:,:,3], rbreak, sortR, cbreak, sortC, "pat", filePath = outputDir * "/pat.nonSort.png")
        _view_matrix_heatmap(smax.H[sortR,sortC], rbreak, sortR, cbreak, sortC, "MAP.H", filePath = outputDir * "/MAP_H.png", clim = [1,2])
        _view_f_heatmap(smax.f, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.f", filePath = outputDir * "/MAP_f.png")
        _view_u_heatmap(smax.u, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.u", filePath = outputDir * "/MAP_u.png")
        _view_a_heatmap(smax.a, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.a", filePath = outputDir * "/MAP_a.png")
        _view_blocktype_heatmap(smax.B, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.B", filePath = outputDir * "/MAP_B.png")
    end

    function viewdata(sampled, outputDir::String)
        pyplot()

        cutoff = 200
        lnprobs = []; iters   = []; for (s,c) in sampled; push!(lnprobs, s.lnProb); push!(iters, c); end; plot(iters[cutoff:length(iters)], lnprobs[cutoff:length(iters)], size=(2000,2000))
        savefig(outputDir * "/lnprobs.png")

        smax = _get_mapstate(sampled)
        println("isValid state")
        println(SamplerType.isvalid(smax, debug = true))
        mat, rbreak, sortR, cbreak, sortC = _sort_biclustermatrix(smax.Z, smax.usage_s, smax.usage_v)

        # println("MAP, B")
        # println(smax.B)
        # println("isError valid matrix")
        # println(BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(buffPhyloMatrix.BuffPhyloMatrix{Int64} <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix)
        # println(typeof(smax.treeCache) <: BuffPhyloMatrix{Int64})
        # isTree(smax.treeCache)

        _view_matrix_heatmap(mat, rbreak, sortR, cbreak, sortC, "MAP.Z", filePath = outputDir * "/MAP_Z.png")
        _view_matrix_heatmap(smax.ln_p_data[sortR,sortC,1], rbreak, sortR, cbreak, sortC, "err", filePath = outputDir * "/err.png")
        _view_matrix_heatmap(smax.ln_p_data[sortR,sortC,2], rbreak, sortR, cbreak, sortC, "mat", filePath = outputDir * "/mat.png")
        _view_matrix_heatmap(smax.ln_p_data[sortR,sortC,3], rbreak, sortR, cbreak, sortC, "pat", filePath = outputDir * "/pat.png")

        _view_matrix_heatmap(smax.H, rbreak, sortR, cbreak, sortC, "MAP.H", filePath = outputDir * "/MAP_H.png", clim = [1,2])
        _view_f_heatmap(smax.f, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.f", filePath = outputDir * "/MAP_f.png")
        _view_u_heatmap(smax.u, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.u", filePath = outputDir * "/MAP_u.png")
        _view_a_heatmap(smax.a, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.a", filePath = outputDir * "/MAP_a.png")
        _view_blocktype_heatmap(smax.B, smax.usage_s, rbreak, sortR, smax.usage_v, cbreak, sortC, "MAP.B", filePath = outputDir * "/MAP_B.png")
    end

    function evaluatefvalue(summaryPath::String, scorePath::String, thres::REAL = (REAL)(0.0), rmError::Bool = true)::Tuple{REAL,REAL,REAL}
        if endswith(scorePath, "jld2")
            sampled::Sampler{INT, REAL}, lnprobs::Array{REAL, 1} = readjld(scorePath)
            return _evaluatefvalue(summaryPath, sampled, thres, rmError)
        else
            score::Array{REAL, 2}       = InputParser.parseInputSummary(scorePath)
            return _evaluatefvalue(summaryPath, score, thres)
        end
    end

    function savefvalue(outputPath::String, prescision::REAL, recall::REAL, fvalue::REAL, tag::String)
        @assert endswith(outputPath, "jld2")
        if isfile(outputPath)
            println(load(outputPath))
            precisions::Dict{String, REAL} = load(outputPath, "precisions")
            recalls::Dict{String, REAL}    = load(outputPath, "recalls")
            fvalues::Dict{String, REAL}    = load(outputPath, "fvalues")
            precisions[tag] = prescision
            recalls[tag]    = recall
            fvalues[tag]    = fvalue
            # precisions::Dict{String, REAL} = jldopen(outputPath, "r") do file
            #     read(file, "precisions")
            # end
            # recalls::Dict{String, REAL} = jldopen(outputPath, "r") do file
            #     read(file, "recalls")
            # end
            # fvalues::Dict{String, REAL} = jldopen(outputPath, "r") do file
            #     read(file, "fvalues")
            # end
            # precisions[tag] = prescision
            # recalls[tag] = recall
            # fvalues[tag] = fvalue
        else
            precisions = Dict{String,REAL}(tag=>prescision)
            recalls   = Dict{String,REAL}(tag=>recall)
            fvalues   = Dict{String,REAL}(tag=>fvalue)
        end
        println(precisions)
        println(recalls)
        println(fvalues)
        jldopen(outputPath, "w") do file
            file["precisions"] = precisions
            file["recalls"]    = recalls
            file["fvalues"]    = fvalues
        end
    end


    function _sort_biclustermatrix(matrix, usageR, usageC, B)
        cKeys = []
        for c in keys(usageC)
            num = 0
            for r in keys(usageR)
                num += (B[r,c] == 1)
                num -= (B[r,c] == 4)
            end
            append!(cKeys, tuple([c, num]))
        end
        sort!(cKeys, by = x -> x[2], rev = true)
        sortR = []
        sortC = []
        rBreaks = []
        cBreaks = []
        for x in cKeys
            append!(sortC, usageC[x[1]])
            push!(cBreaks, length(usageC[x[1]]))
        end
        for (key, array) in usageR
            append!(sortR, array)
            push!(rBreaks, length(array))
        end
        for i in 1:(length(rBreaks)-1)
            rBreaks[i+1] += rBreaks[i]
        end
        for i in 1:(length(cBreaks)-1)
            cBreaks[i+1] += cBreaks[i]
        end
        return (matrix[sortR, sortC], rBreaks, sortR, cBreaks, sortC)
    end


    function _sort_biclustermatrix(matrix, usageR, usageC)
        sortR = []
        sortC = []
        rBreaks = []
        cBreaks = []
        for (key, array) in usageR
            append!(sortR, array)
            push!(rBreaks, length(array))
        end
        for (key, array) in usageC
            append!(sortC, array)
            push!(cBreaks, length(array))
        end
        for i in 1:(length(rBreaks)-1)
            rBreaks[i+1] += rBreaks[i]
        end
        for i in 1:(length(cBreaks)-1)
            cBreaks[i+1] += cBreaks[i]
        end
        return (matrix[sortR, sortC], rBreaks, sortR, cBreaks, sortC)
    end

    function _view_matrix_heatmap(matrix, rBreaks, sortR, cBreaks, sortC, title;
                              filePath = "", clim = nothing)
        pyplot()
        (clim == nothing) && heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                             yticks = (1:size(matrix,1),sortR), matrix, title = title, size=(2000,2000), c=:Blues)
        (clim != nothing) && heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                             yticks = (1:size(matrix,1),sortR), matrix, title = title, clim = clim, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function _view_f_heatmap(f, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                            filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        X, Y = size(matrix)
        R = length(sortR)
        C = length(sortC)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = f[y]
            end
        end
        matrix = matrix[sortR, sortC]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (0.7, 1.0), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function _view_u_heatmap(u, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                            filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        X, Y = size(matrix)
        R = length(sortR)
        C = length(sortC)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = (u[y] == x)
            end
        end
        matrix = matrix[sortR, sortC]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (0,1), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function _view_a_heatmap(a, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                            filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        X, Y = size(matrix)
        R = length(sortR)
        C = length(sortC)
        for y in 1:Y
            for x in 1:X
                matrix[x,y] = a[x]
            end
        end
        matrix = matrix[sortR, sortC]
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (1,2), matrix, title = title, size=(2000,2000), c=:Blues)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
        end
    end

    function _view_blocktype_heatmap(B, usageR, rBreaks, sortR, usageC, cBreaks, sortC, title;
                                     filePath = "")
        pyplot()
        matrix = zeros(Float64, length(sortR), length(sortC))
        R = length(sortR)
        C = length(sortC)
        for (r,c) in Iterators.product(keys(usageR),keys(usageC))
            for (i,j) in Iterators.product(usageR[r],usageC[c])
                matrix[i,j] = B[r,c]
            end
        end
        matrix = matrix[sortR, sortC]
        print(matrix)
        heatmap(yflip=true, xticks = (1:size(matrix,2),sortC),
                yticks = (1:size(matrix,1),sortR), clim = (1,4), matrix, title = title, size=(2000,2000), c=:pu_or)
        hline!(rBreaks .+ 0.5)
        vline!(cBreaks .+ 0.5)
        if filePath != ""
            savefig(filePath)
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

    function _evaluatefvalue(summaryPath::String, score::Array{REAL, 2}, thres::REAL = (REAL)(0.0))::Tuple{REAL,REAL,REAL}
        summaryTemp::Array{REAL, 2} = InputParser.parseInputSummary(summaryPath)
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

    function _evaluatefvalue(summaryPath::String, sampled::Sampler{INT, REAL},
                             thres::REAL = (REAL)(0.0),
                             rmError::Bool = true)::Tuple{REAL,REAL,REAL}
        score::Array{REAL, 2}       = convert.(REAL, sampled.Z) .- 1.0
        S::INT, M::INT = size(score)
        if rmError
            for (i,j) in Iterators.product(1:S, 1:M)
                (sampled.er[j] == 2) && (score[i,j] = 0.0)
            end
        end
        return _evaluatefvalue(summaryPath, score, thres)
    end

end
