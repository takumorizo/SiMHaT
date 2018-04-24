using JLD
const REAL = Float64
const INT  = Int64
using Plots
pyplot()


function viewInHeatMap2D(scores::Dict{Tuple{A, A} , REAL},
                         filePath::String; title::String = "",
                         xlabel::String = "", ylabel::String = "", clim = (0.50, 1.0), sorting::Bool = true) where{ A <: Any }
    aSet = Set{A}([])
    bSet = Set{A}([])
    abSet= Set{Tuple{A,A}}([])
    for (a,b) in keys(scores)
        @assert (a,b) ∉ abSet
        push!(aSet, a)
        push!(bSet, b)
        push!(abSet, (a,b))
    end
    as = collect(aSet); (sorting) && (sort!(as))
    bs = collect(bSet); (sorting) && (sort!(bs))
    toAsIndex = Dict{Any,Integer}()
    toBsIndex = Dict{Any,Integer}()
    for i in 1:length(as); toAsIndex[as[i]] = i; end;
    for i in 1:length(bs); toBsIndex[bs[i]] = i; end;

    ansScore::Array{REAL, 2} = zeros(REAL, length(as), length(bs))
    for (a,b) in keys(scores)
        ansScore[ toAsIndex[a], toBsIndex[b] ] = scores[a,b]
    end

    pyplot()
    println(as)
    println(bs)
    heatmap(as, bs, ansScore, xlabel = xlabel, ylabel = ylabel, title = title, c = :lightrainbow, clim = clim )
    savefig(filePath * ".png")
end

function toTupleDic(scores::Dict{String, REAL}, divStr::String = "")::Dict{Tuple{REAL, REAL}, REAL}
    ans::Dict{Tuple{REAL, REAL}, REAL} = Dict{Tuple{REAL, REAL}, REAL}()
    for tag in keys(scores)
        tags::Array{REAL, 1} = parse.(split(tag, divStr))
        ans[ tuple(tags...) ] = scores[tag]
    end
    return ans
end

function selectScores(scores::Dict{String, REAL};
                      prefix::String = "", suffix::String = "",
                      rmPrefix::Bool = true, rmSuffix::Bool = true)
    ansScore::Dict{String, REAL} = Dict{String, REAL}()
    for tag in keys(scores)
        if (startswith(tag, prefix) && endswith(tag, suffix))
            tagStart::INT = 1 + (INT)(rmPrefix) * (length(prefix))
            tagEnd::INT   = length(tag) - (INT)(rmSuffix) * (length(suffix))
            ansScore[tag[tagStart:tagEnd]] = scores[tag]
        end
    end
    return ansScore
end


function main()
    inputfile  = ARGS[1]
    prefix     = ARGS[2]
    suffix     = ARGS[3]
    outputFile = ARGS[4]

    precisions::Dict{String, REAL} = jldopen(inputfile, "r") do file
        read(file, "precisions")
    end
    recalls::Dict{String, REAL} = jldopen(inputfile, "r") do file
        read(file, "recalls")
    end
    fvalues::Dict{String, REAL} = jldopen(inputfile, "r") do file
        read(file, "fvalues")
    end

    println(selectScores(precisions, prefix = prefix, suffix = suffix))
    println(toTupleDic(selectScores(precisions, prefix = prefix, suffix = suffix), "_"))
    viewInHeatMap2D(toTupleDic(selectScores(precisions, prefix = prefix, suffix = suffix), "_"),
                    outputFile * ".precisions", title = "precisions", sorting = true, xlabel = "falseVAF", ylabel="trueVAF")

    println(selectScores(recalls, prefix = prefix, suffix = suffix))
    println(toTupleDic(selectScores(recalls, prefix = prefix, suffix = suffix), "_"))
    viewInHeatMap2D(toTupleDic(selectScores(recalls, prefix = prefix, suffix = suffix), "_"),
                    outputFile * ".recalls", title = "recalls", sorting = true, xlabel = "falseVAF", ylabel="trueVAF")

    println(selectScores(fvalues, prefix = prefix, suffix = suffix))
    println(toTupleDic(selectScores(fvalues, prefix = prefix, suffix = suffix), "_"))
    viewInHeatMap2D(toTupleDic(selectScores(fvalues, prefix = prefix, suffix = suffix), "_"),
                    outputFile * ".fvalues", title = "fvalues", sorting = true, xlabel = "falseVAF", ylabel="trueVAF")


    # viewInHeatMap2D(toTupleDic(selectScores(precisions, prefix, suffix, rmPrefix, rmSuffix)),
    #                 outputFile * ".precisions", title = "precisions", sorting = true)
    #
    # viewInHeatMap2D(toTupleDic(selectScores(recalls, prefix, suffix, rmPrefix, rmSuffix)),
    #                 outputFile * ".recalls", title = "recalls", sorting = true)
    #
    # viewInHeatMap2D(toTupleDic(selectScores(fvalues, prefix, suffix, rmPrefix, rmSuffix)),
    #                 outputFile * ".fvalues", title = "fvalues", sorting = true)
end

main()
