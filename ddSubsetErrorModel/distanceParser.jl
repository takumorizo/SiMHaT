Include("config.jl")
Include("inputFileParser.jl")

module distanceParser
    using ..config
    using ..inputParser
    using ConfParser
    function parseBFHammingDistance(errScorePath::String,
                                    patScorePath::String,
                                    matScorePath::String,
                                    paramFilePath::String,
                                    alphaName::String)::Array{REAL, 2}
        errScore::Array{REAL, 2} = inputParser.parseInputSummary(errScorePath)
        matScore::Array{REAL, 2} = inputParser.parseInputSummary(matScorePath)
        patScore::Array{REAL, 2} = inputParser.parseInputSummary(patScorePath)

        conf = ConfParse(paramFilePath)
        parse_conf!(conf)
        BFThreshold::REAL = parse(REAL, String(retrieve(conf, "distance", "BFThreshold")))
        alpha::REAL       = parse(REAL, String(retrieve(conf, "distance", alphaName)) )

        S::INT, M::INT = size(errScore)
        BFPat::Array{REAL, 2} = patScore .- errScore
        BFMat::Array{REAL, 2} = matScore .- errScore
        BF::Array{REAL, 2}    = max.(BFPat, BFMat)

        mutMat::Array{INT, 2} = map(BF) do x
           if x > BFThreshold
               return (INT)(1)
           else
               return (INT)(0)
           end
       end
       ans::Array{REAL, 2} = zeros(REAL, S, S)
       for (from, to) in Iterators.product((INT)(1):S, (INT)(1):S)
           dist::REAL = (REAL)(0.0)
           for m in (INT)(1):M
               dist += (REAL)( (mutMat[from, m] != mutMat[to, m]) )
           end
           ans[to, from] = dist / M
       end
       for i in (INT)(1):S
           ans[i, i] = alpha
       end
       return ans
    end

    function parsePhysicalDistance(distanceFilePath::String)::Array{REAL, 2}
        distanceMatrix::Array{REAL, 2} = inputParser.parseInputSummary(distanceFilePath)
        return distanceMatrix
    end

    function parseDecayFunction(paramFilePath::String)
        conf = ConfParse(paramFilePath)
        parse_conf!(conf)
        functionType::String = String(retrieve(conf, "distance", "decayFunction"))
        a::REAL      = parse(REAL, String(retrieve(conf, "distance", "decayRate")))
        if      functionType == "Exponential"
            return x->(REAL)(exp(-1.0 * x / a))
        elseif  functionType == "Sigmoid"
            return x->(REAL)(( exp(-d + a) / (1 + exp(-d+a)) ))
        elseif  functionType == "Window"
            return x->( (REAL)( x < a ) )
        elseif  functionType == "Uniform"
            return x->(REAL)(1)
        else
            error("unexpected type of decay functions!")
        end
    end
end
