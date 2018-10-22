Include("InputParser.jl")

module DistanceParser
    using ..InputParser
    using ConfParser
    function parse_bf_hamming_distance(err_score_path::String,
                                       pat_score_path::String,
                                       mat_score_path::String,
                                       param_file_path::String,
                                       alpha_name::String,
                                       INT::Type{<:Integer} = Int32,
                                       REAL::Type{<:Real}   = Float32)::Array{REAL, 2}
        errScore::Array{REAL, 2} = InputParser.parse_input_summary(err_score_path)
        matScore::Array{REAL, 2} = InputParser.parse_input_summary(mat_score_path)
        patScore::Array{REAL, 2} = InputParser.parse_input_summary(pat_score_path)

        conf = ConfParse(param_file_path)
        parse_conf!(conf)
        BFThreshold::REAL = parse(REAL, String(retrieve(conf, "distance", "BFThreshold")))
        alpha::REAL       = parse(REAL, String(retrieve(conf, "distance", alpha_name)) )

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

    function parse_physical_distance(distance_file_path::String, INT::Type{<:Integer} = Int32, REAL::Type{<:Real} = Float32)::Array{REAL, 2}
        distanceMatrix::Array{REAL, 2} = InputParser.parse_input_summary(distance_file_path, INT, REAL)
        return distanceMatrix
    end

    function parse_decayfunction(param_file_path::String;
                                 distance_tag::String       = "distance",
                                 decay_function_tag::String = "decayFunction",
                                 decay_rate_tag::String     = "decayRate",
                                 INT::Type{<:Integer} = Int32,
                                 REAL::Type{<:Real}   = Float32)
        conf = ConfParse(param_file_path)
        parse_conf!(conf)
        functionType::String = String(retrieve(conf, distance_tag, decay_function_tag))
        a::REAL      = parse(REAL, String(retrieve(conf, distance_tag, decay_rate_tag)))
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
