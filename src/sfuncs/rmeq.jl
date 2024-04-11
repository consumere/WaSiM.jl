# 
function rmeq()
    """
    removes empty TS; 
    use with caution!
    """
    #x = pwd()

    # files = filter(file -> (occursin(Regex(x, "i"), file) & 
    # (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file))
    # ), readdir())
    
    files = filter(file -> 
    !occursin(r"tex|pl|sh|csv|html|xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file)
    , readdir())
    
    ms = ["-9999", "lin", "log", "--"]
    for inF in files
        if isfile(inF)
            df = CSV.File(inF; delim="\t", header=1,
            silencewarnings=true, 
                normalizenames=false, 
                missingstring=ms, 
                types=Float64) |> DataFrame
            dropmissing!(df,ncol(df))
            if nrow(df)==0
                println(basename(inF)," removed!")
                rm(inF)
            end
        end
    end
end

