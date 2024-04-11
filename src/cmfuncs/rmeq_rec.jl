# 
function rmeq_rec(; rootdir = ".")
        
        ext_regex = r".R|.py|.jl|.tex|.pl|.sh|.csv|.html|.xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt"i
        ms = ["-9999", "lin", "log", "--"]
        
        files::Vector{String} = []       #String[]
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if !occursin(ext_regex, filename)
                    push!(files, joinpath(looproot, filename))
                end
            end
        end
        
        for inF in files
            if isfile(inF)
                println("reading ", inF, "...")
                df = CSV.File(inF; delim = "\t", header = 1,
                              silencewarnings = true, 
                              normalizenames = false, 
                              missingstring = ms, 
                              types = Float64) |> DataFrame
                if (isempty(df) || nrow(dropmissing(df,ncol(df)))==0 || ncol(df)==0) 
                    rm(inF)
                    println(basename(inF), " removed!")
                end
            end
        end
    end

    