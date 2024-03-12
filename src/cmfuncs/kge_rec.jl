# 
function kge_rec()
        """
        reads recursively
        """
        #fls = rglob("qout")
        needle = r"qout"
        rootdir = pwd()
        results = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if (occursin(needle,filename)
                    && 
                    !occursin(r"yr|mon|grid|scn"i,filename) && 
                    !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|xml|sh|grd|yrly|eps)$", filename)
                    )
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
        # results = filter(x -> isfile(x) && 
        # !occursin(r"yr|mon|grid|scn"i,x) && 
        # !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|xml|sh|grd|yrly|eps)$", x), results)
        out = []
        for x in results
            println("reading $x ...")
            try
                dd = CSV.read(x,DataFrame,
                missingstring="-9999",
                #maxwarnings=1,
                silencewarnings=true,
                ignorerepeated=true,
                delim="\t")
            push!(out,Dict("name"=>x,"KGE"=>kge2(dd)))
            catch
                @warn("$x can't be loaded as a DataFrame, skipping ...")
                continue
            end
            
            
        end
        df = DataFrame(out)
        df.KGE .= replace(df.KGE, NaN => missing)
        dropmissing!(df)
        return sort(df,:KGE;rev = true)
    end

    