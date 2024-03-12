# 
function rec()
        #fls = rglob("qout")
        needle = r"qout"
        rootdir = pwd()
        results = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if (occursin(needle,filename)
                    && 
                    !occursin(r"yr|mon|grid|scn"i,filename) && 
                    !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly|eps)$", filename)
                    )
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
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
            push!(out,Dict("name"=>x,"NSE"=>nse2(dd),"KGE"=>kge2(dd),"VE"=>vef2(dd)))
            catch
                @warn("$x can't be loaded as a DataFrame, skipping ...")
                continue
            end
            
        end
        df = DataFrame(out)
        df.NSE .= replace(df.NSE, NaN => missing)
        df.KGE .= replace(df.KGE, NaN => missing)
        df.VE .= replace(df.VE, NaN => missing)
        dropmissing!(df)
        df.bn .= map(x->splitdir(x)|>last,df.name)
        df.pt .= map(x->splitdir(dirname(x))|>last,df.name)
        return sort(df,:NSE;rev = true)
    end

    """
    tries to read all files in a directory as df
    """
    