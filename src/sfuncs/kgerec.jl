# 
function kgerec(ext::String)
        """
        should be recursive
        see also kge_rec
        """
        path::String = pwd()
        results = []
        v = []
        for (looproot, dirs, filenames) in walkdir(path)
            for filename in filenames
                #if (endswith(filename, ext)) 
                if (occursin(Regex(ext,"i"),filename)) && 
                    (!occursin(r"output|yrly|nc|png|svg|jpg|xml|ctl",filename))
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
        sz=length(results)
        println("found $sz files...\n$results")
        for file in results
                    dd = CSV.read(file,DataFrame,missingstring="-9999",delim="\t")
                    observed  = dd[:,5]
                    simulated = dd[:,6]
                    kge_value = kge2(observed, simulated)
                    nse_value = nse(observed, simulated)
                    nm = basename(file)
                    printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
                    println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
                    push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm,:path=>file))
                    v = DataFrame(v)
        end
        return(v)
    end

    