# 
function kge_df(ext::Union{String,Regex};path=pwd())
        
        if ext isa Regex
            files = filter(file -> occursin(ext,file), readdir())
        else
            files = filter(file -> endswith(file, ext), readdir())
        end
        v = []
        for file in files
            file_path = joinpath(path, file)
            if isfile(file_path) && endswith(file, ext)
                dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
                simulated = dd[:,5]
                observed  = dd[:,6]
                kge_value = kge2(simulated,observed)
                nse_value = nse(simulated, observed)
                ve_value = vef(simulated, observed)
                nm = basename(file_path)
                println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
                println(replace("VE value is $ve_value on $nm", "\\"  => "/"))
                printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
                push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:VE=>ve_value,:name=>nm))
                v = DataFrame(v)
            end
        end
        return(v)
    end

    