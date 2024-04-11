# 
function kge_df3()
        """
        should be non-recursive
        """
        x1=r"qoutjl"
        files = filter(file -> occursin(x1,file),
            readdir()[
                broadcast(x->!endswith(x,
                r"nc|pl|txt|svg|png|jpg|grd|ftz|ftz_0|list|xml|sh|yrly"), 
                readdir())])
        v = []
        for file_path in files
            if isfile(file_path)
                #df = waread(file_path);
                #df = fread(file_path)
                df = waread2(file_path)
                if !(startswith(names(df)[1],"C"))
                    nm=names(df)[1]
                    @warn "sim is $nm .. is it correct?"
                end
                dropmissing!(df)
                simulated = df[:,1]
                observed  = df[:,2]
                kge_value = kge1(simulated,observed)
                nse_value = nse(simulated,observed)
                nm = basename(file_path)
                println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
                printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
                push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
                v = DataFrame(v)
            end
        end
        return(v)
    end

    