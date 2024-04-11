# 
function kge_read(path::AbstractString, ext::AbstractString)
        files = readdir(path)
        for file in files
            file_path = joinpath(path, file)
            if isfile(file_path) && endswith(file, ext) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|tif|jpg|png|svg",file))
                dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
                # observed  = dd[:,5]
                # simulated = dd[:,6]
                # kge_value = kge2(observed, simulated)
                kge_value = kge2(dd)
                println(replace("KGE value is $kge_value on $file_path", "\\"  => "/"))
            elseif isdir(file_path)
                dfs_in_subdir = kge_read(file_path, ext)
            end
        end
    end

    """
    sim, obs
    hydroeval.objective_functions.nse(simulations, evaluation)
    """
    