# 
function readall(path::Union{Vector{Any},Vector{String}})
        files = dfonly(path)
        dfs::Vector{DataFrame} = []
        for file in files
            if isfile(file) && (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
                file_path = file
                println("reading ",file_path,"...")
                #p1 = waread(file_path)
                ##modified waread
                ms = ["-9999","lin","log","--"]
                df = CSV.read(file_path, DataFrame; 
                    delim="\t", header=1, missingstring=ms, 
                    #maxwarnings = 1, 
                    silencewarnings = true,
                    rows_to_check = 100, ignorerepeated = true,
                    normalizenames=true, types=Float64)
                df = dropmissing(df, 1)
                dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
                df.date = dt2
                df = select(df, Not(1:4))
                DataFrames.metadata!(df, "filename", file_path, style=:note)
                for x in names(df)
                    if startswith(x,"_")
                        newname=replace(x,"_"=>"C", count=1)
                        rename!(df,Dict(x=>newname))
                    end
                end
                
                push!(dfs, df)
            end
        end
        return(dfs)
    end

    