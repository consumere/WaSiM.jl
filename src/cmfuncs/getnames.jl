# 
function getnames(dfs::Vector{Any})
        try 
            nms=[]
            for df in dfs
                x=collect(DataFrames.metadata(df))[1][2]|>basename
                println(x)
                push!(nms,x)
            end
            return nms
        catch
            @error "no metadata in $df !"
            return
        end
    end

    