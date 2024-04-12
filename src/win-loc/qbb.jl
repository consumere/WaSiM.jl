@time using DataFrames: DataFrame, Cols, Not, permutedims, dropmissing, rename!, nrow, hcat
@time using CSV: File
#, metadata!
#using DataFrames

function rglob(prefix::String)
    rootdir="."
    results = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (occursin(Regex(prefix,"i"),filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    return results
end

"""
rglob qgko
return a vector of DFs
"""
function qbb()
    files = rglob("qgko")
    dfs = []
    for file in files
        x = file
        try
            df = DataFrame(File(x, #CSV.File(x,
                            header=false, 
                            delim="\t",
                            #ntasks = 1,
                            limit=5000,
                            #skipto=2,
                            silencewarnings=true,
                            ignorerepeated=true,
                            #debug=true,
                            types = String))

            println("check $x ...")
            pattern = r"^[LIN. R]|^[LOG. R]|^CO"
            mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
            dx = df[mask, :]
            dx = permutedims(dx) |>dropmissing
            basins = []
            for i in copy(df[1,5:end])
                push!(basins,string.("B_"*i))
            end
            insert!(basins, 1, "score")
            insert!(basins, 2, "timestep")
            dx[!, "basin"] = basins
            cn = (dx[1,:])
            rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
            dout = dx[3:end,:]
            dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
            #metadata!(dout, "filename", file, style=:note);
            dout.nm .= file
            push!(dfs,dout)
            
        catch e
            @error "$e"
            @warn("error! file $x can not be loaded as a DataFrame! ")
            # Skip files that can't be loaded as a DataFrame
            continue
        end
    end
    return(dfs)
end

@time dfs = qbb(ARGS...)
println(dfs)
println("done!")
