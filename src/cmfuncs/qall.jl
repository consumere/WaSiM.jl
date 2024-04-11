# 
function qall(;recursive=false)
        if recursive
            files = rglob("qgko")
        else
            files = glob("qgko")
        end
        outdf::Vector{DataFrame} = []
        for file in files
            # Load the file into a DataFrame
            #x = "qgkofab.m6.2010"
            x = file
            try
                df = DataFrame(CSV.File(x, header=1, 
                                    delim="\t",
                                    skipto=366,
                                    ntasks=1,
                                    types = String,
                                    silencewarnings=true,
                                    ignorerepeated=true,
                                    ignoreemptyrows=true,
                                    stripwhitespace=true))
                                    
                #df = CSV.read(x,DataFrame;ntasks=1)
                println(x)
                pattern = r"^[LIN. R]|^[LOG. R]|^CO"
                mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
                dx = df[mask, :]
                dx = permutedims(dx) |>dropmissing
                
                #mapcols!(x -> parse(Float64, x), dx)
                #mapcols(x -> x.^2, dx)

                #basins = copy(df[1,5:end])
                #AsTable(basins)
                basins = []
                #for i in copy(df[1,5:end])
                for i in names(df)[5:end]
                    push!(basins,i)
                end
                #size(basins)
                insert!(basins, 1, "score")
                insert!(basins, 2, "timestep")
                dx[!, "basin"] = basins

                cn = (dx[1,:])
                rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
                dout = dx[3:end,:]
                for col in names(dout)[1:end-1]
                    dout[!, col] = parse.(Float64, replace.(dout[!, col], "," => ""))
                end          
                #mapcols!(x -> parse(Float64, x), dout)
                #dout = hcat(dx[!,Cols(r"bas")],dx[:,Not(Cols(r"bas"))])
                dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
                dout.basin = parse.(Int64,dout[!, :basin])
                dout.nm .= file
                push!(outdf,dout)
            catch
                @warn("error! files that can't be loaded as a DataFrame")
                # Skip files that can't be loaded as a DataFrame
                continue
            end
        end
        return(outdf)
    end

    