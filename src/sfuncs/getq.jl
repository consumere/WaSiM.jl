# 
function getq(;prefix::AbstractString="qgko")
    rootdir = pwd()
    results = []
    if any(x -> isdir(x), readdir())
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if occursin(Regex(prefix, "i"), filename) && !occursin(r"txt|yrly|nc|png|svg", filename)
                    push!(results, joinpath(looproot, filename))
                end
            end
        end
    else
        printstyled("eval on: $rootdir !\n", color=:light_red)
        for filename in filter(x -> isfile(x), readdir(; join=false))
            if occursin(Regex(prefix, "i"), filename) && !occursin(r"txt|yrly|nc|png|svg", filename)
                printstyled("collecting $filename...\n", color=:light_yellow)
                push!(results, filename)
            end
        end
    end

    if isempty(results)
        printstyled("no qgk files found!\n", color=:light_red)
        return 
    end

    for file in results
        x = file
        try
            df = CSV.read(x, DataFrame, delim="\t"; header=true, types=String, silencewarnings=true, skipto=364)
            pattern = r"^[LIN. R]|^[LOG. R]|^CO"
            mask = occursin.(pattern, df[!, 1])
            ddd = df[mask, :]
            new = names(ddd)[5:end]
            insert!(new, 1, "basin")
            insert!(new, 2, "timestep")
            ddd = permutedims(ddd)
            dropmissing!(ddd)
            ddd.basin = new
            select!(ddd, :basin, :)
            df = permutedims(ddd)
            columns = uppercasefirst.(getproperty(df, propertynames(df)[1]))
            rename!(ddd, Symbol.(columns))
            kd = ddd[3:end, :]
            return kd
        catch e
            println(e)
            @warn "skipping $x"
        end
    end

    return 
end

"""

tries to parse all columns to float, except of date.
"""
