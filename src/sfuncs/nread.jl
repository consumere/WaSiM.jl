# 
function nread(x::Union{String,Regex};kw...)
        if x isa String
            printstyled("reading $x\n",color=:light_red)
        else x isa Regex
            x = first(dfonly(x))
            printstyled("reading $x\n",color=:light_red)
        end
        ms = ["-9999","-9999.0","lin", "log", "--"]
        df = CSV.read(x,DataFrame;missingstring=ms,kw...)
        if "YY" ∉ names(df)
            println("Column 'YY' not found in the CSV file.")
            @show first(df,5)
            return nothing
        end
        if !all(i -> i isa Int64, df.YY)
            filtermask = broadcast(x->looks_like_number(x),df[!,1])
            df = df[filtermask, :]        
        end
        date_strings = string.(df.YY, "-", df.MM, "-", df.DD, "-", df.HH)
        df.date = DateTime.(date_strings, "yyyy-mm-dd-HH")
        df = select(df, Not(1:4))
        DataFrames.metadata!(df, "filename", x, style=:note)
        for x in names(df)
            if startswith(x,"_")
                newname=replace(x,"_"=>"C", count=1)
                rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    """
    sim obs plot rglob regex
    x::Union{Regex,DataFrame}; 
        yscale::Symbol = :log,
        fun::Function = sum, 
        freq::String="monthly",simcol=1,obscol=2
        freq can also shortened, like D,M,Q,Y
    """
    