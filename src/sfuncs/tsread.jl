# 
function tsread(x::Union{String,Regex};kw...)
    if x isa String
        printstyled("reading $x\n",color=:light_red)
    else x isa Regex
        x = first(dfonly(x))
        printstyled("reading $x\n",color=:light_red)
    end
    #ms = ["-9999","-9999.0","lin", "log", "--"]
    #missingstring=ms,
    #df = CSV.read(x,DataFrame;kw...)
    df = CSV.File(x;kw...)|>DataFrame|>z->dropmissing(z,1)
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
Fastest Reader. is also dfr.
Read the text file, preserve line 1 as header column
"""
