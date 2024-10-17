# 
function dfrdt(x::Union{Regex,String})
    if x isa(Regex)
        try 
            x = first(dfonly(x))
        catch
            @error "no match for $x ! "
        end
    end
    ms = ["-9999","lin","log","--"]
    df = CSV.read(x, DataFrame; 
        delim="\t", header=1, missingstring=ms, 
        maxwarnings = 1, #silencewarnings = true,
        normalizenames=true, types=Float64)
    df = dropmissing(df, 1)
    dt2 = map(row -> Dates.DateTime(Int(row[1]), 
        Int(row[2]), Int(row[3]), 
        Int(row[4]),0,0,0),  #note that: Dates.DateTime,dateformat"yyyy mm dd HH MM SS s"
        eachrow(df))
    df.date = dt2
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
