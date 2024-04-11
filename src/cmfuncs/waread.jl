# 
function waread(x::Regex)
        x = dfonly(x)|>first
        ms = ["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame; delim="\t", header=1, missingstring=ms, normalizenames=true, types=Float64)
        df = dropmissing(df, 1)
        dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
        df.date = dt2
        df = select(df, Not(1:4))
        metadata!(df, "filename", x, style=:note)
        #renamer
        for x in names(df)
            if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    