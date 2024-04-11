# 
function parse_date(row)
        try
            year = parse(Int, row[1])
            month = parse(Int, row[2])
            day = parse(Int, row[3])
            return Date(year, month, day)
        catch
            return missing  # Return missing if parsing fails
        end
    end
    
    df.date = map(parse_date, eachrow(df))
    
    df = select(df, Not(1:4))
    
    for x in names(df)
        if looks_like_number(x)
            newname=replace(x,"$x"=>"C$x", count=1)
            rename!(df,Dict(x=>newname))
        end
    end

    #map(y->typeof(y),eachcol(df))

    # Iterate over column names
    for colname in names(df)
        # Check if the column type is not Date
        if eltype(df[!, colname]) != Date
            df[!, colname] .= tryparse.(Float64, df[!, colname])
        end
    end

    DataFrames.metadata!(df, "filename", x, style=:note)
    return df 
end

