# 
function old_waread2(x::String)
        """
        Read the text file, preserve line 1 as header column
        """
        ms=["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame, 
            delim="\t",
            header=1,
            missingstring=ms,
            normalizenames=true,
            types=Float64)
        dropmissing!(df,1)
        dt2::Vector{Date} = []
        for i in eachrow(df)
            z=i[1:3]|>collect|>transpose
            push!(dt2,Date(z[1],z[2],z[3]))
        end
        df.date = dt2
        df=df[:,Not(1:4)]
        metadata!(df, "filename", x, style=:note);
    end

    