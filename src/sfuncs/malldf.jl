# 
function malldf(files::Vector{DataFrame},on::Symbol)
        "reduces + merges by date"
        df = reduce((left, right) -> 
        innerjoin(left, right, on = on,makeunique=true), 
        files)
        df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        return(df)
    end

    