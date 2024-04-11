# 
function malldf(files::Vector{DataFrame},on::Symbol)
        "reduces + merges by date"
        df = reduce((left, right) -> 
        innerjoin(left, right, on = on,makeunique=true), 
        files)
        return(df)
    end

    