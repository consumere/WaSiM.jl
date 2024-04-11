# 
function mall(left::DataFrame, right::DataFrame)
        "reduces + merges by date"
        df = innerjoin(left, right, on = :date,makeunique=true)
        return(df)
    end

    