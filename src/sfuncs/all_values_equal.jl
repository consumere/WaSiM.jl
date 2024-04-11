# 
function all_values_equal(df::DataFrame)
        for col in eachcol(df)
            if any(col .!= col[1])
                return false
            end
        end
        return true
    end

    """
    reduces + merges by date + plots all
    pall(files::Vector{DataFrame};toyr=true,leg=:outertopright)
    """
    