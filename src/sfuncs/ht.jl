# 
function ht(df::DataFrame)
        mapcols(x -> x[Not(2:end-1)], df)
    end

    """
    Plots a histogram of the values in the DataFrame
    """
    