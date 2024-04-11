# 
function dfilter(df::DataFrame, col, val)
        DataFrames.subset(df,Symbol(col)=> ByRow(==(val)))
    end

    """
    qplot from regex
    """
    