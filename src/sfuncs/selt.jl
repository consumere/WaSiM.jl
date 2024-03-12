# 
function selt(x::DataFrame, col::Any; dtcol=:date)
        df = select(x,Cols(col,dtcol))
        println(names(df))
        return df   #vec(Matrix(df))
    end

    """
    hydrgraph plot, df meta approach
    """
    