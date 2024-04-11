# 
function tovec(x::DataFrame, col::Any)
        df = select(x,col)
        println(names(df))
        return vec(Matrix(df))
    end

    # """
    # prints out the 