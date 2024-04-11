# 
function colsums(df::DataFrame)
        colsum = sum.(eachcol(df[!,Not(Cols(r"date|month|year"))]))
        return colsum
    end

    """
    Subset a DataFrame to exclude only columns with zero sums.
    """
    