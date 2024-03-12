# 
function subsum(df::DataFrame)
        column_sums = sum.(eachcol(df[!, Not(:date)]))
        # Find the column indices with sums equal to 0.0
        zero_sum_columns = findall(==(0.0), column_sums)
        # Subset the DataFrame to include only columns with zero sums
        dout = hcat(df.date, select(df, Not(zero_sum_columns)))
        rename!(dout, 1=>"date")
        return dout
    end

    """
    Returns the total size of all files in the current directory non recursively.
    """
    