# 
function findindf(df::DataFrame, x::Regex)
        """
        like Grep.grep("x",df)
        """
        filter(row -> any(occursin(x, 
            string(value)) for value in row), 
                eachrow(df))
    end

    # findindf(df,"Et")
    # findindf(df,"15")
    # findindf(df,"-")
    # Grep.grep(r"Et",df.name)
    # Grep.grep(r"Et",df)
    # Grep.grep("Et",df)


    