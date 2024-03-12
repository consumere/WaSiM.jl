# 
function findindf(df::DataFrame, x::Regex)
        """
        like Grep.grep("x",df)
        """
        filter(row -> any(occursin(x, 
            string(value)) for value in row), 
                eachrow(df))
    end

    