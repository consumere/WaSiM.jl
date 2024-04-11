# 
function wqplot(file_path::AbstractString)
        data = read_wq(file_path)
        p = @df data plot(data[!,1],cols(propertynames(data)[2:end]))
        println(describe(data))
        return p
    end

    """
    reads recursively
    """
    