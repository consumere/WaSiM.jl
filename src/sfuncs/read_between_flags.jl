# 
function read_between_flags(file::IOStream, flag1::String, flag2::String)
    #     line1 = readuntil(file, flag1)
    #     line2 = readuntil(file, flag2)
    #     return line1[1:end - length(flag1)], line2[1:end - length(flag2)]
    # end

    """
    skips first line after [soil_table] i.e. no of soil types
    now returns a DataFrame
    """
    