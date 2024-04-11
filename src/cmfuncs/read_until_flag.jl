# 
function read_until_flag(file::IOStream, flag::String)
        line = readuntil(file, flag)
        return line[1:end-length(flag)]
    end

    