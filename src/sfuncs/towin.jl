# 
function towin(file_path::Union{String, Nothing})
    if isnothing(file_path)
        lw = split(pwd(), '/')[3]
        rst = join(split(pwd(), '/')[4:end], '/')
        win_path = string(uppercase(lw), ":/", rst)
    else
        lw = split(file_path, '/')[3]
        rst = join(split(file_path, '/')[4:end], '/')
        win_path = string(uppercase(lw), ":/", rst)
    end
    return win_path
end

"""
windows path to wsl
"""
