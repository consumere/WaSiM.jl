# 
function glob(x::Regex)
    filter(file -> occursin(x,file), readdir())
end

"""
fdi(;xm::Regex=r"*")
lists dirs if isdir(dir) & occursin(xm,dir)
"""
