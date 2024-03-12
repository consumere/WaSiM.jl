# 
function tocb(s::Union{String,Regex})
        if isa(s,Regex)
            s = Regex(s,"i")
        end
        y = first(filter(file -> occursin(s,file), readdir()))
        println("abspath of $y in clipboard!")
        abspath(y)|>cb
    end

    # 