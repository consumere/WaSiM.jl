# 
function pew()
    try
        in = clipboard()
        wpath = replace(in, "\\" => "/")
        println("pt=$wpath")
        clipboard("pt=$wpath")
        return wpath
    catch e
        @error "smth errored $e"
    end
end

