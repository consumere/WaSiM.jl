# 
function wcl(x::AbstractString)
        files = filter(file -> occursin(Regex(x,"i"),file), readdir())
        for file in files
            open(file) do f
                ct=(count(_ -> true, eachline(f)))
                #println(file,ct)
                println("$file:\t $ct")
            end
        end
    end

    