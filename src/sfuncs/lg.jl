# 
function lg(path::AbstractString, ext::AbstractString)
        files = readdir(path)
        v=[]
        for file in files
            file_path = joinpath(path, file)
            if isfile(file_path) && endswith(file, ext)
            println(file_path)
        push!(v,file_path)
            end
        end
        return(v)
    end

    