# 
function fdi(cwd::AbstractString)   
        dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
        for dir in dirs
            if isdir(dir)
                @printf("%-8s\t|","$dir");
            end
        end
    end

    