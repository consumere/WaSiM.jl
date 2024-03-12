# 
function fdi(;cwd::AbstractString=pwd(),xm::Regex=r"")
        dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
        for dir in dirs 
            if isdir(dir) & occursin(xm,dir)
                @printf("%-8s\t|","$dir")
            end
        end
    end
    

