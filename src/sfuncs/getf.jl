# 
function getf(ext::AbstractString)
        cwd = pwd() 
        m = []
        for (root, dirs, files) in walkdir(cwd)
        for file in files
        if isfile(file) && occursin(Regex(ext),file)
        nm=joinpath(root, file)
        push!(m,(nm))
        end
        end 
        end 
        return(m)
    end 

    # getf(".*(^th)+.*(nc)+.*")  
    # #SAME
    # getf("^th+.*nc")
    # ###lookbehind	
    # #getf("stack?+.*nc") 
    # #getf("!stack?+.*nc") 

    