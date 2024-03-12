# 
function freaddf(ext::AbstractString)
        cwd = pwd() 
        m = DataFrame[]
        for (root, dirs, files) in walkdir(cwd)
        for file in files
        if isfile(file) && occursin(Regex(ext),file)&&
            (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",file))
        nm=joinpath(root, file)
        push!(m,readdf(nm))
        end
        end 
        end 
        return(m)
    end 

    """
    reads, reduces + merges by date and plots
    pall(glob(r"qges|qbas|qd"))
    """
    