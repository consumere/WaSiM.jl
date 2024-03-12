# 
function fz()
        """
        gets sorted DF by size recursivley
        """
        cwd = pwd() 
        osize = 0
        m = []
        for (root, dirs, files) in walkdir(cwd)
        for file in files
            if isfile(file)
                nm=joinpath(root, file)
                osize = stat(nm).size/1024^2
                push!(m,Dict(:name=>file,
                :size=>osize,
                :fullnaname=>nm))
            end
        end 
        end
        df = DataFrame(m)     
        sort!(df, [order(:size,rev = true), order(:name)])
        return(df)
    end 

    