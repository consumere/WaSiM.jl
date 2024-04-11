# 
function fz()
    cwd = pwd() 
    osize = 0
#    fn = 0
    m = []
    for (root, dirs, files) in walkdir(cwd)
     for file in files
        if isfile(file)
           nm=joinpath(root, file)
           osize = stat(nm).size/1024^2
#	       fn += stat(nm).size
	       #push!(m,(nm))
           push!(m,Dict(:name=>file,
           :size=>osize,
#           :total=>(fn/1024^2),
           :fullnaname=>nm))
        end
    end 
end
    df = DataFrame(m)     
    sort!(df, [order(:size,rev = true), order(:name)])
    return(df)
end 

