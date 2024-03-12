# 
function latx()
    directory = pwd()
    files = readdir(directory)  
    file_times = Dict{String, Dates.DateTime}()
    for file in files
        file_path = joinpath(directory, file)
        stat_info = stat(file_path)
        file_times[file] = Dates.unix2datetime(stat_info.mtime)
    end
    #xf = sort(file_times, by = x -> x[2], rev=true)
    
    file_times = collect(file_times)
    
    xf = sort(file_times, by = x -> x[2], rev=true)
    xf = first(xf,11)
    sz = []                 #::Vector{Float64}
    for i in 1:length(xf)
       push!(sz,round(stat(xf[i][1]).size/1024^2,digits=4))
    end
    xf2 = merge_vectors(sz,xf)
    for (size, file, dt) in xf2
        datetime = Dates.format(dt, "yyyy-mm-dd HH:MM")
        printstyled(rpad(file,35),color=:yellow),
        printstyled("$datetime\t" ,color=:green),
        printstyled(rpad(size,7)," MB\n",color=:yellow)    
    end
end

ll=latx

"""
grabs methods
asin|>getm  
?asin
@code_llvm readf|>getm|>first 
"""
