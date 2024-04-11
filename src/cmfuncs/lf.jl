# 
function lf()
        """
        list_files_sorted_by_last_change
        formerly lat()
        """
        directory = pwd()
        files = readdir(directory)  
        file_times = Dict{String, Dates.DateTime}()
        for file in files
            file_path = joinpath(directory, file)
            stat_info = stat(file_path)
            file_times[file] = Dates.unix2datetime(stat_info.mtime)
        end   
        file_times = collect(file_times)
        xf = sort(file_times, by = x -> x[2], rev=true)
        #xf = first(xf,11)
        sz = []                 #::Vector{Float64}
        for i in 1:length(xf)
           push!(sz,round(stat(xf[i][1]).size/1024^2,digits=4))
        end
        #dt = Dates.unix2datetime(stat(first(xf[1])).mtime)
        #Dates.format(dt, "yyyy-mm-dd HH:M")
        merged_vector = []
        for (item1, item2) in zip(xf,sz)
            merged_item = (item1.first, item1.second,item2)
            push!(merged_vector, merged_item)
        end
        
        # for (file,datetime,size) in merged_vector
        #     printstyled(rpad(file,35),color=:yellow),
        #     printstyled("$datetime\t" ,color=:green),
        #     printstyled(rpad(size,7)," MB\n",color=:yellow)    
        # end
        df = DataFrame(merged_vector);
        rename!(df,["latest_file","modified","size_in_mb"]);
        df.modified=map(x -> Dates.format(x, "yyyy-mm-dd HH:MM"),df.modified)
        return(df)
    end
    
    