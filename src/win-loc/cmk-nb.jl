### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ cdb4bd40-4199-11ee-1eb9-19c27eed6305
begin
		using CairoMakie, Makie    
		using DataFrames, CSV, Statistics, Dates, Distributions
		using DelimitedFiles, Grep , Printf
	end

# ╔═╡ 478b18ee-c9b5-46f0-abd5-fa22a8731bbf
module cmk
    using CairoMakie, Makie    
    using DataFrames, CSV, Statistics, Dates, Distributions
    using DelimitedFiles, Grep , Printf
    # using PrettyTables
    # using SHA
    # using PyCall

    function qgk(;rootdir=".", prefix="qgk")
        """
        filters internal WaSiM stats of routed discharge files
        works recursively
        """
        files = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                    push!(files, joinpath(looproot, filename)) 
                end
            end
        end
        
        for z in files
            println(raw"file:	",basename(z),"...")
            m = filter(line -> occursin(r"^[LIN. R]|^[LOG. R]|^CO",line), readlines(open(z)))
            for l in m
                x = replace(l, r"\s+" => "\t")
                x = replace(x, ".\t" => " ")
                println(x)
            end
        end
        return nothing
    end


    function jdd(;return_string=true)
        """
        """
        cwd = pwd()
        dirs = readdir(".")
        vst = []
        for dir in dirs
            if isdir(dir)
                size = 0
                for (root, dirs, files) in walkdir(dir)
                    for file in files
                        size += stat(joinpath(root, file)).size
                    end
                end
                @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2)
                if return_string
                    #v = string.(cwd,"\\",dir,": ",size/1024^2," MB\n")
                    v = hcat(string.(cwd,"\\",dir,),size/1024^2)
                    push!(vst,v)
                end   
            end
        end
      return(vst)
    end

    function dd()
        cwd = pwd()
        osize = 0
        for (root, dirs, files) in walkdir(cwd)
        for file in files
            osize += stat(joinpath(root, file)).size
        end
        end 
        @printf("%-40s %15.3f GB\n","$(cwd):",osize/1024^3);
    end 

    function ct(ext::AbstractString)
        cwd = pwd() 
        osize = 0
        fz = 0
        m = []
        for (root, dirs, files) in walkdir(cwd)
        for file in files
        if isfile(file) && occursin(Regex(ext),file)
        nm=joinpath(root, file)
        osize = stat(nm).size
        @printf("%-40s %15.2f MB\n","$(nm):",osize/1024^2);
        fz += stat(nm).size
        push!(m,(nm))
        end
        end 
        end 
        n=repeat(" - -",10)
        println(n*" sum of ",ext*n)
        @printf("%-40s %15.2f MB\n","$(cwd):",fz/1024^2);
        println(n,length(m)," matches "*n,"\n")
        return(m)
    end 

    function ct()
        cwd = pwd() 
        osize = 0
        fz = 0
        m = []

        for (root, dirs, files) in walkdir(cwd)
        for file in files
        if isfile(file)
            nm=joinpath(root, file)
            osize = stat(nm).size
            #sizes[file] = stat(nm).size/1024^2
            @printf("%-40s %15.2f MB\n","$(nm):",osize/1024^2);
            fz += stat(nm).size
            push!(m,(nm))
        end     
        end 
        end 
        # sort(df, [order(:a), order(:b, rev = true)]) 
        # sorted_files = sort(collect(keys(sizes)), by=x->sizes[x], rev=false)
        # println(sorted_files)
        n=repeat(" - -",10)
        println(n*" total sum ")
        @printf("%-40s %15.2f MB\n","$(cwd):",fz/1024^2);
        println(n,length(m)," files present "*n,"\n")
        return(m)
    end

    function print_sorted_sizes(dir)
        """
        print_sorted_sizes("../")
        """
        folders = [joinpath(dir, f) for f in readdir(dir)]
        sizes = Dict()
        for f in folders
            if isdir(f)
                sizes[f] = get_folder_size(f)
            end
        end
        sorted_folders = sort(collect(keys(sizes)), by=x->sizes[x], rev=false)
        for f in sorted_folders
            if sizes[f] >= 1000000
            printstyled(rpad(f,60, ' '), rpad(sizes[f] ÷ 10^6, 6, ' '), "MB\n",color=:green)
            end
        end
    end

    function get_folder_size(folder)
        files = readdir(folder)
        size = 0
        for file in files
            path = joinpath(folder, file)
                if isfile(path)
                    size += stat(path).size
                elseif isdir(path)
                    size += get_folder_size(path)
            end
        end
        return size
    end

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

    function readdf(x::Regex;recurse=false)
        """
        readdf(x::Regex)
        reads first match of regex wasim timeseries
        """
        if recurse
            results = []
            rootdir="."
            for (looproot, dirs, filenames) in walkdir(rootdir)
                for filename in filenames
                    if (occursin(x,filename)) && (!occursin(r"txt|yrly|nc|png|svg|grd",filename))
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
            if length(results)==0
                @error "no match for $x !"
                return
            end
            fn=first(results)
        else
            fn=try 
                first(filter(file -> (occursin(x,file) & 
                (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",
                file))), readdir()))
                #That’s a speedup of about 2x, but x is a string
                # regex = Regex(x * "(?!\\.(xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg))","i")
                # files = [file for file in readdir() if occursin(regex, file)]
                # first(files) #1st match
            catch
                    @error "no match for $x !"
                    return
            end
        end
        x = fn
        ms=["-9999","lin","log","--"]
        df = CSV.read(x,DataFrame,
        missingstring = ms, 
        # ntasks=4,
        limit = typemax(Int),
        types = Float64,
        ignorerepeated = true,
        delim="\t",
        silencewarnings=true,
        stripwhitespace=true,
        normalizenames=true,
        drop=(i, nm) -> i == 4)
        dropmissing!(df,1)
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,Not(1:3)]
        metadata!(df, "filename", x, style=:note);
        
        #s = (filter(x->!occursin(r"year|date",x),names(df)))
        #renamer - remove char _   
        for x in names(df)
            if startswith(x,"_")
            #newname=replace(x,"_"=>"C")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        #names(df)
        #s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
        return df 
    end

    function readdf(x::AbstractString)
        """
        --- main reader ---
        delim: if no argument is provided, 
        parsing will try to detect the most consistent delimiter on the 
            first 10 rows of the file
        """
        ms=["-9999","lin","log"]
        df::DataFrame = CSV.read(x,DataFrame,
        missingstring=ms,
        #ignorerepeated = true,
        #delim="\t",
        #skipto=4,
        types = Float64,
        silencewarnings=false,
        normalizenames=true,
        drop=(i, nm) -> i == 4) #|> dropmissing
        dropmissing!(df,1)
        #df.DD  = map(x -> begin val = tryparse(Int, x); 
        #ifelse(typeof(val) == Nothing, missing, val) end, df.DD )
        #map(x ->Int(x),df[!,1:3])
        #map(x ->round(x;digits=0),df.YY)
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        #df.YY=map(x ->Date(x,"yyyy"),df.YY);
        #dropmissing!(df)
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,Not(1:3)]
        DataFrames.metadata!(df, "filename", x, style=:note);
            #s = (filter(x->!occursin(r"year|date",x),names(df)))
        #renamer - remove char _   
        for x in names(df)
            if startswith(x,"_")
            #newname=replace(x,"_"=>"C")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end


    function old_waread2(x::String)
        """
        Read the text file, preserve line 1 as header column
        """
        ms=["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame, 
            delim="\t",
            header=1,
            missingstring=ms,
            normalizenames=true,
            types=Float64)
        dropmissing!(df,1)
        dt2::Vector{Date} = []
        for i in eachrow(df)
            z=i[1:3]|>collect|>transpose
            push!(dt2,Date(z[1],z[2],z[3]))
        end
        df.date = dt2
        df=df[:,Not(1:4)]
        metadata!(df, "filename", x, style=:note);
    end
    function waread2(x::String)
        """
        Read the text file, preserve line 1 as header column
        Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
        This avoids reading the entire file into memory at once, 
        which can be more memory-efficient for large datasets.
        """
        ms = ["-9999", "lin", "log", "--"]
        df = CSV.File(x; delim="\t", header=1, normalizenames=true, missingstring=ms, types=Float64) |> DataFrame
        dropmissing!(df,1)
        dt2 = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
        select!(df, Not(1:4))
        df.date = dt2
        metadata!(df, "filename", x, style=:note)
        return df
    end

    function waread(x::String)
        """
        Fastest Reader. is also dfr.
        Read the text file, preserve line 1 as header column
        """
        ms = ["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame; delim="\t", header=1, missingstring=ms, normalizenames=true, types=Float64)
        df = dropmissing(df, 1)
        dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
        df.date = dt2
        df = select(df, Not(1:4))
        DataFrames.metadata!(df, "filename", x, style=:note)
        for x in names(df)
            if startswith(x,"_")
                newname=replace(x,"_"=>"C", count=1)
                rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    function waread(x::Regex)
        """
        Read the text file, preserve line 1 as header column
        """
        x = dfonly(x)|>first
        ms = ["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame; delim="\t", header=1, missingstring=ms, normalizenames=true, types=Float64)
        df = dropmissing(df, 1)
        dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
        df.date = dt2
        df = select(df, Not(1:4))
        metadata!(df, "filename", x, style=:note)
        #renamer
        for x in names(df)
            if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    function loadalldfs(path::AbstractString)
        files = readdir(path)
        dfs = DataFrame[]
        for file in files
            if isfile(file) && (!occursin(
                r"^wq|ftz_0|tex|year|xml|qgk|fzt|ftz|log|ini|txt|yrly|nc|png|svg",file))
                file_path = joinpath(path, file)
                println("reading ",file_path,"...")
                p1 = waread(file_path)
                push!(dfs, p1)
            end
        end
        return(dfs)
    end

    function loadalldfs(path::Vector{Any})
        files = path
        dfs = DataFrame[]
        for file in files
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = file
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        push!(dfs, p1)
            end
        end
        return(dfs)
    end

    function loadalldfs(files::Vector{String})
        dfs = DataFrame[]
        for file in files
            if isfile(file) && (!occursin(r"tex|xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = file
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        push!(dfs, p1)
            end
        end
        return(dfs)
    end

    function loadalldfs(path::Regex)
        v::Vector{String} = readdir();
        v = v[broadcast(x->!endswith(x,"nc"),v)];
        files = v[(broadcast(x->occursin(path,x),v))];
        dfs::Vector{DataFrame} = []
        for file in files
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = file
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        push!(dfs, p1)
            end
        end
        return(dfs)
    end

    function loadso(path::AbstractString, prefix::AbstractString)
        files = readdir(path)
        dfs = DataFrame[]
        for file in files
            if isfile(file) && occursin(Regex(prefix),file)&& (!occursin(r"fzt|fzs|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = joinpath(path, file)
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        push!(dfs, p1)
            end
        end
        return(dfs)
    end

    function loadalldfs(path::AbstractString)
        files = readdir(path)
        dfs = DataFrame[]
        #nms = []
        for file in files #&& occursin(Regex(prefix),file)
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = joinpath(path, file)
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        push!(dfs, p1)
        #push!(nms, file)
            end
        end
        return(dfs)
        #return(nms)
    end


    function listdfs(path::AbstractString)
        files = readdir(path)
        #dfs = DataFrame[]
        nms = []
        for file in files #&& occursin(Regex(prefix),file)
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = joinpath(path, file)
        println("reading ",file,"...")
        #p1 = readdf(file_path)
        #push!(dfs, p1)
        push!(nms, file)
            end
        end
        #return(dfs)
        return(nms)
    end

    function qpl(df::DataFrame)
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
        end
        s = names(df)[1:2]
        t2 = string.(ti,"\n",s[1],"|",s[2],ti)
        StatsPlots.plot( 
        qqplot(df[!,1],df[!,2], qqline = :fit), 
        qqplot(Cauchy,df[!,2]), 
        qqnorm(df[!,2], qqline = :R),
        title = t2)
    end

    function qpl(x::AbstractString)
        df = readdf(x)
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
        end
        s = names(df)[1:2]
        t2 = string.(ti,"\n",s[1],"|",s[2],ti)
        StatsPlots.plot( 
        qqplot(df[!,1],df[!,2], qqline = :fit), 
        qqplot(Cauchy,df[!,2]), 
        qqnorm(df[!,2], qqline = :R),
        title = t2)
    end

    qqp=qpl

    # function vibx(df::String)
    #     df = readdf(df)
    #     str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    #     ln = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
    #     @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
    #     @df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.5,marker=(:black,stroke(1)),legend=false)
    # end

    # function vibx(df::DataFrame)
    #     str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    #     ln = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    #     @df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
    #     @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
    #     @df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
    # end

    # function vio(mm::Regex)
    #     """
    #     vioplot wasim timeseries 
    #     """
    #     df = readdf(mm)
    #     str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    #     ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #     @warn "No basename in metadata!"
    #     ti = raw""
    #     end    
    #     if (any(x->occursin("year",x),names(df)))
    #         #s = Symbol.(filter(x->!occursin("year",x),names(df)))
    #         df = df[!,Not("year")]
    #         s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #         @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
    #         #@df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
    #         title!(ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
    #     #@df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
    #     title!(ti)
    #     #@df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
    #     end
    # end

    # function vio(df::String)
    #         df = readdf(df)
    #         str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    #         ti = try
    #             DataFrames.metadata(df)|>only|>last|>basename
    #         catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #         end    
    #         if (any(x->occursin("year",x),names(df)))
    #             #s = Symbol.(filter(x->!occursin("year",x),names(df)))
    #             df = df[!,Not("year")]
    #             s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #             @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
    #             #@df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
    #             title!(ti)
    #         else    
    #         s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #         @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
    #         #@df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
    #         title!(ti)
    #         #@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.15,legend=false)
    #         #marker=(:black,stroke(1)),legend=false)
    #     end
    # end

    function ldf(path::AbstractString, prefix::AbstractString)
        files = readdir(path)
        dfs = DataFrame[]
        for file in files
            if isfile(file) && occursin(Regex(prefix),file)&& (!occursin(r"txt|yrly|nc|png|svg",file))
            file_path = joinpath(path, file)
        println("reading",file_path)
        p1 = loaddf(file_path)
        push!(dfs, p1)
            end
        end
        return(dfs)
    end

    function kge2(simulated::Vector{Float64}, observed::Vector{Float64})
        r = cor(simulated, observed)
        α = std(simulated) / std(observed)
        β = mean(simulated) / mean(observed)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    function kge2(df::DataFrame)
        observed, simulated = df[:,6],df[:,5]
        r = cor(observed, simulated)
        α = std(simulated) / std(observed)
        β = mean(simulated) / mean(observed)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    function kge_read(path::AbstractString, ext::AbstractString)
        files = readdir(path)
        for file in files
            file_path = joinpath(path, file)
            if isfile(file_path) && endswith(file, ext) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|tif|jpg|png|svg",file))
                dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
                # observed  = dd[:,5]
                # simulated = dd[:,6]
                # kge_value = kge2(observed, simulated)
                kge_value = kge2(dd)
                println(replace("KGE value is $kge_value on $file_path", "\\"  => "/"))
            elseif isdir(file_path)
                dfs_in_subdir = kge_read(file_path, ext)
            end
        end
    end

    function nse(predictions::Vector{Float64}, targets::Vector{Float64})
        return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
    end

    function nse(df::DataFrame)
        simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
        return (1 - (sum((simulated .- observed).^2) / sum((observed .- mean(observed)).^2)))
    end

    function kge(df::DataFrame)
        simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
        r = cor(observed, simulated)
        α = std(simulated) / std(observed)
        β = mean(simulated) / mean(observed)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    # function lplot(x::Regex)
    #     df=readdf(x)
    #     ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #     ti = raw"" 
    #     end
    #     #o = collect(DataFrames.metadata(df))[1][2] |>basename
    #     ln = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df Plots.plot(:date,cols(ln),yaxis=:log,title=ti)
    # end

    # function dfl(regex::Regex)
    #     "selects first match and plots in log y-axis..."
    #     df=readf(regex)
    #     ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #     ti = raw"" 
    #     end
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
    #         @df df Plots.plot(:year,cols(s),yaxis=:log,legend = :topright, title=ti)
    #     else   
    #         s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #         @df df Plots.plot(:date,cols(s),yaxis=:log, legend = :topright, title=ti)
    #         end
    # end

    # function dfl!(regex::Regex)
    #     "adds first match and plots in log y-axis..."
    #     df=readf(regex)
    #     ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #     ti = raw"" 
    #     end
    #     println("adding $ti")
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
    #         @df df Plots.plot!(:year,cols(s),yaxis=:log,legend = :topright)
    #         Plots.annotate!([(20,5,text(ti, 12, :left, :top, :green))])
    #     else   
    #         s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #         @df df Plots.plot!(:date,cols(s),yaxis=:log, legend = :topright)
    #         Plots.annotate!([(20,5,text(ti, 12, :left, :top, :green))])
    #         end
    # end
   
    function tree()
        dir = pwd()
        dirs = readdir(dir)
        routes::Vector{String} = []
        for directory in dirs
            if isdir("$dir/" * directory)
                push!(routes, "$dir/$directory")
            else
                if ~(directory in routes)
                    push!(routes, "$routes/$directory")
                    #newread = dir * "/$directory"
                    #push!(newread, "$dir/$directory")
                    #[push!(routes, r) for r in newrs]
                end
            end
        end
        routes
    end

    function fdi()
        cwd = pwd()
        dirs = readdir(cwd)
        for dir in dirs
            if isdir(dir)
                @printf("%-8s\t|","$dir");
            end
        end
    end

    function fdi(cwd::AbstractString)   
        dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
        for dir in dirs
            if isdir(dir)
                @printf("%-8s\t|","$dir");
            end
        end
    end

    function ll()
        readdir()
    end

    function ls()
        p=pwd()
        f=readdir()
        dirs=filter(x -> (isdir(x)),f)
        files=filter(x -> (isfile(x)),f)
        if length(files)>0 && length(dirs)>0
            nf=length(files)
            println("$p\ndirs: $dirs\n $nf files:\n$files")
        elseif length(dirs)>0
            println("$p\ndirs: $dirs\n")
        elseif (length(files)>0 && length(files)<=12)
            nf=length(files)
            println("$p\n$nf files:\n $files\n")
        elseif (length(files)>0 && length(files)>12)
            nf=length(files)
            println("$p\n$nf files:\n",first(f,6),"...",last(f,6))
        else
            println("$p\n")
        end
    end

    # function lplot(df::DataFrame)
    #     nm = propertynames(df)[1:end-1];
    #     o = DataFrames.metadata(df)|>collect
    #     ti = basename(o[1][2])
    #     @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=ti)     
    # end

    # function lplot(df::String)
    #     df=readdf(df)
    #     nm = propertynames(df)[1:end-1];
    #     o = collect(DataFrames.metadata(df))[1][2] |>basename
    #     @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=o)     
    # end

    # function lplotf(df::String)
    #     df=readdf(df)
    #     nm = propertynames(df)[1:end-1];
    #     o = collect(DataFrames.metadata(df))[1][2] |>basename
    #     @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=o)     
    # end

    function getnames(dfs::Vector)
        nms = [];
        for i in dfs;	
            x=collect(DataFrames.metadata(i))[1][2]|>basename
            push!(nms, x)
        end
        return(nms)
    end

    function getnames(rg::String,dfs::Vector)
        nms = [];
        for i in dfs;	
            x=collect(DataFrames.metadata(i))[1][2]|>basename
            push!(nms, x)
        end
        return(filter(x->occursin(Regex(rg,"i"),x),nms))
    end

    function getnames(dfs::DataFrame)
        x=collect(DataFrames.metadata(dfs))[1][2]|>basename
        return(x)
    end

    function loadalldfs(path::AbstractString)
        files = readdir(path)
        dfs = DataFrame[]
        #nms = []
        for file in files #&& occursin(Regex(prefix),file)
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = joinpath(path, file)
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        push!(dfs, p1)
        #push!(nms, file)
            end
        end
        return(dfs)
        #return(nms)
    end

    function listdfs(path::AbstractString)
        files = readdir(path)
        #dfs = DataFrame[]
        nms = []
        for file in files #&& occursin(Regex(prefix),file)
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = joinpath(path, file)
        println("reading ",file,"...")
        #p1 = readdf(file_path)
        #push!(dfs, p1)
        push!(nms, file)
            end
        end
        #return(dfs)
        return(nms)
    end

    function vars()
        varinfo()
    end

    function vars(pt::AbstractString)
        #varinfo(Core,r".*field.*")
        #varinfo(Main,r".*load*")
        varinfo(Main,Regex(".*pt*"))
    end

    #if (occursin(Regex(prefix,"i"),filename))

    function regand(v::Vector{String},x1::AbstractString,y1::AbstractString)
        needle=join([x1,y1],"+.*");
        z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))] 
    return(z)
    end

    
    function regand(v::Vector{String},xv::Tuple{String, String})
        needle=join([xv[1],xv[2]],"+.*");
        z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
    return(z)
    end

    #function regand(v::Vector{String},xv::Tuple{Symbol,Symbol})
    function regand(v::Vector{String},xv::Vector{Symbol})
        needle=join(xv,"+.*");
        z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
    return(z)
    end

    function regand(v::Vector{String},xv::Vector{String})
        needle=join(xv,"+.*");
        z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
    return(z)
    end

    function regand(a::String, b::String)
        """
        basic a + b
        """
        needle=join([a,b],"+.*");
        z=Regex(needle,"i")
        return(z)
    end

    function regand(v::Vector{Any},a::String, b::String)
        """
        here you can put any string to filter the Vector
        like regand(getnames(dfs),"scn","ssr")
        """
        needle=Regex(join([a,b],"+.*"),"i")
        z = v[(broadcast(x->occursin(needle,x),v))] 
        if length(z)==1
            return(only(z))
        else
            return(z)
        end
    end


    function regand(v::Vector{Any},xv::Regex)
        """
        here you can put any regex to filter the Vector
        like regand(getnames(dfs),r"tem")
        """
        z = v[(broadcast(x->occursin(xv,x),v))] 
    return(z)
    end


    function nconly(x1::AbstractString)
        v::Vector{String} = readdir();
        v = v[broadcast(x->endswith(x,"nc"),v)];
        z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
        return(z)
    end

    function dfonly(x1::AbstractString)
        v = filter(file -> occursin(Regex(x1,"i"),file), readdir());
        #z = v[broadcast(x->!endswith(x,r"xml|nc|png|svg|jpg"),v)];
        z = v[broadcast(x->!occursin(r"^wq|txt|xml|nc|png|svg|jpg",x),v)];
        return(z)
    end
    
    function dfonly(x1::Regex)
        z = filter(file -> occursin(x1,file), 
        readdir()[broadcast(x->!occursin(r"^wq|xml|nc|png|svg|jpg",x),readdir())]);       
        # readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);       
        return(z)
    end

    function dfon(x1::AbstractString)
        z = filter(file -> occursin(Regex(x1,"i"),file), 
        readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
        return(z)
    end

    function nconly(Any)
        v = readdir();
        z = v[broadcast(x->endswith(x,"nc"),v)];
        return(z)
    end


    # function dfp(df::String)
    #     df=readdf(df)
    #     o = DataFrames.metadata(df)|>collect
    #     ti = basename(o[1][2])
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin("year",x),names(df)))
    #         @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
    #     end
    # end

    # function dfp(regex::AbstractString,dfs::Vector{DataFrame})
    #     "selects first match and plots..."
    #     df = dfs[map(n->occursin(Regex(regex,"i"),n),
    #         map(x->basename(only(DataFrames.metadata(x))[2]),
    #         dfs))] |> first
    #     o = DataFrames.metadata(df)|>collect
    #     ti = basename(o[1][2])
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin("year",x),names(df)))
    #         @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
    #     end
    # end


    # function dfp!(df::DataFrame)
    #     ti = try
    #             DataFrames.metadata(df)|>only|>last|>basename
    #         catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
    #         @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    #     @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
    #     end
    # end

    # function dfp!(df::String)
    #     df=readdf(df)
    #     o = DataFrames.metadata(df)|>collect
    #     ti = basename(o[1][2])
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin("year",x),names(df)))
    #         @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
    #     end
    # end

    # function dfp!(mm::Regex)
    #     """
    #     plots wasim timeseries
    #     """
    #     df=readdf(mm)
    #     o = DataFrames.metadata(df)|>collect
    #     ti = basename(o[1][2])
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin("year",x),names(df)))
    #         @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
    #     end
    # end

    # function dfp!(regex::AbstractString,dfs::Vector{DataFrame})
    #     "selects first match and plots..."
    #     df = dfs[map(n->occursin(Regex(regex,"i"),n),
    #         map(x->basename(only(DataFrames.metadata(x))[2]),
    #         dfs))] |> first
    #     o = DataFrames.metadata(df)|>collect
    #     ti = basename(o[1][2])
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin("year",x),names(df)))
    #         @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #     @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
    #     end
    # end


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

    function getdf(ext::AbstractString)
        cwd = pwd() 
        m = []
        for (root, dirs, files) in walkdir(cwd)
        for file in files
        if isfile(file) && occursin(Regex(ext),file)&&(!occursin(r"txt|yrly|nc|png|svg",file))
        nm=joinpath(root, file)
        push!(m,(nm))
        end
        end 
        end 
        return(m)
    end 


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

    # function plotf(ext::AbstractString)
    #     cwd = pwd() 
    #     m = []
    #     for (root, dirs, files) in walkdir(cwd)
    #     for file in files
    #     if isfile(file) && occursin(Regex(ext),file)&&(!occursin(r"txt|yrly|nc|png|svg",file))
    #     nm=joinpath(root, file)
    #     push!(m,(nm))
    #     end
    #     end 
    #     end 
    #     return(
    #     dfp(readdf(m[1])))
    # end 

    # function plotf(ext::String)
    #     dfp(readdf(ext))
    #     plot!(title=basename(ext))
    # end 

    # function plotf(ext::DataFrame)
    # dfp(ext)
    # end 


    # function homg()
    #     pt="D:/Wasim/Goldbach/revision/"
    #     cd(pt)
    #     println("you are here: ",pwd())
    # end


    # function hombr()
    #     pt="D:/Wasim/Tanalys/DEM/brend_fab/out/m4/"
    #     cd(pt)
    #     println("you are here: ",pwd())
    # end

    # function hometeo()
    #     cd("D:/Wasim/Tanalys/DEM/Input_V2/meteo/")
    #     println("you are here: ",pwd())
    # end

    # function homreg()
    #     cd("D:/Wasim/regio/out/rc200/");
    #     println("you are here: ",pwd())
    #     fdi()
    # end

    function homg()
        cd("D:/Wasim/Goldbach/");
        println("you are here: ",pwd())
        fdi()
    end

    function writewa(file::AbstractString, df::DataFrame)
        """
        old, but working version
        newer see wawrite
        """
        dout = df
        dout.YY = map(x ->year(x),dout.date)
        dout.MM = map(x ->month(x),dout.date)
        dout.DD = map(x ->day(x),dout.date)
        dout[!, "HH"] .= 0
        #df = select!(df,Symbol.(filter(x->!occursin("date",x), names(df))))
        #dout = select(df, Not(:date))
        #dout = dout[!,Cols([:YY,:MM,:HH,:DD],1:end-4)]
        dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
        #cls = propertynames(df)|>sort|>reverse
        #df = df[!,cls[2:end]] 
        CSV.write(file, dout, transform = (col, val) -> something(val, missing),delim="\t")  
        nothing
    end


    function wawrite(df::DataFrame,file::AbstractString)
        """
        newer version with copy df and switched func positions
        """
        dout = copy(df)
        #dout[!,Cols(r"date")]
        #in("date",names(dout))
        if in("year",names(dout))
            @warn "yearcol found!"
            CSV.write(file, dout, 
            transform = (col, val) -> something(val, missing), delim="\t")  
            return
        end
        dout.YY = map(x ->year(x),dout.date)
        dout.MM = map(x ->month(x),dout.date)
        dout.DD = map(x ->day(x),dout.date)
        dout[!, "HH"] .= 0
        dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
        CSV.write(file, dout, 
        transform = (col, val) -> something(val, missing), delim="\t")  
        nothing
    end


    function writedf(file, table)
        CSV.write(file, table, transform = (col, val) -> something(val, missing),delim="\t")  
        nothing
    end

    function writedesc(file, table)
        CSV.write(file, describe(table), transform = (col, val) -> something(val, missing),delim="\t")  
        nothing
    end

    #writedf("tst.csv",describe(df))
    #writedesc("tst.csv",df)

    #wc -l in julia:
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

    function vgrep(regex, file_ending)
        files = filter(file -> endswith(file, file_ending), readdir())
        # loop over each file
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    # check if the line matches the regex
                    if occursin(Regex(regex), line)
                        #m=count(_ -> true, line) #das zählt die linechars
                        println("$file: $counter:\t $line")
                    end
                end
            end
        end
    end

    function vg(snippet::AbstractString, file_ending::AbstractString)
        files = filter(file -> endswith(file, file_ending), readdir())
        # loop over each file
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    # check if the line matches the regex
                    #if occursin(Regex(regex), line)
                    if contains(line,snippet)
    #                    println("$file: $counter:\t $line")
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
    end

    function vgjl(snippet::AbstractString)
        owd=pwd()
        cd("C:/Users/Public/Documents/Python_Scripts/julia")
        files = filter(file -> endswith(file, ".jl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd("C:/Users/Public/Documents/Python_Scripts/julia/win")
        files = filter(file -> endswith(file, ".jl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(owd)
    end



    # for (root, dirs, files) in walkdir(owd)
    #     z = filter(file -> (endswith(file, ".jl")) && isfile(file), files)
    #     println(z)   
    # end

    function vgjlrec(snippet::AbstractString)
        """
        recursive grep
        """
        owd="C:/Users/Public/Documents/Python_Scripts/julia"
        for (root, dirs, files) in walkdir(owd)
            for file in files 
                if (endswith(file, ".jl"))
                    pt=(joinpath(root, file))
                    open(pt) do f
                        counter = 0 # Zähler initialisieren
                        for line in eachline(f)
                            counter += 1 # Zähler erhöhen
                            if contains(line,snippet)
                                printstyled("$counter:\t",color=:light_red) 
                                printstyled("$pt:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                                printstyled("$line\n",color=:green,bold=true) 
                            end
                        end
                    end
                end
            end
        end
    end

    #vgjl("wsl")

    function vgro(snippet::AbstractString)
        owd=pwd()
        cd("D:/Fernerkundungsdaten/Klassifikation/R-Sessions")
        files = filter(file -> endswith(file, ".R"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(owd)
    end

    function vgr(snippet::AbstractString)
        owd="D:/Fernerkundungsdaten/Klassifikation/R-Sessions"
        files = filter(file -> endswith(file, ".R"), readdir(owd,join=true))
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled(rpad("$file:",50),color=:light_magenta)
                        #printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled(lpad("$line\n",30),color=:green,bold=true)
                        #underline = true 
                    end
                end
            end
        end
    end

    function vgpy(snippet::AbstractString)
        owd="C:/Users/Public/Documents/Python_Scripts"
        files = filter(file -> endswith(file, ".py"), readdir(owd,join=true))
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled(rpad("$file:",50),color=:light_magenta)
                        #printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled(lpad("$line\n",30),color=:green,bold=true)
                        #underline = true 
                    end
                end
            end
        end
    end

    function vgctl(snippet::AbstractString)
        """
        hint. vgctl("set \$TS")
        """
        owd=pwd()
        nwd="D:/Wasim/regio/control/"
        nwd2="D:/temp/saale/control/"
        nwd3="D:/Wasim/Tanalys/DEM/brend_fab/control/"
        nwd4="D:/Wasim/regio/control/"
        nwd5="D:/Wasim/streu/control/"
        cd(nwd)
        #println("greps from *ctl from  \n$nwd and \n$nwd2...")
        println("greps from *ctl from  \n$nwd \n$nwd2 \n$nwd3 \n$nwd4 \n$nwd5...")
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd2)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd2",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd3)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd3",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd4)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd4",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd5)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd5",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(owd)
    end

    function rglob(prefix::AbstractString)
        rootdir="."
        results = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                if (occursin(Regex(prefix,"i"),filename))
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
        return results
    end

    function rglob(prefix::Regex)
        rootdir="."
        results = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if (occursin(prefix,filename))
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
        return results
    end

    #go dir up
    function cdb()
        dirname(pwd())|>cd
        pwd()|>println
    end

    function fdd()
        cwd = pwd()
        dirs = readdir(".")
        s = []
        for dir in dirs
            if isdir(dir)
                push!(s,joinpath(cwd, dir))
                size = 0
                for (root, dirs, files) in walkdir(dir)
                    for file in files
                        size += stat(joinpath(root, file)).size
                    end
                end
            @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
    #        else
    #	    @printf("%-40s\n","$(cwd)");
        end
        end
        return(s)
    end

    function vgr(regex, file_ending)
        rootdir=pwd()
        println("starting on: $rootdir...\n searching for >> $regex << with file ending >> $file_ending <<\n")
        files = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                #if (occursin(Regex(prefix,"i"),filename))
                if (endswith(filename, file_ending))
                    push!(files, joinpath(looproot, filename)) 
                end
            end
        end
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if occursin(Regex(regex,"i"), line)
                        println("$file: $counter:\t $line")
                    end
                end
            end
        end
    end

    function vjl(regex)
        # greps jl from julia folder
        pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia";
        file_ending=".jl"
        files = filter(file -> endswith(file, file_ending), readdir(pt,join=true))
        for file in files
            open(file) do f
                counter = 0
                for line in eachline(f)
                    counter += 1
                    if occursin(Regex(regex,"i"), line)
                        println("$file: $counter:\t $line")
                    end
                end
            end
        end
    end
    
    function mall(files::Vector{Any};xcol=:date)
        """reads, reduces + merges by date"""
        #files
        dfs = DataFrame[]
        for file in files
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
            file_path = file
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        push!(dfs, p1)
            end
        end
        df = reduce((left, right) -> 
        innerjoin(left, right, on = xcol,makeunique=true), 
        dfs)
        return(df)
    end

    function mall(files::Vector{DataFrame};xcol=:date)
        "reduces + merges by date"
        df = reduce((left, right) -> 
        innerjoin(left, right, on = xcol,makeunique=true), 
        files)
        return(df)
    end

    function getdf(regex::AbstractString,dfs::Vector{DataFrame})
        "selects first match..."
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),
            dfs))] |> first
        return(df)
    end

    function getdf(dfs::Vector{DataFrame},index::Integer)
        "selects by index.."
        df = getindex(dfs,index)
        return(df)
    end

    function getdf(regex::AbstractString,dfs::Vector{Any})
        "selects by regex.."
        #df = getindex(dfs,index)
        
        #y = filter(x->!occursin("date",x),names(df))
        al = filter(x->occursin(Regex(regex,"i"),x),dfs)
        
        df = filter(x->!occursin(r"txt|yrly|nc|png|svg|ftz_0|ftz",x),al)|>
        first |>readdf
            
        return(df)
    end
 
    function yrsum(x::String)
        df = readdf(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        return(df_yearsum)
    end

    function yrsum(x::Regex)
        df = readdf(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        return(df_yearsum)
    end

    function yrsum(x::DataFrame)
        df = copy(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        return(df_yearsum)
    end

    function yrmean(x::String)
        df = readdf(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
        return(df_yearsum)
    end

    function yrmean(x::Regex)
        df = readdf(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
        return(df_yearsum)
    end

    function yrmean(x::DataFrame)
        df = copy(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
        return(df_yearsum)
    end

    # function bardf(x::String)
    #     "with String"
    #     df = readdf(x)
    #     y = filter(x->!occursin("date",x),names(df))
    #     s = map(y -> Symbol(y),y)
    #     #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
    #         ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #     df[!, :year] = year.(df[!,:date]);
    #     df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
    #     @df df_yearsum Plots.plot(:year,
    #         cols(s),
    #         legend = :topright, 
    #         title=ti,
    #         seriestype=:bar)
    # end

    # function bardf(x::Regex)
    #     "with regex, and new metadata extraction"
    #     df = readdf(x)
    #     y = filter(x->!occursin("date",x),names(df))
    #     s = map(y -> Symbol(y),y)
    #     #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
    #         ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #     df[!, :year] = year.(df[!,:date]);
    #     df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
    #     @df df_yearsum Plots.plot(:year,
    #         cols(s),
    #         legend = :topright, 
    #         title=ti,
    #         seriestype=:bar)
    # end

    # function bardf(x::DataFrame)
    #     "with DataFrame input"
    #         df = x
    #         y = filter(x->!occursin("date",x),names(df))
    #         s = map(y -> Symbol(y),y)
    #             ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #         df[!, :year] = year.(df[!,:date]);
    #         df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
    #         @df df_yearsum Plots.plot(:year,
    #             cols(s),
    #             legend = :topright, 
    #             title=ti,
    #             seriestype=:bar)
    # end

    # function bardfm(x::String)
    #     "with String"
    #     df = readdf(x)
    #     y = filter(x->!occursin("date",x),names(df))
    #     s = map(y -> Symbol(y),y)
    #         ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #     df[!, :year] = year.(df[!,:date]);
    #     df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
    #     @df df_yearsum Plots.plot(:year,
    #         cols(s),
    #         legend = :topright, 
    #         title=ti,
    #         seriestype=:bar)
    # end

    # function bardfm(x::Regex)
    #     "with regex, and new metadata extraction"
    #     df = readdf(x)
    #     y = filter(x->!occursin("date",x),names(df))
    #     s = map(y -> Symbol(y),y)
    #     #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
    #         ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #     df[!, :year] = year.(df[!,:date]);
    #     df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
    #     @df df_yearsum Plots.plot(:year,
    #         cols(s),
    #         legend = :topright, 
    #         title=ti,
    #         seriestype=:bar)
    # end

    # function bardfm(x::DataFrame)
    #     "with DataFrame input"
    #         df = x
    #         y = filter(x->!occursin("date",x),names(df))
    #         s = map(y -> Symbol(y),y)
    #             ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #         df[!, :year] = year.(df[!,:date]);
    #         df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
    #         @df df_yearsum Plots.plot(:year,
    #             cols(s),
    #             legend = :topright, 
    #             title=ti,
    #             seriestype=:bar)
    # end

    # function bardfm!(x::DataFrame)
    #     "with DataFrame input"
    #         df = x
    #         y = filter(x->!occursin("date",x),names(df))
    #         s = map(y -> Symbol(y),y)
    #             ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #         df[!, :year] = year.(df[!,:date]);
    #         df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
    #         @df df_yearsum Plots.plot!(:year,
    #             cols(s),
    #             legend = :topright, 
    #             title=ti,
    #             seriestype=:bar)
    # end

    function cnt()
        return(length(readdir(pwd())))
    end

    function du()
        cwd = pwd()
        n = length(readdir(cwd))
    #    dirs = readdir(cwd)
        osize = 0
        for (root, dirs, files) in walkdir(cwd)
        for file in files
            osize += stat(joinpath(root, file)).size
        end
        end 
        println("$(n) files in directory")
        @printf("%-40s %15.2f MB\n","$(cwd):",osize/1024^2)
    end 

    function fdf(regex::AbstractString,dfs::Vector{DataFrame},f::Function)
        "selects first match and applies function..."
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
        map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
        )] |> first
        f(df)

        #like: fdf("win",dfs,yrsum) 
        #like: fdf("win",dfs,describe) 
        #indexin(1:length(dfs),
        #map(x->basename(only(DataFrames.metadata(x))[2]),dfs))
    end

    function route_from_dir(dir::String)
        dirs = readdir(dir)
        routes::Vector{String} = []
        for directory in dirs
            if isfile("$dir/" * directory)
                push!(routes, "$dir/$directory")
            else
                if ~(directory in routes)
                    newread = dir * "/$directory"
                    newrs = route_from_dir(newread)
                    [push!(routes, r) for r in newrs]
                end
            end
        end
        routes
    end

    function filterdf(dfs::Vector{Any})
        "selects presumably dfs from vector..."
        df = filter(x->!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",x),dfs)
        df = loadalldfs(df) #geht
        return(df)
    end

    function filterdf(regex::AbstractString,dfs::Vector{DataFrame})
        """
        for namestring, see as dfilter
        selects df from dfvector...
        same as getdf
        example:
        filterdf("clou",dfs)|>bardfm 
        """
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),
            dfs))] |> first
        
        # filter(n->occursin(Regex(regex,"i"),n),
        # map(x->basename(only(DataFrames.metadata(x))[2]),
        # dfs)
        # )
    end

    function llf()
        cwd = pwd()
        #n = length(readdir(cwd))
        v = []
        v_in_subdir = []
        for (root, dirs, files) in walkdir(cwd)
            for file in files
                file_path = joinpath(root, file)
                if isfile(file_path) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|tif|jpg|png|svg",file))
                    push!(v,file_path)
                elseif isdir(file_path)
                    v_in_subdir = llf()
                    println("is dir ",file_path)
                end
            end
            #return(vcat(v,v_in_subdir))
            return(v)
        end
    end


    function monsum(x::String)
        df = readf(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(groupby(df, :month), y .=> sum .=> y);
        return(df_monthsum)
    end

    function monsum(x::DataFrame)
        df = x
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(groupby(df, :month), y .=> sum .=> y);
        return(df_monthsum)
    end

    function monmean(x::String)
        df = readdf(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(groupby(df, :month), y .=> mean .=> y);
        return(df_monthsum)
    end

    function monmean(x::DataFrame)
        df = x
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(groupby(df, :month), y .=> mean .=> y);
        return(df_monthsum)
    end

    # function barp(x::DataFrame)
    #     "with DataFrame input"
    #         df = x
    #             ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #         if any(x->occursin("year",x),names(df))
    #             ln = Symbol.(filter(x->!occursin("year",x),names(df)))
    #             @df df Plots.plot(:year,
    #                 cols(ln),
    #                 legend = :topright, 
    #                 title=ti,
    #                 seriestype=:bar) #color=:lightrainbow
    #         elseif any(x->occursin("month",x),names(df))
    #             ln = Symbol.(filter(x->!occursin("month",x),names(df)))
    #             @df df Plots.plot(:month,
    #                 cols(ln),
    #                 legend = :topright, 
    #                 title=ti,
    #                 seriestype=:bar)
    #         elseif (
    #             any(x->occursin("month",x),names(df)) & 
    #             any(x->occursin("year",x),names(df))            
    #             )
    #             ln = (filter(x->!occursin("month",x),names(df)))
    #             ln = Symbol.(filter(x->!occursin("year",x),ln))
    #             @df df Plots.plot(:month,
    #                 cols(ln),
    #                 legend = :topright, 
    #                 title=ti,
    #                 seriestype=:bar)
    #         else
    #             dfp(df)        
    #         end
    # end

    function vg2(regex::AbstractString, ending::AbstractString)
        cmd = `grep --color=always -C2 -rIHn -E "$regex" --include="*.$ending"`
        run(cmd)
    end

    function so_read(x::AbstractString)
        "--- reader with drop exept of first col ---"
        ms=["-9999","lin","log"]
        df::DataFrame = CSV.read(x,DataFrame,
        missingstring=ms,
        types = Float64,
        delim="\t",
        silencewarnings=true,
        normalizenames=true,
        drop=(i, nm) -> i == 4) |> dropmissing
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        #df=df[:,Not(1:3)]
        df=df[:,Cols(4,end)]
        DataFrames.metadata!(df, "filename", x, style=:note);
    end


    function readmhm(x::AbstractString)
        "--- main reader ---"
        ms=["-9999","lin","log"]
        #x=lk
        df::DataFrame = CSV.read(x,DataFrame,
        header=false,
        missingstring = ms,
        types =  Dict(6=>Float64,
        1=>Int64,
        2=>Int64,
        3=>Int64,
        4=>Int64,
        5=>Int64), 
        skipto = 6,
        delim=" ",
        ignorerepeated=true,
        silencewarnings=false,
        normalizenames=true)
        DataFrames.metadata!(df, "filename", x, style=:note);
        #DataFrames.metadata!(df, "basename", basename(x), style=:note);
        #split(basename(x),".")[1]
        #    drop=(i, nm) -> i == 4) |> dropmissing
        df = df[!,Cols(1:3,end)]      #subset to dayres only
        newnames = ["YY","MM","DD",split(basename(x),".")[1]]
        rename!(df, newnames)
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,Not(1:3)]    
    end

    # function dfsplog(dfs::Vector{DataFrame};save="")
    #     "plots and adds"
    #     df = dfs[1]
    #     s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    #     p = @df df Plots.plot(:date,
    #             cols(s),
    #             yaxis = :log,
    #             legend = false)
    #             #legend = :bottom)
    #     for i in 2:length(dfs)
    #         nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
    #         println("adding $nm")
    #         s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
    #         @df dfs[i] Plots.plot!(:date,cols(s),
    #         label="$nm") #geht, wenn oben legend true ist.
    #         # label="$nm",
    #         # legend = false)
    #         # Plots.annotate!(0.5, 0.5, text(nm, 14))
    #     end
    #     return p
    #     if !isempty(save) 
    #         Plots.savefig(p,save*".png")
    #         printstyled("$save saved as $save*.png! \n",color=:green)
    #     end
    # end

    # function dfsp(dfs::Vector{DataFrame};save="")
    #     "plots and adds"
    #     df = dfs[1]
    #     s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    #     p = @df df Plots.plot(:date,
    #             cols(s),
    #             legend = false)    #legend = :bottom)
    #     for i in 2:length(dfs)
    #         nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
    #         println("adding $nm")
    #         s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
    #         @df dfs[i] Plots.plot!(:date,cols(s),
    #         label="$nm") #geht, wenn oben legend true ist.
    #     end
    #     return p
    #     if !isempty(save) 
    #         Plots.savefig(p,save*".png")
    #         printstyled("$save saved as $save*.png! \n",color=:green)
    #     end
    # end

    function vgrepl(snippet::AbstractString)
        """
        greps from repl_history
        """
        #file = raw"/home/ubu/.julia/logs/repl_history.jl"
        file = raw"C:\Users\chs72fw\.julia\logs\repl_history.jl"
        #files = filter(file -> endswith(file, ".py"), readdir())
        #for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        #end
        #cd(owd)
    end

    function glob(x::AbstractString)
        """
        greps from current dir iRegex
        """
        filter(file -> occursin(Regex(x,"i"),file), readdir())
    end

    function glob(x::Regex)
        """
        greps from current dir Regex
        """
        filter(file -> occursin(x,file), readdir())
    end

    function reorder_df(df::DataFrame)
        """
        date to last position
        """
        df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        return(df)
    end


    function tree(dir::AbstractString = pwd(); indent::String = "    ")
        println(dir)
        for (root, dirs, files) in walkdir(dir)
            # for file in files
            #     println(indent, "├── ", file)
            # end
            for d in dirs
                println(indent, "├── ", d)
            end
        end
    end
    
    function qall_num()
        files = glob("qgko")
        for file in files
            # Load the file into a DataFrame
            #x = "qgkofab.m6.2010"
            x = file
            try
                df = DataFrame(CSV.File(x, header=false, 
                                    delim="\t",
                                    ignorerepeated=true,
                                    types = String
                                    )) 
                                    #types=[String, String, String])                                )
                println(x)
                pattern = r"^[LIN. R]|^[LOG. R]|^CO"
                #m=match(r".*[.]",s)
                #outfile = string(m.match,"png")
                #string.(df[i,:])
                # first(eachrow(df[!,1]))
                # m=[]
                # for i in eachrow(df[!,1])
                #     k=i
                #     n=(filter(line -> occursin(pattern,line),k))
                #     push!(m,m)
                # end
                # m
                mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
                dx = df[mask, :]
                dx = permutedims(dx) |>dropmissing
                
                #mapcols!(x -> parse(Float64, x), dx)
                #mapcols(x -> x.^2, dx)

                #basins = copy(df[1,5:end])
                #AsTable(basins)
                basins = []
                for i in copy(df[1,5:end])
                    push!(basins,i)
                end
                #size(basins)
                insert!(basins, 1, "score")
                insert!(basins, 2, "timestep")
                dx[!, "basin"] = basins


                cn = (dx[1,:])
                rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
                dout = dx[3:end,:]
                for col in names(dout)[1:end-1]
                    dout[!, col] = parse.(Float64, replace.(dout[!, col], "," => ""))
                end          
                
                #mapcols!(x -> parse(Float64, x), dout)
                #dout = hcat(dx[!,Cols(r"bas")],dx[:,Not(Cols(r"bas"))])
                dout.basin = map(x -> parse(Int, x), dout.basin)
                dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
                return(dout)
            catch
                @warn("error! files that can't be loaded as a DataFrame")
                # Skip files that can't be loaded as a DataFrame
                continue
            end
        end
    end

    function qall(;recursive=false)
        if recursive
            files = rglob("qgko")
        else
            files = glob("qgko")
        end
        outdf::Vector{DataFrame} = []
        for file in files
            # Load the file into a DataFrame
            #x = "qgkofab.m6.2010"
            x = file
            try
                df = DataFrame(CSV.File(x, header=1, 
                                    delim="\t",
                                    skipto=366,
                                    ntasks=1,
                                    types = String,
                                    silencewarnings=true,
                                    ignorerepeated=true,
                                    ignoreemptyrows=true,
                                    stripwhitespace=true))
                                    
                #df = CSV.read(x,DataFrame;ntasks=1)
                println(x)
                pattern = r"^[LIN. R]|^[LOG. R]|^CO"
                mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
                dx = df[mask, :]
                dx = permutedims(dx) |>dropmissing
                
                #mapcols!(x -> parse(Float64, x), dx)
                #mapcols(x -> x.^2, dx)

                #basins = copy(df[1,5:end])
                #AsTable(basins)
                basins = []
                #for i in copy(df[1,5:end])
                for i in names(df)[5:end]
                    push!(basins,i)
                end
                #size(basins)
                insert!(basins, 1, "score")
                insert!(basins, 2, "timestep")
                dx[!, "basin"] = basins

                cn = (dx[1,:])
                rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
                dout = dx[3:end,:]
                for col in names(dout)[1:end-1]
                    dout[!, col] = parse.(Float64, replace.(dout[!, col], "," => ""))
                end          
                #mapcols!(x -> parse(Float64, x), dout)
                #dout = hcat(dx[!,Cols(r"bas")],dx[:,Not(Cols(r"bas"))])
                dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
                dout.basin = parse.(Int64,dout[!, :basin])
                push!(outdf,dout)
            catch
                @warn("error! files that can't be loaded as a DataFrame")
                # Skip files that can't be loaded as a DataFrame
                continue
            end
        end
        return(outdf)
    end

    # function theplot(x::AbstractString)
    #     df = DataFrame(CSV.File(x))
    #     nm=names(df)[end-1] #lastbefore column (qobs)
    #     ##subset DF by value (all positive vals..)
    #     df = filter(nm => x -> x > 0, df)
    #     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    #     ndf = df[!,Not(1:4)]
    #     rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
    #     overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    #     r2 = overall_pearson_r^2
    #     #nse(simulations, evaluation)
    #     nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
    #     kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
    #     #ti = "Time Series of $(uppercase(first(split(basename(x), '-'))))"
    #     ti = first(split(basename(x),"_"))
    #     #subs = "Pearson r: $(round(overall_pearson_r, digits=2))\nPearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    #     subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    #     #p = plot(title=[ti, subs], ylabel="[unit/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    #     p = Plots.plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    #     Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
    #     Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    #     Plots.annotate!(
    #     #nrow(ndf), 0.95*maximum(ndf.Observed),
    #     :bottomright,
    #     text("$subs", 10, :black, :right))
    #     return p
    # end

    # function theplot(x::DataFrame)
    #     ndf = copy(x)
    #     @info "renaming to :Simulated,:Observed,:Date !"
    #     rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
    #     overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    #     r2 = overall_pearson_r^2
    #     nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
    #     kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
    #     ti = try
    #         basename(last(collect(DataFrames.metadata(ndf)))[2])
    #     catch
    #         @warn "No basename in metadata!"
    #         raw""
    #     end 
    #     subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    #     p = Plots.plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    #     Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
    #     Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    #     Plots.annotate!(
    #     :bottomright,
    #     text("$subs", 10, :black, :right))
    #     return p
    # end

    # function theplot(x::Regex)
    #     """
    #     sim obs plot rglob regex
    #     """
    #     x = first(rglob(x))
    #     df = DataFrame(CSV.File(x))
    #     nm=names(df)[end-1] #lastbefore column (qobs)
    #     ##subset DF by value (all positive vals..)
    #     df = filter(nm => x -> x > 0, df)
    #     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    #     ndf = df[!,Not(1:4)]
    #     rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
    #     overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    #     r2 = overall_pearson_r^2
    #     #nse(simulations, evaluation)
    #     nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
    #     kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
    #     ti = first(split(basename(x),"_"))
    #     subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    #     p = Plots.plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    #     Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
    #     Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    #     Plots.annotate!(
    #     :bottomright,
    #     text("$subs", 10, :black, :right))
    #     return p
    # end

    # ftp(z::AbstractString) = theplot(first(
    #     filter(x->occursin(Regex(z,"i"),x),
    #     filter(x->endswith(x,"qoutjl"),readdir())))
    #     )

    # function ftp(df::DataFrame)
    #     if size(df)[2]!=3
    #         #throw(@warn "wrong number of columns - using dfp!")
    #         @warn "wrong number of columns - using dfp!"
    #         @warn "need :Simulated,:Observed,:Date !"
    #         display(dfp(df))
    #         return
    #     end
    #     #ndf = copy(df)
    #     # rename!(ndf,3=>"date")
    #     # @warn "last col renamed! !"
    #     #reorder
    #     #ndf = hcat(ndf[!,Not(Cols(r"date"i))],ndf[:,Cols(r"date"i)])
    #     ndf = reorder_df(df)
    #     #nm=names(ndf)[end-2] #lastbefore date column (qobs)
    #     ##subset DF by value (all positive vals..)
    #     #ndf = filter(nm => x -> x > 0, ndf)
    #     #filter([1]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), ndf)
    #     ndf = filter(:date=> x -> !any(f -> f(x), (ismissing, isnothing)), ndf)
    #     rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
    #     dropmissing!(ndf) ##hmm sketchy..
    #     overall_pearson_r = cor(ndf[!, :Simulated],ndf[!, :Observed])
    #     r2 = overall_pearson_r^2
    #     #nse(simulations, evaluation)
    #     nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
    #     kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
    #     ti = try
    #         basename(last(collect(DataFrames.metadata(ndf)))[2])
    #     catch
    #     @info "No basename in metadata!"
    #         raw""
    #     end 
    #     subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    #     p = plot(ndf[!, :Date], ndf[!, :Simulated], 
    #     title=ti, 
    #     line=:dash, color=:blue, label="Modeled",
    #     ylabel="[mm/day]", xlabel="modeled time", 
    #     yscale=:log2, legend=:topleft)
    #     plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    #     annotate!(
    #     :bottomright,
    #     text("$subs", 10, :black, :right))
    #     return p
    # end

    # function ftp(z::Regex)
    #     """
    #     first match of regex and qoutjl
    #     """
    #     theplot(first(
    #     filter(x->occursin(z,x),
    #     filter(x->endswith(x,"qoutjl"),readdir()))
    #     )
    #     )
    # end

    function grep_KGE(path::AbstractString)
        #@printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
        for file in filter(file -> endswith(file, "_output.txt"), readdir(path))
            #output = read(file, String)
            output = DelimitedFiles.readdlm(file,'\t', String)
            #match = occursin(r"KGE.*[0-9].[3-9]", output)
            match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
            if !isempty(match)
                #@printf("%s: %s\n", file,match)
                #@printf("%s:", first(split(file,"_qout")))
                fn = first(split(file,"_qout"))
                for line in match
                    #@printf("\t%s\n", line)
                    line = strip(line)  # remove leading and trailing whitespace
                    line = join(split(line), " ")  ##remove inner whitespaces
                    printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
                end
            end
        end
    end

    function kgegrepr()
        path = pwd()
        files = rglob(r"_output.txt|_outputjl") #recurse
        @printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file,'\t', String)
            match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
            if !isempty(match)
                fn = first(split(file,"_qout"))
                for line in match
                    line = strip(line)  # remove leading and trailing whitespace
                    line = join(split(line), " ")  ##remove inner whitespaces
                    printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
                end
            end
        end
    end

    function nsegrepr()
        path = pwd()
        files = rglob(r"_output.txt|_outputjl") #recurse
        @printf("Searching for NSE values > 0.3 in files matching pattern %s\n", path)
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file,'\t', String)
            match = Grep.grep(r"NSE.*[0-9].[3-9]",output)
            if !isempty(match)
                fn = first(split(file,"_qout"))
                for line in match
                    line = strip(line) 
                    line = join(split(line), " ")  ##remove inner whitespaces
                    printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
                end
            end
        end
    end

    function odfr(x::AbstractString)
        """
        --- reader with fewer constrains ---
        no |> dropmissing 
        df[!,Cols(r"^Col|date")] |>dfp  
        """
        ms=["-9999","lin","log"]
        df::DataFrame = CSV.read(x,    
        DataFrame,    
        missingstring=ms,
        #delim="\t",    
        normalizenames=true,
        types=Float64,
        drop=(i, nm) -> i == 4) 
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,4:end]
        df=df[!,Not(Cols(r"^Column"))] #drops names starting with Column, usually all missing
        DataFrames.metadata!(df, "filename", x, style=:note);
    end

    function odfr(x::Regex)
        """
        --- reader with fewer constrains ---
        with |> dropmissing on 2nd col
        df[!,Cols(r"^Col|date")] |>dfp  
        """
        x=globdf(x)|>first
        println("reading $x ...")
        ms=["-9999","lin","log","--","A-z"]
        df = CSV.read(x,    
        DataFrame,    
        missingstring=ms,
        #delim="\t",    
        types = Float64,
        normalizenames=true,
        drop=(i, nm) -> i == 4) 
        dropmissing!(df , 2) #2nd column
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,4:end]
        df=df[!,Not(Cols(r"^Column"))] #drops names starting with Column, usually all missing
        #renamer
        for x in 1:size(df,2)-2
            rename!(df,x=>"C"*names(df)[x])
        end
        DataFrames.metadata!(df, "filename", x, style=:note);
    end

    # function readalloutput()
    #     """
    #     reads timeseries and stores to Vector{DataFrame}
    #     reads NetCDFs and stores to Vector{Any}

    #     usage: 
    #     dfs,ncs = readalloutput()

    #     errors if not rmeq()


    #     """
    #     cwd = "."
    #     dfs=loadalldfs(cwd)
    #     ncs=readallras(cwd)
    #     return(dfs,ncs)
    # end

    # loadalloutput = readalloutput


    function gofbatch()
        println("batch R Script for GOF")
        arr = filter(x -> isfile(x) && endswith(x, "_qout") && !occursin(r"\.(png|svg|txt|html|ftz|ftz_0|list|nc|xml|sh|grd|yrly)$", x), readdir())
        for i in arr
            println(i)
        end
        if length(arr) < 1
            println("No match!\nneed qout files")
            return 1
        end
        for i in arr
            x = basename(i)
            run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof3.R" $x`)
        end
    end

    function grep_KGE(path::AbstractString)
        #@printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
        try
            for file in filter(file -> endswith(file, "_output.txt"), readdir(path))
                #output = read(file, String)
                output = DelimitedFiles.readdlm(file,'\t', String)
                #match = occursin(r"KGE.*[0-9].[3-9]", output)
                match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
                if !isempty(match)
                    #@printf("%s: %s\n", file,match)
                    #@printf("%s:", first(split(file,"_qout")))
                    fn = first(split(file,"_qout"))
                    for line in match
                        #@printf("\t%s\n", line)
                        line = strip(line)  # remove leading and trailing whitespace
                        line = join(split(line), " ")  ##remove inner whitespaces
                        printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
                    end
                end
            end    
        catch
            @error "no files present!"
            return
        end
        
    end

    function globdf(x::Regex)
        """
        greps ts from current dir Regex
        """
        filter(file -> (occursin(x,file) & 
        (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
        ), readdir())
    end

    function globdf(x::AbstractString)
        """
        greps ts from current dir Regex
        """
        # Use a regex to match x and exclude the unwanted extensions
        #regex = Regex(x * "(?!\\.(xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg))","i")
        #regex = Regex(x * "(?!\\.(xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg))")
        filter(file -> (occursin(
            Regex(x,"i"),file) & 
        (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
        ), readdir())
        #regex = Regex("\\Q$x*?!\\.(xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg)","i")
        # Use a list comprehension to filter the files
        #files = [file for file in readdir() if occursin(regex, file)]
    end

    # function ldfpall(x::Regex)
    #     "reads, reduces + merges by date and plots log y-axis"
    #     files = rglob(x)
    #     dfs = DataFrame[]
    #     for file in files
    #         if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
    #         file_path = file
    #     println("reading ",file_path,"...")
    #     p1 = readdf(file_path)
    #     ##renamer
    #     for x in 1:size(p1,2)-1
    #         rename!(p1,x=>basename(file_path)*names(p1)[x])
    #     end
    #     push!(dfs, p1)
    #         end
    #     end
    #     df = reduce((left, right) -> 
    #     innerjoin(left, right, on = :date,makeunique=true),dfs)
    #     ##to preserve column order and names (date at last position)
    #     #df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
    #     y = filter(x->!occursin("date",x), names(df))
    #     s = map(y -> Symbol(y),y)
    #     @df df Plots.plot(:date,
    #             cols(s), yaxis = :log,
    #             legend = :outertopright)
    # end

    # function dfl(x::DataFrame)
    #     "selects first match and plots in log y-axis..."
    #     df=x
    #     ti = try
    #         DataFrames.metadata(df)|>only|>last|>basename
    #     catch
    #     @warn "No basename in metadata!"
    #     ti = raw"" 
    #     end
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
    #         @df df Plots.plot(:year,cols(s),yaxis=:log,legend = :topright, title=ti)
    #     else   
    #         s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #         @df df Plots.plot(:date,cols(s),yaxis=:log10, 
    #         legend = :topright, title=ti)
    #         end
    # end

    
    function qbb()
        """
        return a vector of DFs
        """
        files = rglob("qgko")
        dfs = []
        for file in files
            x = file
            try
                df = DataFrame(CSV.File(x, header=false, 
                                    delim="\t",
                                    silencewarnings=true,
                                    ignorerepeated=true,
                                    types = String))
                                    #types=[String, String, String])                                )
                println(x)
                pattern = r"^[LIN. R]|^[LOG. R]|^CO"
                mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
                dx = df[mask, :]
                dx = permutedims(dx) |>dropmissing
                basins = []
                for i in copy(df[1,5:end])
                    push!(basins,string.("B_"*i))
                end
                insert!(basins, 1, "score")
                insert!(basins, 2, "timestep")
                dx[!, "basin"] = basins
                cn = (dx[1,:])
                rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
                dout = dx[3:end,:]
                dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
                DataFrames.metadata!(dout, "filename", file, style=:note);
                push!(dfs,dout)
                
            catch
                @warn("error! ")
                # Skip files that can't be loaded as a DataFrame
                continue
            end
        end
        return(dfs)
    end

    function readf(x::String)
        """
        Read the text file, preserve line 1 as header column
        """
        #ms=["-9999","lin","log","--"] #missingstring=ms,
        df = CSV.read(x, DataFrame, delim="\t",header=1,
        silencewarnings=true,
        types=Float64)
        # filter numeric only / subsetting on first column
        #df = filter([1]=> x -> !any(f -> f(x), (ismissing, isnothing)), df)
        dropmissing!(df,1) #better than above
        #df = filter([6]=> x -> !any(-9999.0,x), df)
        # replace to missing...inplace doesnt work!
        #df = filter([5]=> x -> !any(f -> f(x),replace(-9999.0 => missing)), df)
        for i in 5:size(df,2)
            #replace!(df[!,i],-9999.0 => missing)
            df[!,i]=replace(df[!,i],-9999.0 => missing)
        end 
        for i in 5:size(df,2)
            replace!(df[!,i],-9999.0 => missing)
        end
        # map to int for dates
        for i in 1:3
            df[!,i]=map(x ->Int(x),df[!,i])
        end
        #and parse dates...
        df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
        df=df[:,Not(1:4)]
        metadata!(df, "filename", x, style=:note);
    end

    function routeg(input_file::String, output_file::String)
        open(output_file, "w") do output
            line_num = 0
            in_range = false

            for line in eachline(input_file)
                line_num += 1

                if line_num > 50 && contains(line, "routing_model")
                    in_range = true
                    line = replace(line, r"ß" => "ss")
                    # line = replace(line, r"[\/]" => "_")
                    # line = replace(line, r"_" => "-")
                    line = replace(line, r"[,,]" => "")
                    line = replace(line, r"\xc4" => "Ae")
                    line = replace(line, r"\xd6" => "Oe")
                    line = replace(line, r"\xdc" => "Ue")
                    line = replace(line, r"\xe4" => "ae")
                    line = replace(line, r"\xf6" => "oe")
                    line = replace(line, r"\xfc" => "ue")
                    line = replace(line, r"\xdf" => "ss")
                    println(output, line)
                elseif in_range && contains(line, "timeoffset")
                    in_range = false
                    line = replace(line, r"ß" => "ss")
                    # line = replace(line, r"[\/]" => "_")
                    # line = replace(line, r"_" => "-")
                    line = replace(line, r"[,,]" => "")
                    line = replace(line, r"\xc4" => "Ae")
                    line = replace(line, r"\xd6" => "Oe")
                    line = replace(line, r"\xdc" => "Ue")
                    line = replace(line, r"\xe4" => "ae")
                    line = replace(line, r"\xf6" => "oe")
                    line = replace(line, r"\xfc" => "ue")
                    line = replace(line, r"\xdf" => "ss")
                    println(output, line)
                elseif in_range
                    line = replace(line, r"ß" => "ss")
                    # line = replace(line, r"[\/]" => "_")
                    # line = replace(line, r"_" => "-")
                    line = replace(line, r"[,,]" => "")
                    line = replace(line, r"\xc4" => "Ae")
                    line = replace(line, r"\xd6" => "Oe")
                    line = replace(line, r"\xdc" => "Ue")
                    line = replace(line, r"\xe4" => "ae")
                    line = replace(line, r"\xf6" => "oe")
                    line = replace(line, r"\xfc" => "ue")
                    line = replace(line, r"\xdf" => "ss")
                    println(output, line)
                end
            end
        end
    end

    
    function process_file2(input_file::String)
        """
        removes " inside textfile
        """
        lines = readlines(input_file)
        output_file = input_file * ".tmp"
        
        open(output_file, "w") do file
            for line in lines
                line = replace(line, "\"" => "")
                if length(line) > 0
                    println(file, line)
                end
            end
        end
        
        # Rename the temporary file to replace the original file
        mv(output_file, input_file; force=true)
    end

    function ggofbatch()
        println("batch R Script for GOF")
        arr = filter(x -> isfile(x) && endswith(x, "qoutjl") && !occursin(r"\.(png|svg|txt|html|ftz|ftz_0|list|nc|xml|sh|grd|yrly)$", x), readdir())
        for i in arr
            println(i)
        end
        if length(arr) < 1
            println("No match!\nneed qout files")
            return 1
        end
        for i in arr
            x = basename(i)
            run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof4.R" $x`)
        end
    end

    function te(x::String)
        df = dfr(x)
        dy = yrsum(df)
        for row in eachrow(dy)
            pot=row[2]
            real=row[3]
            year=row[1]
            ref=row[end]
            ref=round((ref/365)*100; digits=2)
        #printstyled("Tdiff for 100 days in $year is $ref [mm]\t|\t",color=:magenta)
        printstyled(rpad("Tdiff for 100 days in $year is $ref [mm]", 45)        
        ,color=:magenta)
        if ref <= 0
            println("\n")
        elseif ref <= 5
            println("conditions are very moist")
        elseif ref <= 10
            println("conditions are moist")
        elseif ref <= 15
            println("conditions are rather moist")
        elseif ref <= 20
            println("conditions are quite moist")
        elseif ref <= 30
            println("conditions are quite dry")
        elseif ref <= 40
            println("conditions are rather dry")
        elseif ref <= 50
            println("conditions are dry")
        elseif ref <= 70
            println("conditions are very dry")
        elseif ref <= 7e10
            println("conditions are exceptionally dry")
        end
    end
    end

    function te(td::DataFrame)
        df = copy(td)
        y = filter(x->!occursin(r"date",x),names(df))
        df[!, :year] = year.(df[!,:date]);
        df = df[!,Not(:date)]
        dy = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        for row in eachrow(dy)
            pot=row[2]
            real=row[3]
            year=row[1]
            ref=row[end]
            ref=round((ref/365)*100; digits=2)
        printstyled(rpad("Tdiff for 100 days in $year is $ref [mm]", 45)        
        ,color=:magenta)
        if ref <= 0
            println("\n")
        elseif ref <= 5
            println("conditions are very moist")
        elseif ref <= 10
            println("conditions are moist")
        elseif ref <= 15
            println("conditions are rather moist")
        elseif ref <= 20
            println("conditions are quite moist")
        elseif ref <= 30
            println("conditions are quite dry")
        elseif ref <= 40
            println("conditions are rather dry")
        elseif ref <= 50
            println("conditions are dry")
        elseif ref <= 70
            println("conditions are very dry")
        elseif ref <= 7e10
            println("conditions are exceptionally dry")
        end
    end
    end

    function tdiff()
        npot = glob("so_pot_trans")|>first|>dfread
        nreal =  glob("so_real_trans*")|>first|>dfread
        td = innerjoin(npot,nreal,on=:date)
        td = reorder_df(td)
        #Tdiff = hcat(td,td[!,1] .- td[!,2])
        td.Tdiff = td[!,1] .- td[!,2]
        te(td)
        return(reorder_df(td))
        #wawrite(td,"tdiff-jl.txt")
    end


    function merge_vectors(vector1::Vector{Any}, vector2::Vector{Pair{String, DateTime}})
        merged_vector = []
        for (item1, item2) in zip(vector1, vector2)
            merged_item = (item1, item2.first, item2.second)
            push!(merged_vector, merged_item)
        end
        return merged_vector
    end

    function merge_vectors(vector1::Vector{String}, vector2::Dict{String, DateTime})
        merged_vector = []
        for (item1, item2) in zip(vector1, vector2)
            merged_item = (item1, item2.first, item2.second)
            push!(merged_vector, merged_item)
        end
        return merged_vector
    end

    function latx()
        """
        list_files_sorted_by_last_change
        """
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
    
    function globf(x::Regex)
        """
        greps first match current dir Regex
        """
        first(filter(file -> occursin(x,file), readdir()))
    end
    
    function globf(x::AbstractString)
        """
        greps first match current dir Regex
        """
        first(filter(file -> occursin(Regex(x,"i"),file), readdir()))
    end
    
    function npp(fl::String)
        opener="c:/Program Files (x86)/Notepad++/notepad++.exe"
        run(`$opener $fl`)
    end


    readall = loadalldfs
    dfread = readf
    readmeteo = waread
    loaddf = readdf

    function rename_duplicates(df::DataFrame)
        # loop over the column names
        for (i, name) in enumerate(names(df))
            # check if the name contains any repeated words
            if occursin(r"(\w+)[ _]\1", name)
                # replace the duplicates with a single word
                new_name = replace(name, r"(\w+)[ _]\1" => s"\1")
                # rename the column
                rename!(df, name => new_name)
            end
        end
        # return the modified dataframe
        return df
    end


    function tff2(x::Vector{String})
        for filename in x
            b = Dict{String, Float64}()
            m = Dict{String, Float64}()
            h = Dict{String, Float64}()
            cnte = Dict{String, Int64}()
            ncols = 0

            open(filename) do file
                first_line = true
                

                for (i, line) in enumerate(eachline(file))
                    fields = split(line)
                    if first_line
                        println("filename: $filename")
                        ncols = length(fields)
                        println("no of fields: $ncols")
                        println("year\t$(join(fields[5:end], "\t"))")
                        first_line = false
                        continue
                    end

                    if match(r"^\d{4}$", fields[1]) != nothing
                        cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
                        b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
                        m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
                        h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
                    end
                end
            end

            if ncols <= 5
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols == 6
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols >= 7
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", m[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", m[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            end
        end
    end


    function ssup()
        include("C:/Users/Public/Documents/Python_Scripts/julia/win/smallfuncs.jl")
    end

    # function filterplot(regex::AbstractString,ncs::Vector{Raster})
    #     "selects first match and plots..."
    #     r = ncs[map(n->occursin(Regex(regex,"i"),n),
    #     map(x->string.(name(x)),ncs)
    #     )] |> first
    #     plot(trim(r),xlab="",ylab="",title=name(r))
    #     cz=r.metadata.val|>collect|>last|>last
    #     Plots.annotate!(:bottomright,
    #     text("cellsize: $cz", 8, :black, :right))
    # end

    # function filterplot!(regex::AbstractString,ncs::Vector{Raster})
    #     "selects first match and add to plot..."
    #     r = ncs[map(n->occursin(Regex(regex,"i"),n),
    #     map(x->string.(name(x)),ncs)
    #     )] |> first
    #     plot!(trim(r),xlab="",ylab="",title=name(r))
    # end

 
    function kge_fread()
        kge_read(pwd(),"outjl");
    end

    # function baryrsum(df::Regex)
    #     """
    #     automatically sums if only datecolumn is available
    #     """
    #     #getm(waread)
    #     df = waread(df)
    #     v = map(
    #         (x->occursin(r"date", x) & !occursin(r"year", x)),
    #         (names(df))
    #         )

    #     if any(v)
    #         df = yrsum(df) 
    #     end
        
    #     s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
    #     ti = try
    #         z=DataFrames.metadata(df)|>only|>last|>basename
    #         basename(pwd())*" $z"
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = "Series of "*basename(pwd())
    #     end
    #     @df df groupedbar(df.year,cols(s), 
    #     legend = :outertopright,
    #     xticks = df.year,
    #     xrotation = 45,
    #     xlabel = "", ylabel = "[mm]", title = ti)
    # end

    # function baryrsum(df::DataFrame)
    #     """
    #     automatically sums if only datecolumn is available
    #     """
    #     v = map(
    #         (x->occursin(r"date", x) & !occursin(r"year", x)),
    #         (names(df))
    #         )

    #     #v = map(x->occursin(r"date", x),(names(df)))
        
    #     if any(v)
    #         df = yrsum(df) 
    #         # inplace yrsum:
    #         # y = filter(x->!occursin("date",x),names(df))
    #         # df[!, :year] = year.(df[!,:date])
    #         # df = DataFrames.combine(groupby(df, :year), 
    #         #     y .=> sum .=> y)
    #     end
        
    #     s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
    #     ti = try
    #         z=DataFrames.metadata(df)|>only|>last|>basename
    #         basename(pwd())*" $z"
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = "Series of "*basename(pwd())
    #     end
    #     @df df groupedbar(df.year,cols(s), 
    #     legend = :outertopright,
    #     xticks = df.year,
    #     xrotation = 45,
    #     xlabel = "", ylabel = "[mm]", title = ti)
    # end

    # function baryrmean(df::DataFrame)
    #     """
    #     automatically sums if only datecolumn is available
    #     """
    #     v = map(
    #         (x->occursin(r"date", x) & !occursin(r"year", x)),
    #         (names(df))
    #         )

    #     if any(v)
    #         # #df = yrsum(df), but inplace
    #         # y = filter(x->!occursin("date",x),names(df))
    #         # df[!, :year] = year.(df[!,:date])
    #         # df = DataFrames.combine(groupby(df, :year), 
    #         #     y .=> sum .=> y)
    #         df = yrmean(df)
    #     end
        
    #     s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
    #     ti = try
    #         z=DataFrames.metadata(df)|>only|>last|>basename
    #         basename(pwd())*" $z"
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = "Series of "*basename(pwd())
    #     end
    #     @df df groupedbar(df.year,cols(s), 
    #     legend = :outertopright,
    #     xticks = df.year,
    #     xrotation = 45,
    #     xlabel = "", ylabel = "[mm]", title = ti)
    # end

    # function baryr(df::DataFrame)
    #     s = filter(x -> !(occursin(r"year|date", x)), names(df))
    #     for x in s
    #         newname = replace(x, "_1" => "")
    #         rename!(df, Dict(x => newname))
    #     end
    #     s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
    #     # Create a grouped bar plot for each column    
    #     # StatsPlots.groupedbar(df.year,[df[!, col] for col in s], 
    #     # #group = s,#repeat(s, outer = size(df, 1)),
    #     # xlabel = "Year", ylabel = "Value", title = "Grouped Bar Plot")
    #     ti = try
    #         z=DataFrames.metadata(df)|>only|>last|>basename
    #         basename(pwd())*" $z"
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = "Series of "*basename(pwd())
    #     end
    #     @df df groupedbar(df.year,cols(s), 
    #     #group = s,#repeat(s, outer = size(df, 1)),
    #     legend = :outertopright,
    #     xticks = df.year,
    #     xrotation = 45,
    #     xlabel = "", ylabel = "[mm]", title = ti)
    # end

    function readf(x::Regex)
        """
        Read the text file, preserve line 1 as header column
        """
        x=glob(x)|>first
        df = CSV.read(x, DataFrame, delim="\t",header=1,
        silencewarnings=true,types=Float64)
        dropmissing!(df,1)
        for i in 5:size(df,2)
            df[!,i]=replace(df[!,i],-9999.0 => missing)
        end 
        # map to int for dates
        for i in 1:3
            df[!,i]=map(x ->Int(x),df[!,i])
        end
        #and parse dates...
        df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
        df=df[:,Not(1:4)]
        metadata!(df, "filename", x, style=:note);
    end

    function waread2(x::Regex)
        """
        waread2 on regex
        Read the text file, preserve line 1 as header column
        Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
        This avoids reading the entire file into memory at once, 
        which can be more memory-efficient for large datasets.
        """
        #x=glob(x)|>first
        #glob==filter(file -> occursin(Regex(x,"i"),file), readdir())

        #@time filter(file -> occursin(Regex(x,"i"),file), readdir()) #slower
        #@time filter(file -> occursin(r"sb",file), readdir()) #slower
        
        
        # @time Grep.grep(r"sb",readdir()) #0.000541 
        # x="sb"
        # @time glob(x) #faster
        # @time Grep.grep(Regex(x),readdir()) #equally fast.

        #filter(x->!occursin(r"yrly|nc|png|svg|grd",x),readdir("."))
        #@time x = first(filter(x->!occursin(r"yrly|nc|png|svg|grd",x),glob(x)))
        #slightly faster.
        #@time x = first(filter(x->!occursin(r"yrly|nc|png|svg|grd",x),Grep.grep(x,readdir())))
        inF = first(filter(x->!occursin(r"yrly|nc|png|svg|grd",x),Grep.grep(x,readdir())))

        ms = ["-9999", "lin", "log", "--"]
        df = CSV.File(inF; delim="\t", header=1, normalizenames=true, missingstring=ms, types=Float64) |> DataFrame
        dropmissing!(df,1)
        dt2 = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
        select!(df, Not(1:4))
        df.date = dt2
        metadata!(df, "filename", x, style=:note)
        return df
    end

    # function kge_df3()
    #     """
    #     should be non-recursive
    #     """
    #     x1=r"qoutjl"
    #     files = filter(file -> occursin(x1,file),
    #         readdir()[
    #             broadcast(x->!endswith(x,
    #             r"nc|pl|txt|svg|png|jpg|grd|ftz|ftz_0|list|xml|sh|yrly"), 
    #             readdir())])
    #     v = []
    #     for file_path in files
    #         if isfile(file_path)
    #             df = waread(file_path);
    #             dropmissing!(df)
    #             simulated = df[:,1]
    #             observed  = df[:,2]
    #             kge_value = kge1(simulated,observed)
    #             nse_value = nse(simulated,observed)
    #             nm = basename(file_path)
    #             println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
    #             printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
    #             push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
    #             v = DataFrame(v)
    #         end
    #     end
    #     return(v)
    # end

    # function dsbar(ds::DataFrame)
    #     ds.name=map(x->replace(x,r"-qoutjl*" => ""),ds.name)
    #     ann = map(x->string.(round(x;sigdigits=1)),ds.KGE)
    #     bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
    #         title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
    #         fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
    #         annotations = (ds.name,ds.KGE, ann, :top),
    #         bar_width = 0.6)
    # end

    function vef(obs, sim)
    return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
    end

    function ggofjl()
        """
        with lapply on all qoutjl files...
        """
        println("batch R Script for GOF")
            run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof5.R"`)
    end

    function ggofjl_nosvg()
        """
        with lapply on all qoutjl files...
        """
        println("batch R Script for GOF")
            run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof6.R"`)
    end

    # function dpr(x::Regex)
    #     """
    #     correlation plots on dataframe
    #     """
    #     df = globdf(x)|>first|>waread
    #     df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, 
    #     xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    #     r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    #     # annotate!(last(df.date), 0.85*maximum(df[!,1]),
    #     # text("R² = $r2", 10, :black, :right))
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         #:topright,
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    # function dpr!(x::Regex)
    #     """
    #     correlation plots on dataframe
    #     """
    #     df = globdf(x)|>first|>waread
    #     df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, 
    #     xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    #     r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    #     #annotate!(last(df.date), 0.85*maximum(df[!,1]),
    #     #text("R² = $r2", 10, :black, :right))
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         #:topright,
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    # function dpr(x::String)
    #     df=readdf(x)
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    #     r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    #     # annotate!(last(df.date), 0.85*maximum(df[!,1]),
    #     # text("R² = $r2", 10, :black, :right))
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         #:topright,
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    # function dpr!(x::String)
    #     df=readdf(x)
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    #     r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    #     # annotate!(last(df.date), 0.85*maximum(df[!,1]),
    #     # text("R² = $r2", 10, :black, :right))
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         #:topright,
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    # function dpr(x::DataFrame)
    #     """
    #     correlation plots on dataframe
    #     """
    #     df = copy(x)
    #     if any(map(x->occursin("year",x),names(df)))
    #         df = df[!,Not(:year)]
    #     end
        
    #     if any(map(x->occursin("month",x),names(df)))
    #         df = df[!,Not(:month)]
    #     end
        
    #     if propertynames(df)[end]!=:date
    #         df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
    #     end
        
    #     dropmissing!(df)
        
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, 
    #     xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    #     dropmissing!(df)
    #     r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    #     # annotate!(last(df.date), 0.85*maximum(df[!,1]),
    #     # text("R² = $r2", 10, :black, :right))
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         #:topright,
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    # function dpr!(x::DataFrame)
    #     """
    #     correlation plots on dataframe
    #     """
    #     df = copy(x)
    #     if any(map(x->occursin("year",x),names(df)))
    #         df = df[!,Not(:year)]
    #     end
        
    #     if any(map(x->occursin("month",x),names(df)))
    #         df = df[!,Not(:month)]
    #     end
        
    #     if propertynames(df)[end]!=:date
    #         df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
    #     end
        
    #     dropmissing!(df)
        
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    #     r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    #     # annotate!(last(df.date), 0.85*maximum(df[!,1]),
    #     # text("R² = $r2", 10, :black, :right))
    #     dropmissing!(df)
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         #:topright,
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    
    # function dpr(a::Regex,b::Regex)
    #     """
    #     correlation plots on dataframe
    #     """
    #     a = waread(a)
    #     b = waread(b)
    #     # colA = ncol(a)-1
    #     # colB = ncol(b)-1

    #     # a = a[!,Cols(colA,:date)]
    #     # b = b[!,Cols(colB,:date)]
        
    #     a = a[!,Cols(1,:date)]
    #     b = b[!,Cols(1,:date)] 

    #     df = mall(a,b)
    #     dropmissing!(df)
        
    #     df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)

    #     Plots.plot(df.date,[df[!,1], df[!,2]],  label=a, 
    #     #    seriestype = :bar,
    #         xlabel="Date", ylabel="[mm/day]",
    #         legend = :topleft)

    #     r2 = round(cor(df[!,1], df[!,2])^2, digits=2)
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)

    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)

    #     annotate!(
    #         #:topright,
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end


    function jlt(x::Vector{String})
        """
        tff2
        """
        for filename in x
            b = Dict{String, Float64}()
            m = Dict{String, Float64}()
            h = Dict{String, Float64}()
            cnte = Dict{String, Int64}()
            ncols = 0

            open(filename) do file
                first_line = true
                

                for (i, line) in enumerate(eachline(file))
                    fields = split(line)
                    if first_line
                        #println("filename: $filename")
                        printstyled("filename: $filename\n",color=:yellow)
                        ncols = length(fields)
                        println("no of fields: $ncols")
                        println("year\t$(join(fields[5:end], "\t"))")
                        first_line = false
                        continue
                    end

                    if match(r"^\d{4}$", fields[1]) != nothing
                        cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
                        b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
                        m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
                        h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
                    end
                end
            end

            if ncols <= 5
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols == 6
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols >= 7
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", m[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", m[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            end
        end
    end


    function jldf(x::Vector{String})
        result = Vector{DataFrame}()

        for filename in x
            b = Dict{String, Float64}()
            m = Dict{String, Float64}()
            h = Dict{String, Float64}()
            cnte = Dict{String, Int64}()
            ncols = 0

            open(filename) do file
                first_line = true

                for (i, line) in enumerate(eachline(file))
                    fields = split(line)
                    if first_line
                        ncols = length(fields)
                        first_line = false
                        continue
                    end

                    if match(r"^\d{4}$", fields[1]) != nothing
                        cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
                        b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
                        m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
                        h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
                    end
                end
            end

            df = DataFrame(
                Year = String[],
                B = Float64[],
                M = Float64[],
                H = Float64[],
                B_Mean = Float64[],
                M_Mean = Float64[],
                H_Mean = Float64[],
                Counts = Int64[]
            )
            metadata!(df, "filename", filename, style=:note);

            for key in sort(collect(keys(b)))
                push!(df, (
                    Year = key,
                    B = b[key],
                    M = m[key],
                    H = h[key],
                    B_Mean = b[key] / cnte[key],
                    M_Mean = m[key] / cnte[key],
                    H_Mean = h[key] / cnte[key],
                    Counts = cnte[key]
                ))
            end

            push!(result, df)
        end

        return result
    end


    function malldf(files::Vector{DataFrame},on::Symbol)
        "reduces + merges by date"
        df = reduce((left, right) -> 
        innerjoin(left, right, on = on,makeunique=true), 
        files)
        return(df)
    end

    function jldfnm(x::Vector{String})
        """
        with names
        """
        result = Vector{DataFrame}()

        for filename in x
            b = Dict{String, Float64}()
            m = Dict{String, Float64}()
            h = Dict{String, Float64}()
            cnte = Dict{String, Int64}()
            ncols = 0

            open(filename) do file
                first_line = true

                for (i, line) in enumerate(eachline(file))
                    fields = split(line)
                    if first_line
                        ncols = length(fields)
                        first_line = false
                        continue
                    end

                    if match(r"^\d{4}$", fields[1]) != nothing
                        cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
                        b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
                        m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
                        h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
                    end
                end
            end

            df = DataFrame(
                Year = String[],
                B = Float64[],
                M = Float64[],
                H = Float64[],
                B_Mean = Float64[],
                M_Mean = Float64[],
                H_Mean = Float64[],
                Counts = Int64[],
                Name = String[]
            )
            metadata!(df, "filename", filename, style=:note);
            for key in sort(collect(keys(b)))
                push!(df, (
                    Year = key,
                    B = b[key],
                    M = m[key],
                    H = h[key],
                    B_Mean = b[key] / cnte[key],
                    M_Mean = m[key] / cnte[key],
                    H_Mean = h[key] / cnte[key],
                    Counts = cnte[key],
                    Name = filename
                ))
            end
            push!(result, df)
        end

        return result
    end

    function ctl()
        # Print a message
        println("looks for control file in all xmls")
        # Loop through the current directory and its subdirectories
        matches::Vector{Any} = []
        for (root, dirs, files) in walkdir(".")
        # Loop through each file name
        for file in files
            # If the file name ends with .xml
            if endswith(file, ".xml")
            # Join the root and file name to get the full path
            path = joinpath(root, file)
            # Open the file for reading
            open(path) do f
                # Loop through each line of the file
                for line in eachline(f)
                # If the line contains 'compiling symbols in control file '
                if occursin("compiling symbols in control file ", line)
                    # Split the line by whitespace and get the fields from index 9 to 15
                    fields = split(line)[8:end] #," "
                    # Join the fields by space and print them
                    println(join(fields, " "))
                    out=join(fields, " ")
                    push!(matches,out)
                end
                end
            end
            end
        end
        end
        return(matches)
    end


    macro bash_str(s) open(`bash`,"w",stdout) do io; print(io, s); end;end
    #bash""" which python """

    macro pwrs_str(s) open(`powershell -noprofile`,"w",stdout) do io; print(io, s); end;end
    #pwrs""" which python """
    #pwrs""" pwd """
    function op()
        #pwrs""" explorer . """
        open(`powershell -noprofile explorer . `,"w",stdout)
    end 

    macro pwp_str(s) open(`powershell`,"w",stdout) do io; print(io, s); end;end

    macro cmd_str(s) open(`cmd \c`,"w",stdout) do io; print(io, s); end;end

    function rsqgrep()
        path = pwd()
        files = glob(r"_output.txt|_outputjl")
        @printf("Searching for R² values > 0.3 in files matching pattern %s\n", path)
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file, '\t', String)
            match = Grep.grep(r"R2.*[0-9].[3-9]", output)
            if !isempty(match)
                fn = first(split(file, "_qout"))
                for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                    line = strip(line)
                    line = join(split(line), " ")
                    printstyled(rpad("$fn:", 30), lpad("$line\n", 10), color = :green)
                end
            end
        end
    end

    function rsq(x::AbstractVector, y::AbstractVector)
        cor(x, y)^2
    end

    function kgegrep()
        path = pwd()
        files = glob(r"_output.txt|_outputjl") #non-recurse
        #@printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file, '\t', String)
            #match = Grep.grep(r"KGE.*[0-9].[3-9]", output)
            match = Grep.grep(r"KGE", output)
            if !isempty(match)
                fn = first(split(file, "_qout"))
                for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                    line = strip(line)  # remove leading and trailing whitespace
                    line = join(split(line), " ")  ##remove inner whitespaces
                    printstyled(rpad("$fn:", 45), lpad("$line\n", 10), color = :green)
                end
            end
        end
    end


    function nsegrep()
        path = pwd()
        files = glob(r"_output.txt|_outputjl") #non-recurse
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file, '\t', String)
            match = Grep.grep(r"mNSE", output)
            if !isempty(match)
                fn = first(split(file, "_qout"))
                for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                    line = strip(line)  # remove leading and trailing whitespace
                    line = join(split(line), " ")  ##remove inner whitespaces
                    printstyled(rpad("$fn:", 45), lpad("$line\n", 10), color = :green)
                end
            end
        end
    end

    # function baryrsum(df::Regex)
    #     """
    #     automatically sums if only datecolumn is available
    #     """
    #     df = globdf(df)|>first|>readf
    #     v = map(
    #         (x->occursin(r"date", x) & !occursin(r"year", x)),
    #         (names(df))
    #         )

    #     if any(v)
    #         df = yrsum(df) 
    #     end
        
    #     s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
    #     ti = try
    #         z=DataFrames.metadata(df)|>only|>last|>basename
    #         basename(pwd())*" $z"
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = "Series of "*basename(pwd())
    #     end
    #     @df df groupedbar(df.year,cols(s), 
    #     legend = :outertopright,
    #     xticks = df.year,
    #     xrotation = 45,
    #     xlabel = "", ylabel = "[mm]", title = ti)
    # end

    # function baryrmean(df::Regex)
    #     """
    #     automatically sums if only datecolumn is available
    #     """
    #     df = globdf(df)|>first|>readf
    #     v = map(
    #         (x->occursin(r"date", x) & !occursin(r"year", x)),
    #         (names(df))
    #         )

    #     if any(v)
    #         df = yrmean(df)
    #     end
        
    #     s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
    #     ti = try
    #         z=DataFrames.metadata(df)|>only|>last|>basename
    #         basename(pwd())*" $z"
    #     catch
    #         @warn "No basename in metadata!"
    #         ti = "Series of "*basename(pwd())
    #     end
    #     @df df groupedbar(df.year,cols(s), 
    #     legend = :outertopright,
    #     xticks = df.year,
    #     xrotation = 45,
    #     xlabel = "", ylabel = "[mm]", title = ti)
    # end

    # function kgeval()
    #     """
    #     kge barplot 
    #     """
    #     ds = kge_df3()
    #     ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.name)
    #     ann = map(x->string.(round(x;sigdigits=2)),ds.KGE)
    #     p1 = Plots.bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
    #     title = splitpath(pwd())|>last, 
    #     xrotation = 45, 
    #     fmt = :png, 
    #     size = (800, 600), 
    #     fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
    #     #annotations = (ds.name,ds.KGE, ann, :top),
    #     xaxis = "",
    #     left_margin = 10mm,
    #     bottom_margin = 15mm, 
    #     bar_width = 0.6);

    #     for i in 1:length(ds.name)
    #         Plots.annotate!(ds.name[i],ds.KGE[i],(ann[i],11,
    #             :center,:top,:black))
    #         #println(ann[i]*" added")
    #     end
    #     return p1
    # end

    # function nseval()
    #     """
    #     nse barplot with values > 0
    #     """
    #     ds = kge_df3()
    #     ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.name)
    #     dfi = filter(row -> row.NSE .> 0, ds)
    #     ann = map(x->string.(round(x;sigdigits=2)),dfi.NSE)
    #     Plots.bar(dfi.name, dfi.NSE, 
    #     xlabel = "Name", ylabel = "NSE", legend = false,
    #     title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600),
    #     fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
    #     annotations = (dfi.name,dfi.NSE, ann, :top),
    #     left_margin = 10mm,
    #     bottom_margin = 15mm, 
    #     bar_width = 0.6)        
    # end

    # function nsevalraw()
    #     """
    #     nse barplot with all values
    #     """
    #     ds = kge_df3()
    #     ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.name)
    #     dfi = ds
    #     ann = map(x->string.(round(x;sigdigits=2)),dfi.NSE)
    #     Plots.bar(dfi.name, dfi.NSE, xlabel = "Name", ylabel = "NSE", legend = false, 
    #         title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
    #         fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
    #         annotations = (dfi.name,dfi.NSE, ann, :top),
    #         bar_width = 0.6)        
    # end

    function waba()
        wpth="C:/Users/Public/Documents/Python_Scripts/julia/water-balance.jl"
        include(wpth)
        #img=load("waba-jl.png")
        yd=waread("waba-input.wa")|>yrsum
        @warn "try baryr(yd) !"
        baryr(yd)|>display
    end

    # function waba2()
    #     """
    #     calculates water balance and displays
    #     """
    #     begin
    #         af = filter(x -> occursin(r"^so_", x), readdir(pwd()))
    #         if length(af) <= 2 || any(ismissing.(af))
    #             error("match failed \n ... abort.\nno special output files present!")
    #             exit(86)
    #         end
        
    #         println("calculating yearly water balance of WaSiM special output data..\n
    #         bwvr = rain + snow + uprs - perc - qb - qd - qi - etr_ - ei_ - etrs_\n")
        
    #         re = Regex("preci|snow_storage_tota|Capi|Perc|baseflo|directflo|interflo|real_evap|real_tran|interception_evaporatio|snow_evaporatio","i")
    #         my = filter(x -> occursin(re, x), af)
    #         printstyled("loading...\n $my\n",color=:green)
        
    #         if (length(my) .!= 11)==true
    #             lng=length(my)
    #             printstyled("found only $lng files...\n $my\n",color=:yellow)
    #             error("\nfiles are missing!\ncheck so files...")
    #         end
        
    #         rain = filter(x -> occursin("precip", x), af)|>only|>so_read
    #         snow = filter(x -> occursin("snow_storage_total", x), af)|>only|>so_read
    #         uprs = filter(x -> occursin("Capil", x), af)|>only|>so_read
    #         perc = filter(x -> occursin("Perco", x), af)|>only|>so_read
    #         qb = filter(x -> occursin("baseflow", x), af)|>only|>so_read
    #         qd = filter(x -> occursin("directflow", x), af)|>only|>so_read
    #         qifl = filter(x -> occursin("interflow", x), af)|>only|>so_read
    #         etr = filter(x -> occursin("real_evapo", x), af)|>only|>so_read
    #         etrans = filter(x -> occursin("real_trans", x), af)|>only|>so_read
    #         ei = filter(x -> occursin("interception_evaporation", x), af)|>only|>so_read
    #         etrs = filter(x -> occursin("snow_evaporation",x), af)|>only|>so_read
        
    #         l = [rain, snow, uprs, perc, qb, qd, qifl, etr, etrans, ei, etrs]
    #         #typeof(l)
    #         nm = [names(l[i])[1] for i in 1:size(l, 1)]
    #         #same:
    #         #map(x->names(x)[1],l)
    #         println("loaded dataframes:\n$nm")
        
    #         if (length(l) .!= 11)==true
    #             error("files are missing!\ncheck so files...")
    #         end
        
    #         d = mall(l)
    #         xd=copy(d)
    #         writewa("waba-input.wa",xd)
        
    #         #dyr = yrsum(d)
        
    #         pos = d[!,Cols(r"date|^(prec)|^(snow_stora)|^(Capi)")]
    #         pos = yrsum(pos)
    #         # calculate the sum of each row
    #         psum = DataFrame(
    #             possums = [sum(eachrow(pos)[i]) for i in 1:size(pos, 1)],
    #             year=pos[!,:year]
    #         )
        
        
    #         neg = d[!,Not(Cols(r"^(prec)|^(snow_stora)|^(Capi)"))]
    #         neg = sum.(yrsum(neg))
        
        
    #         nsum = DataFrame(
    #             negsums = [sum(eachrow(neg)[i]) for i in 1:size(neg, 1)],
    #             year=neg[!,:year]
    #         )
        
        
    #         bw = innerjoin(psum, nsum, on=:year)
    #         bw.bw = bw[!,:possums] .- bw[!,:negsums]
    #     end
    #     ti="water-balance of "*basename(pwd())
    #     #theme(:vibrant)
    #     #theme(:wong)
    #     #theme(:dao)       #latex fonts.
    #     theme(:mute)
    #     #theme(:sand)
    #     ann = map(x->string.(round(x;sigdigits=3))*" [mm]",bw.bw)
    #     fact=.60
    #     plotsize = (1600*fact,800*fact)
    #     p1 = @df bw Plots.plot(
    #         :year,:bw,
    #         #annotations =(bw.year, bw.bw, ann, :top),
    #         #annotations = (bw.year,bw.bw,(ann,10,:center,:top,:black)),
    #         #annotations = (bw.year, bw.bw, ann, 8, :left, :top, :black),
    #         legend = false, 
    #         seriestype=:bar,
    #         xticks = bw.year,
    #         xtickfont = 12,
    #         xlabel = "",
    #         ylabel = "[mm]",
    #         ytickfont = 12,
    #         title = ti,
    #         fillcolor = ifelse.(bw.bw .> 0, "cornflowerblue", "coral2"),
    #         size=plotsize,
    #         #xrotation = 60);
    #         left_margin = 10mm,
    #         bottom_margin = 10mm, 
    #         #bottom_margin = 10px, 
    #         xrotation = 45)
            
    #     #Plots.annotate!(bw.year,bw.bw,(ann,10,:left,:top,:black))
    #                 #text("hey", 14, :left, :top, :green)
    #     for i in 1:length(bw.year)
    #         Plots.annotate!(bw.year[i],bw.bw[i],(ann[i],10,:center,:top,:black))
    #         #println(ann[i]*" added")
    #     end
    #     #display(p1)
    #     #return(bw)
    #     return p1
    # end

    ##very nice meta programming... 
    macro vv(s) vgjl(s);end
    #@vv "unter"
    macro vpy(s) vgpy(s);end
    #@vpy "climate"
    macro vr(s) vgr(s);end
    #@vr "climate"
    macro vct(s) vgctl(s);end
    #@vct "das ist"
    macro rg(s) rglob(s);end
    macro glb(s) glob(s);end
    macro gl(s) glob(s)|>first;end
    #fastplot
    macro fp(s) dfp(Regex(s));end
    macro flog(s) dfl(Regex(s));end
    macro ncrm() ncrem="C:/Users/Public/Documents/Python_Scripts/julia/ncremover.jl";include(ncrem);end
    macro rasrm() remer="C:/Users/Public/Documents/Python_Scripts/julia/raster_remover.jl";include(remer);end
    macro nco(s) nconly(s);end
    
    macro wajs() pt="C:\\Users\\Public\\Documents\\Python_Scripts\\julia\\wajs.jl";include(pt);end

    function readbetween(io::IO, start::String, stop::String)
        output = Vector{String}()
        while !eof(io)
            line = readline(io)
            if contains(line, start)
                push!(output, line)
                while !eof(io)
                    line = readline(io)
                    if contains(line, stop)
                        skip(io, length(line))
                        break
                    end
                    push!(output, line)
                end
                break
            end
        end
        return output
    end


    function rds(filename::String)
        """
        skips first line after [soil_table] i.e. no of soil types
        now returns a DataFrame
        """
        data = open(filename) do io
            a = readbetween(io, "soil_table", "substance_transport")
            return(a)
        end
        data = broadcast(x -> replace(x,    
        r"^#.*" => "",
        r"^[[].*" => "",
        r"method" => "",
        r"MultipleHorizons" => "",
        r"}" => "",
        r" = " => "=" ), data)


        filter!(s -> !isempty(s), data)
        lines = data[2:end]
        result = []
        # iterate over each line
        for line in lines
            # skip empty lines
            if isempty(line)
                continue
            end
            
            # split the line into number and fields
            number, fields = split(line, " {", limit=2)
            
            # convert the number to an integer
            number = parse(Int, number)
            
            # split the fields into individual fields
            fields = split(fields, ';')
            
            # initialize a dictionary to store the data for this line
            data = Dict{String, Any}()
            
            # store the number in the dictionary
            data["number"] = number
            
            # iterate over each field
            
            for field in fields
                # check if the field contains the " = " substring
                if occursin("=", field)
                    # split the field into key and value
                    key, value = split(field, "=")
                    
                    # check if the key is "Name"
                    if key == "Name"
                        # keep the value as a string
                    else
                    try # convert the value to an array of floats
                        
                      
                        # check if the value is a number
                        if occursin(r"^-?\d+(\.\d+)?$", value)
                            # convert the value to a float
                            value = parse(Float64, value)
                        elseif occursin(r"^-?\d+(\.\d+)?\d+?$", value)
                            # convert the value to a float
                            value = parse(Float64, value)
                        elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                            # convert the value to a float (scientific notation)
                            value = parse(Float64, value)
                        elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                            # convert the value to a float (scientific notation)
                            value = parse(Float64, value)
                        elseif occursin(r"^\d+ \d+", value)
                            # convert the value to an array of integers
                            value = parse.(Int, split(value))
                            #value = parse.(Float64, split(value))
                        elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)

                            value = parse.(Float64, split(value))
                                                    
                        end
                    catch
                        @warn "could not parse $value"
                        continue
                    end

                    end

                    
                    # store the key-value pair in the dictionary
                    data[key] = value
                end
            end

            
            # append the dictionary to the result array
            push!(result, data)
        end
    
        return DataFrame(result)
    end
        
    function read_soildata(filename::String)
        """
        skips first line after [soil_table] i.e. no of soil types
        now returns a DataFrame
        """
        data = open(filename) do io
            a = readbetween(io, "soil_table", "substance_transport")
            return(a)
        end
        data = broadcast(x -> replace(x,    
        r"^#.*" => "",
        r"^[[].*" => "",
        r"method" => "",
        r"MultipleHorizons" => "",
        r"}" => "",
        r" = " => "=" ), data)


        filter!(s -> !isempty(s), data)
        lines = data[2:end]
        result = []
        # iterate over each line
        for line in lines
            # skip empty lines
            if isempty(line)
                continue
            end
            
            # split the line into number and fields
            number, fields = split(line, " {", limit=2)
            
            # convert the number to an integer
            number = parse(Int, number)
            
            # split the fields into individual fields
            fields = split(fields, ';')
            
            # initialize a dictionary to store the data for this line
            data = Dict{String, Any}()
            
            # store the number in the dictionary
            data["number"] = number
            
            # iterate over each field
            try
            for field in fields
                # check if the field contains the " = " substring
                if occursin("=", field)
                    # split the field into key and value
                    key, value = split(field, "=")
                    
                    # check if the key is "Name"
                    if key == "Name"
                        # keep the value as a string
                    else
                      
                        # check if the value is a number
                        if occursin(r"^-?\d+(\.\d+)?$", value)
                            # convert the value to a float
                            value = parse(Float64, value)
                        elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                            # convert the value to a float (scientific notation)
                            value = parse(Float64, value)
                        elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                            # convert the value to a float (scientific notation)
                            value = parse(Float64, value)
                        elseif occursin(r"^\d+ \d+", value)
                            # convert the value to an array of integers
                            value = parse.(Int, split(value))
                            #value = parse.(Float64, split(value))
                        elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)

                            # convert the value to an array of floats
                            value = parse.(Float64, split(value))
                                                    
                        end

                    end
                    
                    # store the key-value pair in the dictionary
                    data[key] = value
                end
            end
        catch
            @warn "could not parse $value"
            continue
        end
            
            # append the dictionary to the result array
            push!(result, data)
        end
    
        return DataFrame(result)
    end

    function fsize(m::String)
        """
        names and size in MB via readdir
        """
        files = filter(x -> occursin(Regex(m,"i"), x), readdir())
        fzs = [(file, filesize(file) / 2^20) for file in files]
        tot = round(sum(map(x->x[2],fzs));digits=3)
        #println("Total Size: $tot MB")
        printstyled("Total Size: $tot MB\n",color=:green)
        return(fzs)
    end

    function cntcols(x::String)
        # Get a list of all files in the directory
        #x = "."
        files = filter(file -> (occursin(Regex(x, "i"), file) & 
                            (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file))
                            ), readdir())

        files = filter(inF->isfile(inF),files)
                            #if isfile(inF)
        file_columns = []
        
        for file in files
            # Open the file for reading
            open(file, "r") do io
                # Read the first line of the file
                line = readline(io)
                # Split the line on tabs
                columns = split(line, '\t')
                # Count the number of columns
                num_columns = length(columns)
                # Store the file name and number of columns as a tuple
                push!(file_columns, (file, num_columns))
            end
        end
        
        # Sort the file_columns array based on the number of columns in descending order
        sorted_files = sort(file_columns, by = x -> x[2], rev = true)
        
        for (file, num_columns) in sorted_files
            #println("File: ", file, " | Columns: ", num_columns)
            #printstyled(rpad("$file:",50),color=:light_magenta)
            printstyled(
                rpad("File: $file",45),
            lpad(" | Columns: $num_columns\n",10),color=:green,bold=true)
        end
    end

    # cntcols("qoutjl")
    # cntcols("v0")
    function rmeq()
        """
        removes empty TS; 
        use with caution!
        """
        #x = pwd()

        # files = filter(file -> (occursin(Regex(x, "i"), file) & 
        # (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file))
        # ), readdir())
        
        files = filter(file -> 
        !occursin(r"tex|pl|sh|csv|html|xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file)
        , readdir())
        
        ms = ["-9999", "lin", "log", "--"]
        for inF in files
            if isfile(inF)
                df = CSV.File(inF; delim="\t", header=1,
                silencewarnings=true, 
                    normalizenames=false, 
                    missingstring=ms, 
                    types=Float64) |> DataFrame
                dropmissing!(df,ncol(df))
                if nrow(df)==0
                    println(basename(inF)," removed!")
                    rm(inF)
                end
            end
        end
    end

    function getnames(dfs::Vector{Any})
        try 
            nms=[]
            for df in dfs
                x=collect(DataFrames.metadata(df))[1][2]|>basename
                println(x)
                push!(nms,x)
            end
            return nms
        catch
            @error "no metadata in $df !"
            return
        end
    end

    function kgedf()
        path = pwd()
        files = glob(r"_output.txt|_outputjl") #non-recurse
        results = DataFrame(File = String[], KGE = Float64[])  # Create an empty DataFrame to store the results
        
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file, '\t', String)
            match = Grep.grep(r"KGE", output)
            
            if !isempty(match)
                fn = first(split(file, "_qout"))
                parsed_lines = Float64[]
                
                for line in match
                    line_parts = split(line)
                    kge_value = parse(Float64, line_parts[end])
                    push!(parsed_lines, kge_value)
                end
                
                sort!(parsed_lines, rev = true)  # Sort the parsed KGE values in descending order
                
                for kge_value in parsed_lines
                    push!(results, (File = fn, KGE = kge_value))  # Add the result to the DataFrame
                end
            end
        end
        
        sort!(results, :KGE, rev = true)  # Sort the DataFrame by the 'KGE' column in descending order
        return results
    end

    # function nsx(dfi::DataFrame)
    #     """
    #     nse barplot with all values
    #     """
    #     dfi.name=map(x->replace(x,r"-qoutjl.*" => ""),dfi.name)
    #     ann = map(x->string.(round(x;sigdigits=3)),dfi.NSE)
    #     Plots.bar(dfi.name, dfi.NSE, xlabel = "Name", ylabel = "NSE", legend = false, 
    #         title = "x", xrotation = -25, fmt = :png, size = (800, 600), 
    #         fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
    #         annotations = (dfi.name,dfi.NSE, ann, :top),
    #         xtickfont = font(6),    
    #         bar_width = 0.77)        
    # end

    function ctlg(dir_path::String, match::String;file_ending="ctl")
        """
        dir_path::String, match::String;file_ending="ctl"
        """
        for file in readdir(dir_path)
            if occursin(file_ending, file)
                fx = joinpath(dir_path, file)
                prevline = ""
                for (i, line) in enumerate(eachline(fx))
                    if findfirst(match, line) !== nothing
                        printstyled("File: $file\n",color=:red)
                        println("$(prevline)")
                        printstyled("$line\n",color=:green)
                        nextline = readline(fx)
                        println("$(nextline)")
                    end
                    prevline = line
                end
            end
        end
    end

    macro sdjl() fx="C:/Users/Public/Documents/Python_Scripts/julia/sd2.jl";include(fx);end
    #@sdjl
    macro sf(s) glob(s);end
    macro listdir() ls();end

    function rmqout()
        map(x->rm(x),glob("qoutjl"))
    end

    function renamer(df::DataFrame)
        """ 
        renamer - remove beginning char _   
        and replace it with C
        """
        for x in names(df)
            if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    function read_until_flag(file::IOStream, flag::String)
        line = readuntil(file, flag)
        return line[1:end-length(flag)]
    end


    function findindf(df::DataFrame, x)
        """
        like Grep.grep("x",df)
        Regex(x, "i"), 
        """
        filter(row -> any(occursin(x,
            string(value)) for value in row), 
                eachrow(df))
    end

    function findindf(df::DataFrame, x::Regex)
        """
        like Grep.grep("x",df)
        """
        filter(row -> any(occursin(x, 
            string(value)) for value in row), 
                eachrow(df))
    end

    # findindf(df,"Et")
    # findindf(df,"15")
    # findindf(df,"-")
    # Grep.grep(r"Et",df.name)
    # Grep.grep(r"Et",df)
    # Grep.grep("Et",df)


    function wslpath()
        """
        import InteractiveUtils.clipboard
        wslpath()|>clipboard
        cb(wslpath())
        """
        # Run the `wslpath` command to convert the current directory to a WSL path
        #wsl_cmd = `wsl wslpath -a $(pwd)`
        wsl_cmd = `wsl wslpath -a .`
        wsl_path = readchomp(pipeline(wsl_cmd))
        # Return the WSL path
        return wsl_path
    end

    function wslp(winpt::AbstractString)
        """
        """
        #printstyled("needRAW like raw",color=:red)
        #return(raw$winpt)
        #quote($winpt);end
        #winpt=raw"D:\Wasim\regio\out\lowres\c5\loc5"
        winpt = replace(winpt, "\\"=> "/")
        #winpt = join(winpt,"'")
        #print("\"",winpt,"\"")
        #join(wsl_cmd,winpt)
        wsl_cmd = `wsl wslpath -m $winpt`
        run(wsl_cmd,winpt)
        #wsl_path = readchomp(pipeline(wsl_cmd))
        # Return the WSL path
        return wsl_path
    end

    function xread(filename::String)
        """
        takes first line as header and drops r"[A-z]"
        """
        ms=["-9999","lin","log","--"]
        df = CSV.read(filename,DataFrame, 
            delim="\t",
            missingstring=ms,
            header=1, limit=10^5)
        if typeof(df[1, 1]) != Int64
            df = df[map(x->!occursin(r"[A-z]", x),df[:, 1]), :]
            # map to int for dates
            for i in 1:3
                df[!,i]=map(x ->parse(Int,x),df[!,i])
            end
            #and parse dates...
            df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
            df=df[:,Not(1:4)]
            metadata!(df, "filename", filename, style=:note); 
        else
            df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
            df=df[:,Not(1:4)]
            metadata!(df, "filename", filename, style=:note);  
        end
    end  

    function getm(s::Any)
        """
        grabs methods
        asin|>getm  
        ?asin
        @code_llvm readf|>getm|>first 
        """
        methods(s);
    end

    function sdf()
    
        function calculate_folder_size(directory)
            size = 0
            count = 0
            for (root, dirs, files) in walkdir(directory)
                for file in files
                    size += stat(joinpath(root, file)).size
                    count += 1
                end
            end
            return size, count
        end
    
        function print_folder_size(directory, size, count)
            size_gb = round(size / 1024^3, digits=3)
            printstyled(rpad("$directory: $size_gb GB", 40), color=:green)
            printstyled(lpad("($count files)", 20), "\n", color=:green)
        end
    
        printstyled("folder sizes on Julia...\n", color=:red)
    
        cwd = pwd()
        dirs = readdir(cwd)
    
        rs, rcnt = calculate_folder_size(cwd)
    
        print_folder_size(cwd,rs,rcnt)
    
        n = repeat(" - -", 10)
        printstyled(n * " subfolders of " * basename(cwd) * n, "\n", color=:yellow)
    
        for dir in dirs
            if isdir(dir)
                size, count = calculate_folder_size(joinpath(cwd, dir))
                print_folder_size(dir, size, count)
            end
        end
    end

    function pfix(m::String)
        """
        names and size in MB via readdir
        sorts and stores in df -> biggest files
        """
        files = filter(x -> occursin(Regex(m,"i"), x), readdir())
        fzs = [(file, filesize(file) / 2^20) for file in files]
        tot = round(sum(map(x->x[2],fzs));digits=3)
        printstyled("Total Size: $tot MB\n",color=:green)
        df = DataFrame(fzs)
        sort!(df,2,rev=true)
        return(df)
    end
    
    function pfix()
        """
        names and size in MB via readdir
        sorts and stores in df -> biggest files
        """
        files = readdir()
        fzs = [(file, filesize(file) / 2^20) for file in files]
        tot = round(sum(map(x->x[2],fzs));digits=3)
        printstyled("Total Size: $tot MB\n",color=:green)
        df = DataFrame(fzs)
        sort!(df,2,rev=true)
        return(df)
    end

    function mvwasim2() 
        
        println("\nmoves all wq, xml and log files to from
        c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05
        to current pwd")
        ta=pwd()
        pt="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05";
        println("target dir is $ta");
        #@vv "af"
    
        af = filter(x -> occursin(r"wq", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i)))
            println(basename(i)," --> ", ta)
        end
        #@rg "wq"
    
        af = filter(x -> occursin(r"xml", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i)))
            println(basename(i)," --> ", ta)
        end
        af = filter(x -> occursin("modell", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i)))
            println(basename(i)," --> ", ta)
        end
    
    end
    
    # function qplot(df::DataFrame)
    #     """
    #     takes first two cols of df and plots r2 QQ
    #     """
    #     if any(map(x->contains(x,"date"),names(df)))
    #         df = df[!,Not(:date)]
    #     end
    #     if ncol(df)>2
    #         df = df[!,1:2]
    #     end
    #     r2 = round(cor(df[!,1], df[!,2])^2, digits=3)
    #     p = df |> x -> qqplot(x[!,1],x[!,2], 
    #         #title = "R² = "*string(ti),
    #         qqline = :fit)
    #         #color = :grays) # erstellt ein QQ-Diagramm <- black
    #     xlabel!(p,names(df)[1])
    #     ylabel!(p,names(df)[2])
    #     annotate!(p,:bottomright, text("R² = "*string(r2), :black))
    #     #xr = maximum(df[!,1])
    #     #annotate!(xr-0.5, 2.0, text("R² = $ti", 12))
    # end

    # function dfbar(df::DataFrame)
    #     ti = try
    #             DataFrames.metadata(df)|>only|>last|>basename
    #         catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #     end
    #     if (any(x->occursin("year",x),names(df)))
    #         s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
    #         #@df df Plots.bar(:year,cols(s),legend = :topright, title=ti)
    #         @df df groupedbar(df.year,cols(s), legend = :outertopright, title=ti)
    #     else    
    #     s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    #     #@df df Plots.bar(:date,cols(s),legend = :topright, title=ti)
    #     @df df groupedbar(df.date,cols(s), legend = :outertopright, title=ti)
    #     end
    # end

    function qba()
        """
        return 1st DFs of qgko
        """
        x = first(glob("qgko"))
        println("loading $x ...")
            try
                df = DataFrame(CSV.File(x, header=false, 
                                    delim="\t",
                                    ignorerepeated=true,
                                    silencewarnings=true,
                                    limit=10^4,
                                    types = String))
                
                pattern = r"^[LIN. R]|^[LOG. R]|^CO"
                mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
                dx = df[mask, :]
                dx = permutedims(dx) |>dropmissing
                basins = []
                for i in copy(df[1,5:end])
                    push!(basins,parse.(Int,i))
                    #push!(basins,string.("B_"*i))
                end
                insert!(basins, 1, "score")
                insert!(basins, 2, "timestep")
                dx[!, "basin"] = basins
                cn = (dx[1,:])
                rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
                dout = dx[3:end,:]
                #reorder columns
                dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
                DataFrames.metadata!(dout, "filename", x, style=:note);
                return(dout)
            catch
                @error("error! ")
        end
    end

    function wread(x::String;skip=3)
        """
        Read wasim ts with DelimitedFiles.readdlm, skipto line 3 
        no header column
        """
        df = DelimitedFiles.readdlm(x, '\t', Float64, '\n';
            header=false,skipstart=skip)
        df = DataFrame(df,:auto)
        for i in 5:size(df,2)
            df[!,i]=replace(df[!,i],-9999.0 => missing)
        end 
        for i in 5:size(df,2)
            replace!(df[!,i],-9999.0 => missing)
        end
        for i in 1:3
            df[!,i]=map(x ->Int(x),df[!,i])
        end
        #and parse dates...
        df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
        df=df[:,Not(1:4)]
        metadata!(df, "filename", x, style=:note);
    end
    
    dfr = waread

    function grec(dir_path::String, file_ending::String, match::String)
        """
        grep recursive and looks for matches in files
        """
        for (looproot, dirs, filenames) in walkdir(dir_path)
            for file in filenames
                fx = joinpath(looproot, file)
                if isfile(fx) && occursin(Regex(file_ending*"\$","i"), fx)
                    #printstyled("check $file\n",color=:yellow)
                    #fx = joinpath(looproot, file)
                    prevline = ""
                    for (i, line) in enumerate(eachline(fx))
                        if findfirst(match, line) !== nothing
                            printstyled("File: $fx\n",color=:red)
                            println("$(prevline)")
                            printstyled("$line\n",color=:green)
                            nextline = readline(fx)
                            println("$(nextline)")
                        end
                        prevline = line
                    end
                end
            end
        end
    end

    # function vio(df::DataFrame)
    #     """
    #     with mon abbr. see ovio for numbered months
    #     """   
    #         ti = try
    #             DataFrames.metadata(df)|>only|>last|>basename
    #         catch
    #         @warn "No basename in metadata!"
    #         ti = raw""
    #         end    
    #         str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    #         month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    #         if (any(x->occursin("year",x),names(df)))
    #             df = df[!,Not("year")]
    #             s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #             @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
    #             xticks!(0.5:11.5 , month_abbr)
    #             title!(ti)
    #         else    
    #         s = Symbol.(filter(x->!occursin("date",x),names(df)))
    #         @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
    #         xticks!(0.5:11.5 , month_abbr)
    #         title!(ti)
    #     end
    # end
    
    # function mbx(df::DataFrame)
    #     """
    #     annotated boxplot with mean
    #     """
    #     str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    #     ln = Symbol.(filter(x->!occursin(r"date|year|month"i,x),names(df)))
    #     month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    #     p = @df df StatsPlots.boxplot(str,cols(ln),
    #         fillalpha=0.75, 
    #         linewidth=0.25,
    #         notch = true,
    #         whisker_width = :match,
    #         legend=false)
    #     xticks!(0.5:11.5 , month_abbr)
    #     df.Month = month.(df.date)
    #     means = DataFrames.combine(groupby(df,:Month ), ln[1] => mean)
    #     #mean(means[!,2])
    #     for i in eachrow(means)
    #         m = i[2]
    #         annotate!(i.Month - 0.5, m, #+ 1 
    #         text(round(m; digits=2), 6, :center, :top))
    #     end
    #     return p
    # end

    # function dprbig(x::DataFrame)
    #     """
    #     correlation plots on dataframe
    #     size=(1200,800)
    #     has to be sim,obs
    #     """
    #     df = copy(x) #<:DataFrame

    #     if any(map(x->occursin("year",x),names(df)))
    #         df = df[!,Not(:year)]
    #     end
        
    #     if any(map(x->occursin("month",x),names(df)))
    #         df = df[!,Not(:month)]
    #     end
        
    #     if propertynames(df)[end]!=:date
    #         df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
    #     end
        
    #     dropmissing!(df)

    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, 
    #     xlabel="Date", ylabel="[mm/day]",
    #     legend = :topleft,
    #     size=(1200,800)
    #     )
    #     r2 = round(cor(df[!,1], df[!,2])^2, digits=2)
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    # function dprbig(x::Regex)
    #     """
    #     correlation plots on dataframe
    #     size=(1200,800)
    #     """
    #     df = globdf(x)|>first|>waread
    #     df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
    #     v = (names(df[!,1:2]))
    #     a = reshape(v, 1, 2)
    #     Plots.plot(df.date,[df[!,1], df[!,2]], 
    #     label=a, 
    #     xlabel="Date", ylabel="[mm/day]",
    #     legend = :topleft,
    #     size=(1200,800)
    #     )
    #     r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    #     kge = round(kge2(df[!,1], df[!,2]), digits=2)
    #     nse_value = round(nse(df[!,1], df[!,2]), digits=2)
    #     annotate!(
    #         last(df.date), 0.95*maximum(df[!,1]),
    #     text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
    #     )
    # end

    # function qplot(x::Vector{Float64},y::Vector{Float64})
    #     # if any(map(x->contains(x,"date"),names(df)))
    #     #     df = df[!,Not(:date)]
    #     # end
    #     # if ncol(df)>2
    #     #     df = df[!,1:2]
    #     # end
    #     r2 = round(cor(x, y)^2, digits=3)
    #     p = qqplot(x,y,
    #         qqline = :fit)
    #     # xlabel!(p,names(df)[1])
    #     # ylabel!(p,names(df)[2])
    #     annotate!(p,:bottomright, text("R² = "*string(r2), :black))
    # end

    function rowsums(df::DataFrame)
        return DataFrames.combine(df, names(df) .=> sum)
    end
    
    function rowmeans(df::DataFrame)
        return DataFrames.combine(df, names(df) .=> mean)
    end
    
    function addname(indf::DataFrame,nmdf::DataFrame)
        """
        addname(indf::DataFrame,nmdf::DataFrame)
        indf: input dataframe
        nmdf: name dataframe

        indf = dfonly(r"gwst")|>first|>waread
    
        ofl="route.txt"
        begin
            df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
            rename!(df,1=>"sim",2=>"obs",3=>"name")
            df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
            df.name=map(x->replace(x,r"_>.*" => ""),df.name)
            sort!(df, :sim)
        end
        
        addname(indf,nmdf::df)
        select(indf,Not(Cols(r"^[0-9]")))|>dfp
    
        """
        rename!(indf,map(x->replace(x,r"^C" => ""),names(indf)))
        for row in eachrow(nmdf)
            old_name = string(row.sim)
            new_name = string(row.name)
            if old_name in names(indf)
                println("renaming $old_name")
                rename!(indf, (old_name) => (new_name*"_"*old_name))
            end
        end
        return indf
    end

    function rmdub(;directory=pwd())
        """
        recursively removes duplicates 
        uses SHA
        """
        
        # Create a dictionary to store file hashes as keys and file paths as values
        #hash_dict = Dict{String, String}()
        hash_dict = Dict{Vector{UInt8}, String}()
    
        # Get a list of all files in the directory
        #files = readdir(directory)
        for (root, dirs, files) in walkdir(directory)
            for file in files
                filepath = joinpath(directory, file)
                if isfile(filepath)
                    # Calculate the SHA-256 hash of the file's contents
                    #hash = string(sha256(open(filepath, "r")))
                    io = open(filepath, "r")
                    #filehash = string(sha256(io))
                    #filehash = string(hash(io))
                    filehash = sha256(io)
                    #println(filehash)
                    close(io)
    
                    # If the hash is not already in the dictionary, add it
                    if !haskey(hash_dict, filehash)
                        hash_dict[filehash] = filepath
                    else
                        # If a file with the same hash is found, delete it
                        println("Deleting duplicate file: $filepath")
                        rm(filepath)
                    end
                end
            end
        end
    end

    function irfan(fl::String)
        opener="C:/Program Files/IrfanView/i_view64.exe"
        run(`$opener $fl`)
    end

    function extract_duration_from_xml(xml_file::AbstractString)
        finished_timestamp = ""
        start_timestamp = ""
    
        # Read the content of the XML file
        content = readlines(xml_file)
    
        # Find the lines containing "WaSiM finished" and "WaSiM start"
        finished_lines = filter(x -> occursin(r"WaSiM finished", x), content)
        start_lines = filter(x -> occursin(r"WaSiM start", x), content)
    
        # Extract the timestamps from the lines if found
        if !isempty(finished_lines)
            finished_line = first(finished_lines)
            finished_timestamp = match(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}", finished_line).match
        end
    
        if !isempty(start_lines)
            start_line = first(start_lines)
            start_timestamp = match(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}", start_line).match
        end
    
        # Calculate the duration in hours and minutes
        if !isempty(finished_timestamp) && !isempty(start_timestamp)
            finished_time = DateTime(finished_timestamp)
            start_time = DateTime(start_timestamp)
            duration_milliseconds = Dates.value(finished_time - start_time)
    
            # Convert duration to hours and minutes
            duration_hours = div(duration_milliseconds, 3_600_000)  # 1 hour = 3,600,000 milliseconds
            duration_minutes = div(rem(duration_milliseconds, 3_600_000), 60_000)  # 1 minute = 60,000 milliseconds
    
            return (duration_hours, duration_minutes)
        else
            return nothing
        end
    end

    function process_folders_and_subfolders(root_dir::AbstractString, log_file::AbstractString = "output.log")
        open(log_file, "w") do io
            for (looproot, dirs, filenames) in walkdir(root_dir)
                for filename in filenames
                    if occursin(r"xml", filename)
                        entry = joinpath(looproot, filename)
                        result = extract_duration_from_xml(entry)
    
                        if result !== nothing
                            folder_name = dirname(entry)
                            println(io, "XML File: $entry")
                            hours, minutes = result
                            println(io, "Duration: $hours hours and $minutes minutes")
                            println(io, "="^50)
                        end
                    end
                end
            end
        end
    end

    function read_log_file(log_file::AbstractString)
        timestamps = []
        durations = []
    
        open(log_file) do io
            while !eof(io)
                line = readline(io)
                if startswith(line, "XML File: ")
                    filename = split(line, "XML File: ")[2]
                    push!(timestamps, filename)
                elseif startswith(line, "Duration: ")
                    duration_str = split(line, "Duration: ")[2]
                    duration_parts = split(duration_str, " hours and ")
                    hours = parse(Int, duration_parts[1])
                    minutes = parse(Int, split(duration_parts[2], " minutes")[1])
                    push!(durations, hours * 60 + minutes)
                end
            end
        end
    
        return timestamps, durations
    end
    
    function plot_duration_graph(log_file::AbstractString)
        timestamps, durations = read_log_file(log_file)
    
    
        # Sort data by timestamps
        sorted_indices = sortperm(timestamps)
        sorted_timestamps = timestamps[sorted_indices]
        sorted_durations = durations[sorted_indices]

        sorted_timestamps = map(x->
            replace(basename(x), "xml" => splitdir(dirname(x))[2]), 
                sorted_timestamps)
    
        # Create bar plot
        p = plot(
            sorted_timestamps,
            sorted_durations,
            seriestype = :bar,
            #xlabel = "XML Files",
            xlabel = "",
            ylabel = "Duration (minutes)",
            title = "Duration of WaSiM runs",
            legend = false,
            xrotation = -35,  # Rotate x-axis labels for better readability
            grid = true,
            right_margin = 10mm,
            left_margin = 5mm,
            top_margin = 5mm,
            bottom_margin = 15mm, 
            bar_width = 0.6);
        
    
        # Add annotations for duration time in minutes
        for (x, y) in zip(sorted_timestamps, sorted_durations)
            annotate!(x, y, Plots.text("$y min",7,:black,:bottom,
            rotation=0            
            ))
        end
    
        return p
    end

    function plot_duration_graph(log_file::Tuple{Vector{Any}, Vector{Any}})
        timestamps, durations = log_file
        # Sort data by timestamps
        sorted_indices = sortperm(timestamps)
        sorted_timestamps = timestamps[sorted_indices]
        sorted_durations = durations[sorted_indices]

        sorted_timestamps = map(x->
            replace(basename(x), "xml" => splitdir(dirname(x))[2]), 
                sorted_timestamps)
    
        # Create bar plot
        p = plot(
            sorted_timestamps,
            sorted_durations,
            seriestype = :bar,
            #xlabel = "XML Files",
            xlabel = "",
            ylabel = "Duration (minutes)",
            title = "Duration of WaSiM runs",
            legend = false,
            xrotation = -35,  # Rotate x-axis labels for better readability
            grid = true,
            right_margin = 10mm,
            left_margin = 5mm,
            top_margin = 5mm,
            bottom_margin = 15mm, 
            bar_width = 0.6);
        
    
        # Add annotations for duration time in minutes
        for (x, y) in zip(sorted_timestamps, sorted_durations)
            annotate!(x, y, Plots.text("$y min",7,:black,:bottom,
            rotation=0            
            ))
        end
    
        return p
    end

    #function kgewrite(;output_file="kge-table_"*last(splitpath(pwd()))*".txt")
    function kgewrite(;output_file="kge-table.txt")
        path = pwd()
        files = glob(r"_output.txt|_outputjl") #non-recurse
        
        # Create an empty vector to store the matched lines
        matched_lines = Vector{String}()
        
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file, '\t', String)
            match = Grep.grep(r"KGE", output)
            if !isempty(match)
                fn = first(split(file, "_qout"))
                for line in sort(match, by = x -> parse(Float64, split(x)[end]); rev = true)
                    line = strip(line)  # remove leading and trailing whitespace
                    line = join(split(line), " ")  # remove inner whitespaces
                    push!(matched_lines, "$fn: $line")  # collect the matched line
                end
            end
        end
        
        #output_file = "matched_results.csv"
        writedlm(output_file, matched_lines, '\t')
        println("Matched results saved to $output_file")
    end

    function subset_dataframe_by_mask(df::DataFrame, msk::DataFrame)
        """
        map(typeof, eachcol(df)) #check types of cols
        msk = broadcast(x->typeof(x)==Vector{Float64},df)
        """
        # Get column names that satisfy the condition
        columns_to_keep = names(df)[collect(msk[1, :])]
        # Subset DataFrame using the mask
        subset_df = select(df, columns_to_keep)
        return subset_df
    end

    function getmoduleorder(file1::AbstractString)
        # Run the awk command and capture the STDOUT
        result = read(`awk '$0 ~ /^[[]/{c=1} c&&c--' $file1`, String)
    
        # Split the result into lines and strip \r and \t characters
        lines = replace(split(result, "\n"), r"[\r\t]" => "")
    
        # Extract the names inside square brackets using regular expressions
        pattern = r"\[(.*?)\]"
        names_inside_brackets = [match(pattern, line) === nothing ? "" : match(pattern, line).captures[1] for line in lines]
    
        # Create a DataFrame with the extracted names inside the square brackets
        df = DataFrame(Line = names_inside_brackets)
    
        return df
    end

    function ctg(match::String ; dir=".", file_ending=".ctl")
        for file in readdir(dir)
            if occursin(file_ending, file)
                fx = joinpath(dir, file)
                prevline = ""
                for (i, line) in enumerate(eachline(fx))
                    if findfirst(match, line) !== nothing
                        printstyled("File: $file\n",color=:red)
                        println("$(prevline)")
                        printstyled("$line\n",color=:green)
                        nextline = readline(fx)
                        println("$(nextline)")
                    end
                    prevline = line
                end
            end
        end
    end

    function all_values_equal(df::DataFrame)
        for col in eachcol(df)
            if any(col .!= col[1])
                return false
            end
        end
        return true
    end

    function mall(left::DataFrame, right::DataFrame)
        "reduces + merges by date"
        df = innerjoin(left, right, on = :date,makeunique=true)
        return(df)
    end

    function cpinto(src::Vector{String}, dst::AbstractString;force=false)
        """
        mkdir("route-bak")
        cpinto(glob("so_inf"), "route-bak")
        rglob("so_inf")
        force=true will first remove an existing dst.
        """
        map(x->cp(x,"$dst/$x";force=force),src)
    end

    function cntcolv(x::String)
        # Get a list of all files in the directory
        #x = "."
        files = filter(file -> (occursin(Regex(x, "i"), file) & 
                            (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file))
                            ), readdir())

        files = filter(inF->isfile(inF),files)
                            #if isfile(inF)
        file_columns = []
        
        for file in files
            # Open the file for reading
            open(file, "r") do io
                # Read the first line of the file
                line = readline(io)
                # Split the line on tabs
                columns = split(line, '\t')
                # Count the number of columns
                num_columns = length(columns)
                push!(file_columns, (file, num_columns))
            end
        end
        
        # Sort the file_columns array based on the number of columns in descending order
        sorted_files = sort(file_columns, by = x -> x[2], rev = true)
        
        for (file, num_columns) in sorted_files
            printstyled(
                rpad("File: $file",45),
            lpad(" | Columns: $num_columns\n",10),color=:green,bold=true)
        end
        return file_columns
    end


    function tovec(x::DataFrame, col::Any)
        df = select(x,col)
        println(names(df))
        return vec(Matrix(df))
    end


    function print_lines_between_patterns(filename::AbstractString, start_pattern::AbstractString, end_pattern::AbstractString)
        in_range = false
        
        for line in eachline(filename)
            if occursin(start_pattern, line)
                in_range = true
            end
            
            if in_range
                println(line)
            end
            
            if occursin(end_pattern, line)
                in_range = false
            end
        end
    end

    function read_soildata_2(filename::String)
        """
         soil types
         returns a DataFrame
         needs rework
        """
        #data = readlines(filename)
        # data = open(fl) do io
        #     a = readbetween(io,"{","}")
        #     return(a)
        # end
        data = open(filename) do io
            a = readbetween(io, "soil_table", "special_output")
            return(a)
        end
        
        #output = DelimitedFiles.readdlm(filename,';', String)
        data = Grep.grep(r"^th.*|^[0-9]",output)
        data = broadcast(x -> replace(x,    
                r"^#.*" => "",
                r"[[]].*" => "",
                r"[{].*" => "",
                r"method" => "",
                r"MultipleHorizons" => "",
                r"}" => ""), data)
        bks = broadcast(x -> split(x, " ",limit=3), data)
        nos = broadcast(x -> if length(x)==2 x[1] end, bks)
        nos = filter(x -> !isnothing(x),nos)
        nos = broadcast(x -> parse.(Int, x),nos)
        
        #tck = broadcast(x -> if (length(x)==3  & startswith(x[1],"thickness")) x[3] end, bks)
        tck = broadcast(x -> if (startswith(x[1],"thickness")) x[3] end, bks)
        filter!(x -> !isnothing(x),tck)
        tck = broadcast(x -> split(x),tck)

        Z=[]
        for m in (tck)
            i = [parse(Float64, x) for x in m]
            push!(Z,i)
        end
        
        #[cumsum(arr) for arr in Z]
        
        #zd = DataFrame(permutedims(Z),:auto)
        zd = DataFrame((Z),:auto)
        for x in 1:length(nos)
            newname = Symbol.(nos[x])
            rename!(zd,Dict(names(zd)[x]=>newname))
        end
        # col_sums = sum.(eachcol(zd))
        # hcat(col_sums*.1, nos)
       
    
        return zd
    end

    function read_soildata_raw(filename::String)
        """
         soil types
         returns a DataFrame
         needs rework
        """
        #data = readlines(filename)
        # data = open(fl) do io
        #     a = readbetween(io,"{","}")
        #     return(a)
        # end
        # data = open(filename) do io
        #     #a = readbetween(io, "soil_table", "special_output")
        #     a = readlines(io)
        #     return(a)
        # end
        
        data = DelimitedFiles.readdlm(filename,';', String)
        data = Grep.grep(r"^th.*|^[0-9]",data)
        data = broadcast(x -> replace(x,    
                r"^#.*" => "",
                r"[[]].*" => "",
                r"[{].*" => "",
                r"method" => "",
                r"MultipleHorizons" => "",
                r"}" => ""), data)
        bks = broadcast(x -> split(x, " ",limit=3), data)
        nos = broadcast(x -> if length(x)==2 x[1] end, bks)
        nos = filter(x -> !isnothing(x),nos)
        nos = broadcast(x -> parse.(Int, x),nos)
        
        #tck = broadcast(x -> if (length(x)==3  & startswith(x[1],"thickness")) x[3] end, bks)
        tck = broadcast(x -> if (startswith(x[1],"thickness")) x[3] end, bks)
        filter!(x -> !isnothing(x),tck)
        tck = broadcast(x -> split(x),tck)
        
        #foreach(x -> parse.(Int, x[1]),tck)
        
        #collect(tck)
        
        Z=[]
        for m in (tck)
            i = [parse(Float64, x) for x in m]
            push!(Z,i)
        end
        
        #[cumsum(arr) for arr in Z]
        
        #zd = DataFrame(permutedims(Z),:auto)
        zd = DataFrame((Z),:auto)
        for x in 1:length(nos)
            newname = Symbol.(nos[x])
            rename!(zd,Dict(names(zd)[x]=>newname))
        end
        # col_sums = sum.(eachcol(zd))
        # hcat(col_sums*.1, nos)
       
    
        return zd
    end   

    import InteractiveUtils: clipboard

    function mywd()
        wd = pwd()
        wd = replace(wd,"\\"=>"/")
        println("qouted $wd in clipboard!")
        #println("\"",wd,"\"")
        # wd|>clipboard
        return "\""*wd*"\"" |>clipboard
        #InteractiveUtils.clipboard()
    end

    function dfilter(df::DataFrame, col, val)
        DataFrames.subset(df,Symbol(col)=> ByRow(==(val)))
    end

    # function nqp(a::Regex,b::Regex;)
    #     #col="tot_average"
    #     #col = Symbol(col) 
    #     # a = select(waread(a),Cols(col,:date))
    #     # b = select(waread(b),Cols((col,:date))
    #     a = waread(a)
    #     b = waread(b)
    #     col = ncol(a)-1
    
    #     a = a[!,Cols(col,:date)]
    #     b = b[!,Cols(col,:date)]  
    
    #     df = mall(a,b)
    #     #qplot(x::Vector{Float64},y::Vector{Float64})
    #     return qplot(df)
    # end

    function cntcolread(x::Vector{Any})
        # Get a list of all files in the directory
        files = x
    
        files = filter(inF->isfile(inF),files)
                            #if isfile(inF)
        file_columns = []
        
        for file in files
            # Open the file for reading
            open(file, "r") do io
                # Read the first line of the file
                line = readline(io)
                # Split the line on tabs
                columns = split(line, '\t')
                # Count the number of columns
                num_columns = length(columns)
                push!(file_columns, (file, num_columns))
            end
        end
        
        # Sort the file_columns array based on the number of columns in descending order
        sorted_files = sort(file_columns, by = x -> x[2], rev = true)
        
        for (file, num_columns) in sorted_files
            printstyled(
                rpad("File: $file",45),
            lpad(" | Columns: $num_columns\n",10),color=:green,bold=true)
        end
        return file_columns
    end

    function wqsum()
        #vw = glob("^wq") #same as
        vw::Vector{String} = filter(f->occursin(Regex("^wq"),f),readdir())
        out::Vector{Float64} = []
        for w in vw
        m = Grep.grep("catchment area",readlines(w))|>first
        m = split(m,":")|>last|>strip
        push!(out,parse.(Float64,String((m))))
        end
        out|>println
        A=sum(out)
        println("total sum is $A km²")
    end

    function wqlen()
        vw = glob("^wq")
        out::Vector{Float64} = []
        for w in vw
        m = Grep.grep("channel length",readlines(w))|>first
        m = split(m,":")|>last|>strip
        push!(out,parse.(Float64,String((m))))
        end
        out|>println
        A=round(sum(out)/1000,digits=2)
        println("total channel length is $A km")
    end

    function wqpand(x::AbstractString)
        """
        wqpand("catchment area")
        wqpand("spec")
        wqpand("max. spec")
        """
        vw::Vector{String} = filter(f->occursin(Regex("^wq"),f),readdir())
        out::Vector{Float64} = []
        for w in vw
            m = Grep.grep(x,readlines(w))|>first
            m = split(m,":")|>last|>strip
            push!(out,parse.(Float64,String((m))))
        end
        u = Grep.grep(x,readlines(vw[end]))|>first
        println("match of ",vw[end]," is: \n$u")
        u = split(u," [")|>last|>n -> split(n,"]")|>first
        out|>println
        A=sum(out)
        println("total sum of $x is $A $u")
    end

    function read_wq(file_path::AbstractString)
        data = CSV.File(file_path,header=false,skipto=24,delim=" ",
        maxwarnings=2,ignorerepeated=true,
        debug=true) |> DataFrame
        col_names = CSV.File(file_path,skipto=21,limit=2,
        maxwarnings=2,ignorerepeated=true,delim=" ",header=21) |> DataFrame
        rename!(data, propertynames(col_names))
        return data    
    end

    # function wqplot(file_path::AbstractString)
    #     data = read_wq(file_path)
    #     p = @df data plot(data[!,1],cols(propertynames(data)[2:end]))
    #     println(describe(data))
    #     return p
    # end

    function kge_rec()
        """
        reads recursively
        """
        #fls = rglob("qout")
        needle = r"qout"
        rootdir = pwd()
        results = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if (occursin(needle,filename)
                    && 
                    !occursin(r"yr|mon|grid|scn"i,filename) && 
                    !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|xml|sh|grd|yrly|eps)$", filename)
                    )
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
        # results = filter(x -> isfile(x) && 
        # !occursin(r"yr|mon|grid|scn"i,x) && 
        # !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|xml|sh|grd|yrly|eps)$", x), results)
        out = []
        for x in results
            println("reading $x ...")
            try
                dd = CSV.read(x,DataFrame,
                missingstring="-9999",
                #maxwarnings=1,
                silencewarnings=true,
                ignorerepeated=true,
                delim="\t")
            push!(out,Dict("name"=>x,"KGE"=>kge2(dd)))
            catch
                @warn("$x can't be loaded as a DataFrame, skipping ...")
                continue
            end
            
            
        end
        df = DataFrame(out)
        df.KGE .= replace(df.KGE, NaN => missing)
        dropmissing!(df)
        return sort(df,:KGE;rev = true)
    end

    function grep_with_context(pattern, filename, context)
        lines = readlines(filename)   
        for (i, line) in enumerate(lines)
            if occursin(pattern, line)
                    #println("$(filename):$(i): $(line)")
                    for j in max(1, i - context):min(length(lines), i + context)
                        #printstyled("$(filename):$j: $(lines[j])\n",color=:red)
                        printstyled("$(filename):$j:",color=:red)
                        printstyled("$(lines[j])\n",color=:green)
                    end
                    println("-" ^ 80)
            end
        end
    end
        
    function grep_files(pattern,file_paths, context)
        for file_path in file_paths
            grep_with_context(pattern, file_path, context)
        end
    end

    function dfonly(x::Vector{Any})
        
        z = try filter(f->!occursin(r"^wq|xml|nc|png|svg|jpg",f),x)
        catch
            @warn "vector should be of type string.."
            return
        end
        # z = filter(file -> occursin(x,file), 
        # readdir()[broadcast(x->!occursin(r"^wq|xml|nc|png|svg|jpg",x),readdir())]);
        return(z)
    end

    function rmeq_rec(; rootdir = ".")
        """
        removes empty TS recursively; 
        use with caution!
        """
        
        ext_regex = r".R|.py|.jl|.tex|.pl|.sh|.csv|.html|.xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt"i
        ms = ["-9999", "lin", "log", "--"]
        
        files::Vector{String} = []       #String[]
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if !occursin(ext_regex, filename)
                    push!(files, joinpath(looproot, filename))
                end
            end
        end
        
        for inF in files
            if isfile(inF)
                println("reading ", inF, "...")
                df = CSV.File(inF; delim = "\t", header = 1,
                              silencewarnings = true, 
                              normalizenames = false, 
                              missingstring = ms, 
                              types = Float64) |> DataFrame
                if (isempty(df) || nrow(dropmissing(df,ncol(df)))==0 || ncol(df)==0) 
                    rm(inF)
                    println(basename(inF), " removed!")
                end
            end
        end
    end

    function stp(fn::String)
        fl = CSV.read(fn,DataFrame;limit=4)
        xc = fl[2,5:end]|>collect
        yc = fl[3,5:end]|>collect
        pts = ArchGDAL.IGeometry[]
        for i in 1:length(xc)
            pt = ArchGDAL.createpoint([xc[i],yc[i]])
            #pt = AG.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
            push!(pts,pt)
        end
        nd = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
        return nd
    end

    function ctsum(xx, filename)
        """
        usage: 
        
        nd = ctsum("thickness",infile)
        nd = ctsum("ksat",infile)
        nd = ctsum("Par_n",infile)
        nd = ctsum("theta_sat",infile)
    
        """
        #M = regand("^[0-9]+.*MultipleHorizons",xx)
        #xx="thickness"
        M = Regex("^[0-9]+.*MultipleHorizons+.*"*xx)
        #data = Grep.grep(r"^[0-9]+.*MultipleHorizons+.*thickness", readlines(infile))
        data = Grep.grep(M, readlines(filename))
        # Initialize variables to store the control file information
        controlfile_info = "controlfile: $filename"
        no = Int64[]
        lck = Float64[]
        for ln in data
            num,flds = split(ln," {")
            num = parse.(Int64,strip(num))
            push!(no,num)
            #tck = Grep.grep(r"thickness",split(flds,";"))
            tck = Grep.grep(Regex(xx),split(flds,";"))
            ts = split(tck[1],"=")[2]
            #ts = parse.(Float64,split(ts," ")[2:end])
            #length.(ts)
            ts = split(ts," ")    #strip also possible
            ts = filter(x->x!="",ts)
            ts = parse.(Float64,ts)
            push!(lck,sum(ts))
        end
    
        df = DataFrame(bfid=no,res=lck)
        # Print the control file information and the extracted values
        println(controlfile_info)
        return df
    end

    function ctov(df::DataFrame)
        """
        col_names, col_vectors = ctov(df)
        a,b = cmk.ctov(df[!,Not(:date)]|>dropmissing)
        typeof(a) : Vector{Any}
        typeof(b) : Vector{Float64}
        
        rainclouds(a,b)
        """
        df = df[!,Not(Cols(r"date|year|time"))]|>dropmissing   
        nr = DataFrames.nrow(df)
        
        col_names = String[]
        # for str in names(df)
        #     a = [str for i in 1:nr]
        #     push!(col_names, a)
        # end
    
        for str in names(df)
            a = [str for i in 1:nr]
            append!(col_names, a)
        end
    
        col_vectors = Float64[]
    
        for col in eachcol(df)
            append!(col_vectors, col)
        end
        return col_names, col_vectors
    end

    #https://docs.makie.org/stable/examples/plotting_functions/rainclouds/
    
    function cloudplot(df::DataFrame)

        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw" "
        	end
        

        a,b = ctov(df)
        # colors = Makie.wong_colors()
        # Int.(size(unique(a)))[1]
        colors = Makie.wong_colors()
        # Get unique values
        unique_values = unique(a)
        # Map unique values to colors
        value_to_color = Dict(unique_values[i] => colors[i % length(colors) + 1] for i in 1:length(unique_values))
        # Map values in a to colors
        selected_colors = [value_to_color[value] for value in a]
   
        rainclouds(a,b;
        xlabel = "Basins",
        ylabel = " ", title = ti,
        plot_boxplots = true, cloud_width=0.5, 
        side = :right, violin_limits = extrema,
        #clouds=hist,
        color = selected_colors)



    end
    

end 

# ╔═╡ dc050c7c-c861-4224-bb99-9aa78224e839
#s = "D:/Wasim/regio/out/rc200/x22/f1/tout"
#s = "D:/Wasim/regio/out/e0/Wolfsmuenster-qoutjl"
#s = "D:/Wasim/regio/out/rc200/x22/f1/allrad"
#s = "D:/Wasim/regio/out/rc200/x22/f1/so_Soil_Temperature_Stack.x22.2017"
#s = "D:/Wasim/regio/out/rc200/x22/spin/so_Soil_Temperature_Stack.x22.2017"
#s = "D:/Wasim/regio/out/rc200/x3/c3/so_Soil_Temperature_Stack.x3.2017"
s = "D:/Wasim/regio/out/rc200/x3/loc2/so_Soil_Temperature_Stack.x3.2017"

# ╔═╡ 60464c94-c11c-4289-86ee-7e634dc4d5e9
begin
	df = cmk.waread(s)
	for x in names(df)
	    if startswith(x, "S")
	        parts = split(x, "_")
	        layer_number = parse(Int, parts[end])
	        new_string = "Soil Temperature Layer $layer_number"
	        rename!(df, Symbol(x) => Symbol(new_string))
	    end
	end
	describe(df)
end

# ╔═╡ 0916da1d-f20d-4c4c-9272-9fdb7001a0a2
begin
	#π
	dates = df.date
	#
	vals = select(df,1)|>Matrix|>vec
	#tempo = string.(timestamp(ta))
	tempo = string.(dates)
	lentime = nrow(df)
	slice_dates = range(1, lentime, step=lentime ÷ 8)
	
	fig = Figure(resolution=(600, 400), fonts=(;regular = "consolas"))

	tit = DataFrames.metadata(df)|>only|>last|>basename

	ax3 = Axis(fig[1, 1], title = replace(tit,r"_|so"=>" "),
		#title = tit*" lyr 1",
		xlabel="Date", ylabel=first(names(df)))

	vals2 = select(df,2)|>Matrix|>vec
	line1 = lines!(ax3, 1:lentime, vals; color=:black, linewidth=0.85)
	line2 = lines!(ax3, 1:lentime, vals2; color=:red, linewidth=0.85)
	#
	ax3.xticks = (slice_dates, tempo[slice_dates])
	ax3.xticklabelrotation = π / 4
	ax3.xticklabelalign = (:right, :center)
	fig
end

# ╔═╡ 333d056e-671e-4b7f-87c6-d08eca687775
cmk.cloudplot(df)

# ╔═╡ 0f8ba3e2-2810-4d83-8a05-287cc9f8b721
begin
	        a,b = cmk.ctov(df)
	        colors = Makie.wong_colors()
	        # Get unique values
	        unique_values = unique(a)
	        # Map unique values to colors
	        value_to_color = Dict(unique_values[i] => colors[i % length(colors) + 1] for i in 1:length(unique_values))
	        # Map values in a to colors
	        selected_colors = [value_to_color[value] for value in a]
	   
	        rainclouds(a,b;
	        xlabel = "Categories of Distributions",
	        ylabel = "[°C]", title = DataFrames.metadata(df)|>only|>last|>basename,
	        plot_boxplots = true, cloud_width=0.5, 
	        side = :right, clouds=hist,
	        color = selected_colors)
end

# ╔═╡ 1f798c20-61e5-406d-9c72-7d7c7500af43
begin
	dropmissing!(df)	
	x = cmk.tovec(df,1)|>unique|> x-> x[1:1000]
	y = cmk.tovec(df,2)|>unique|> x-> x[1:1000]
	z = cmk.tovec(df,ncol(df)-1)|>unique|> x-> x[1:1000]
end

# ╔═╡ 78216131-ccf1-4be2-8d33-48ad98312ac7
begin	
	f, ax, tr = tricontourf(x, y, z, mode = :relative, levels = 0.2:0.1:1)
	scatter!(x, y, color = z, strokewidth = 1, strokecolor = :black)
	Colorbar(f[1, 2], tr)
	f
end

# ╔═╡ c307c09b-90f1-40cf-8698-8485c06f656d
begin
	#volcano = readdlm(Makie.assetpath("volcano.csv"), ',', Float64)
	#lnk="D:/Wasim/regio/rcm200/v11/rcm.art-bfid"
	lnk="D:/Wasim/regio/rcm200/v11/rcm.dhm"
	#lnk = "D:/Wasim/regio/rcm3/rcm.dhm"
	#lnk="D:/Wasim/regio/rcm200/v11/rcm.slp"
	volcano = readdlm(lnk, ' ';use_mmap=true, skipstart=7,skipblanks = true)
	#volcano = readdlm(lnk, '\t';use_mmap=true, skipstart=7,skipblanks = true)
	
	replace!(volcano,"" => NaN64,-9999 => NaN64)
	volcano = permutedims(volcano)
	#volcano = transpose(volcano) #|>permutedims
	#volcano = transpose(volcano)
	#reverse!(volcano, dims=(2, 1))
	reverse!(volcano, dims=(1, 2))
end

# ╔═╡ ca1f5526-e217-40be-a940-e4a7deef41e5
begin
	f2 = Figure(resolution = (800, 400))
	
	Axis(f2[1, 1], title = "Relative mode, drop lowest 30%")
	contourf!(volcano, levels = 0.3:0.1:1, mode = :relative)
	
	Axis(f2[1, 2], title = "Normal mode")
	contourf!(volcano, levels = 10)
	
	f2
end

# ╔═╡ cd777948-8555-4ae0-9633-b269e36ba9ef
begin
	#	str=raw"C:\Users\chs72fw\seaborn-data\diamonds.csv"
	#	kd = CSV.read(str,DataFrame)
	cd(dirname(s))
	str=r"Wolf"
	kd = cmk.waread(str)
	describe(kd)
end

# ╔═╡ 1b9e2279-0834-4efa-8c83-3a17e5bd8cac
cmk.cloudplot(kd)

# ╔═╡ 4d2bf624-ce49-4db4-9621-cb448bce2929
# ╠═╡ disabled = true
#=╠═╡
begin
	ar = volcano
	n1,n2 = size(ar)
	zv = vec(volcano)
	n1 = [x for x in 1:(n1)] 
	n2 = [x for x in 1:(n2)] 
		
end
  ╠═╡ =#

# ╔═╡ 96e302da-3cd9-4221-8e63-45c84f396c52
#=╠═╡
begin
	lvl = 10.0.^range(0.3, 3.5; length=10)
	cmap = Makie.sampler(:hsv, 100; scaling=Makie.Scaling(x2 -> x2^(1 / 10), nothing))
	f0, ax0, ct0 = contour(n1, n2, zv; labels=true, lvl, cmap)
	f0
end
  ╠═╡ =#

# ╔═╡ 823dd256-a9e7-46a3-89bc-d3446845d471
# ╠═╡ disabled = true
#=╠═╡
begin
	himmelblau(x, y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
	x2 = y2 = range(-6, 6; length=100)
	z2 = himmelblau.(x2, y2')
	
	levels = 10.0.^range(0.3, 3.5; length=10)
	colormap = Makie.sampler(:hsv, 100; scaling=Makie.Scaling(x2 -> x2^(1 / 10), nothing))
	f3, ax2, ct2 = contour(x2, y2, z2; labels=true, levels, colormap)
	f3
end
  ╠═╡ =#

# ╔═╡ 9f2b3985-743a-4ff8-95ad-6d6a0c2c3009
# ╠═╡ disabled = true
#=╠═╡
begin
	# Calculate rolling mean and standard deviation
	dn = kd
	#rename!(dn, :Wolfsmuenster => :value)
	window_size = 30
	rolling_mean = [mean(dn.value[i-window_size+1:i]) for i in window_size:length(dn.date)]
	rolling_std = [std(dn.value[i-window_size+1:i]) for i in window_size:length(dn.date)]
	
	# Plot time series with rolling mean and standard deviation
	fig = Figure(resolution = (800, 400))
	axm = Axis(fig[1,1])
	c,d = cmk.ctov(dn)
	#lines!(axm, c,d, color = :blue)
	#lines!(axm, dn.date[window_size:end], rolling_mean, color = :red)
	#band!(axm, dn.date[window_size:end], rolling_mean .- dn.rolling_std, rolling_mean .+ rolling_std, color = (:red, 0.2))
	#xlabel!(axm, "Date")
	#ylabel!(axm, "Value")
	#title!(axm, "Time Series Analysis")
	#fig
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Grep = "775960c6-9b90-5df0-b405-1e337feb71e5"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.11"
CairoMakie = "~0.10.8"
DataFrames = "~1.6.1"
DelimitedFiles = "~1.9.1"
Distributions = "~0.25.100"
Grep = "~0.2.0"
Makie = "~0.19.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "39ad9906cec4607eb533d5047381ff92079a7d33"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractLattices]]
git-tree-sha1 = "f35684b7349da49fcc8a9e520e30e45dbb077166"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.2.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "ef9997b3d5547c48b41c7bd8899e812a917b409d"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.4"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "44dbf560808d49041989b8a96cae4cffbeb7966a"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.11"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools", "SHA"]
git-tree-sha1 = "30562a68ded3dabe80109caf6b4de73a48ac27bc"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.8"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "d9a8f86737b665e15a9641ecbac64deef9ce6724"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.23.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "e460f044ca8b99be31d35fe54fc33a5c33dd8ed7"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.9.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fe2838a593b5f776e1597e086dcd47560d94e816"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.3"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "a1d8532de83f8ce964235eff1edeff9581144d02"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.7.2"
weakdeps = ["MakieCore"]

    [deps.DelaunayTriangulation.extensions]
    DelaunayTriangulationMakieCoreExt = "MakieCore"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "938fe2981db009f531b6332e31c58e9584a2f9bd"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.100"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArraysCore", "Test"]
git-tree-sha1 = "276e83bc8b21589b79303b9985c321024ffdf59c"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.5"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "048dd3d82558759476cff9cff999219216932a08"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.6.0"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "bb198ff907228523f3dee1070ceee63b9359b6ab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "424a5a6ce7c5d97cca7bcc4eac551b97294c54af"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.9"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grep]]
deps = ["Test"]
git-tree-sha1 = "011aa67826e8b02bbd82c13946421bb563836128"
uuid = "775960c6-9b90-5df0-b405-1e337feb71e5"
version = "0.2.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "f57a64794b336d4990d90f80b147474b869b1bc4"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "FastRounding", "LinearAlgebra", "Markdown", "Random", "RecipesBase", "RoundingEmulator", "SetRounding", "StaticArrays"]
git-tree-sha1 = "5ab7744289be503d76a944784bac3f2df7b809af"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.20.9"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "8e59ea773deee525c99a8018409f64f19fb719e6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.7"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "327713faef2a3e5c80f96bf38d1fa26f7a6ae29e"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "90442c50e202a5cdf21a7899c66b240fdef14035"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.7"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Permutations", "Primes", "SimplePolynomials"]
git-tree-sha1 = "558a338f1eeabe933f9c2d4052aa7c2c707c3d52"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.1.12"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Setfield", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "e81675589ba7199a82443e87fc52e17eeceac2e8"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.8"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "f56b09c8b964919373d61750c6d8d4d2c602a2be"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.5"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "8f52dbaa1351ce4cb847d95568cb29e62a307d93"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.6"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mods]]
git-tree-sha1 = "61be59e4daffff43a8cec04b5e0dc773cbb5db3a"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "1.3.3"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.Multisets]]
git-tree-sha1 = "8d852646862c96e226367ad10c8af56099b4047e"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e78db7bd5c26fc5a6911b50a47ee302219157ea8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.10+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "963b004d15216f8129f6c0f7d187efa136570be0"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.7"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "9b02b27ac477cad98114584ff964e3052f656a0f"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.0"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "6e6cab1c54ae2382bcc48866b91cf949cea703a1"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.16"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield"]
git-tree-sha1 = "6ded5b759921314670b726dc6ce479675046bc04"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "ee094908d720185ddbdc58dbe0c1cbe35453ec7a"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.7"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "9712ebc42e91850f35272b48eb840e60c0270ec0"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.7"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "04bdff0b09c65ff3e06a05e3eb7b120223da3d39"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "0d15c3e7b2003f4451714f08ffec2b77badc2dc4"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.3.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "b608903049d11cc557c45e03b3a53e9260579c19"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.4"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "dcc02923a53f316ab97da8ef3136e80b4543dbf1"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.0"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "d537c31cf9995236166e3e9afc424a5a1c59ff9d"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.14"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "3a6fb395e37afab81aeea85bae48a4db5cd7244a"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["CRC32c", "Compat", "Dates", "SHA", "Tables", "TupleTools", "UUIDs"]
git-tree-sha1 = "0b8b801b8f03a329a4e86b44c5e8a7d7f4fe10a3"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "0.3.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "9cabadf6e7cd2349b6cf49f1915ad2028d65e881"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.2"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "8621f5c499a8aa4aa970b1ae381aae0ef1576966"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.4"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═cdb4bd40-4199-11ee-1eb9-19c27eed6305
# ╟─478b18ee-c9b5-46f0-abd5-fa22a8731bbf
# ╠═dc050c7c-c861-4224-bb99-9aa78224e839
# ╠═60464c94-c11c-4289-86ee-7e634dc4d5e9
# ╠═0916da1d-f20d-4c4c-9272-9fdb7001a0a2
# ╠═333d056e-671e-4b7f-87c6-d08eca687775
# ╠═0f8ba3e2-2810-4d83-8a05-287cc9f8b721
# ╠═1f798c20-61e5-406d-9c72-7d7c7500af43
# ╠═78216131-ccf1-4be2-8d33-48ad98312ac7
# ╠═c307c09b-90f1-40cf-8698-8485c06f656d
# ╠═ca1f5526-e217-40be-a940-e4a7deef41e5
# ╠═cd777948-8555-4ae0-9633-b269e36ba9ef
# ╠═1b9e2279-0834-4efa-8c83-3a17e5bd8cac
# ╠═4d2bf624-ce49-4db4-9621-cb448bce2929
# ╠═96e302da-3cd9-4221-8e63-45c84f396c52
# ╠═823dd256-a9e7-46a3-89bc-d3446845d471
# ╠═9f2b3985-743a-4ff8-95ad-6d6a0c2c3009
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
