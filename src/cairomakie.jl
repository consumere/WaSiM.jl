# adapt cairomakie instead of plots 
if !(isdefined(Main, :src_path) && ispath(src_path))
    include("smallfuncs.jl")
end

module cmk
    using Makie, CairoMakie
    import GeoMakie
    using DataFrames, CSV, Statistics, Dates, Distributions
    using DelimitedFiles, Grep, Printf
    # using PrettyTables, SHA
    # using PyCall #for old mkheat
    using Rasters
    using LaTeXStrings
    import NCDatasets
    import ArchGDAL
    import GeoInterface

    #very useful
    CairoMakie.set_theme!(theme_latexfonts())
    # p=Rasters.rplot
    # p(A;title="TIT",colormap = :thermal, ylabel="ylab",
        #plottype=Contourf,colorbar_position=Makie.Bottom())

    #rasters
    export mkecont, mkecont2, mkemon, mkestream, 
    mkpoly, mkrheat, mkh, mkfzs, gmwas
    #timeseries
    export dfpl, ctov, cloudplot, cloudplot2, 
    tsbar, tsdf, tsp, tsp2, tsp3, tspblack

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

    """
    gets sorted DF by size recursivley
    """
    function fz()
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

    """
    enhanced Reader.
    Read the text file, preserve line 1 as header column \n
    `function waread(x::String;station=false,proj=false,pts_and_data=false,src=EPSG(25832),dst=EPSG(4326))`
    """
    function waread(x::String;station=false,proj=false,pts_and_data=false,src=EPSG(25832),dst=EPSG(4326))
        ms = ["-9999","lin","log","--"]
        if station
            fl = CSV.read(x,DataFrame;limit=4)
            ez = fl[1,5:end]|>collect
            xc = fl[2,5:end]|>collect
            yc = fl[3,5:end]|>collect
            no = fl[4,5:end]|>collect
            pts = ArchGDAL.IGeometry[]
            
            #for i in 1:length(xc)
            for i in 1:lastindex(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                if proj
                    pt = ArchGDAL.reproject(pt,src,dst)
                end
                push!(pts,pt)
            end
            
            nd = DataFrame(geometry=pts, 
                name=names(fl)[5:end],
                ez = ez, no = no,
                xc=xc, yc=yc)

            if !pts_and_data
                return nd
            end
        end

        df = CSV.read(x, DataFrame; 
            delim="\t", header=1, 
            missingstring=ms, 
            maxwarnings = 1, #silencewarnings = true,
            normalizenames=true, types=Float64)
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
            
        if (pts_and_data && !station)
            fl = CSV.read(x, DataFrame; limit=4)
            ez, xc, yc, no = [fl[i, 5:end] |> collect for i in 1:4]
            pts = [ArchGDAL.createpoint([xc[i], yc[i]]) for i in 1:lastindex(xc)]
            if proj
                pts = ArchGDAL.reproject.(pts, Ref(src), Ref(dst))
            end
            nd2 = DataFrame(geometry=pts, name=names(fl)[5:end], ez=ez, no=no, xc=xc, yc=yc)
            return (nd2,df)
        elseif (pts_and_data && station)
            return (nd,df)
        else
            return df
        end
        
    end

    function waread(x::AbstractString)
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


    """
    Read the text file, preserve line 1 as header column
    Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
    This avoids reading the entire file into memory at once, 
    which can be more memory-efficient for large datasets.
    kwargs are passed to CSV.read
    """
    function waread2(x::String;kwargs...)
        ms = ["-9999","-9999.0","lin", "log", "--"]
        df = CSV.File(x; delim="\t", header=1, normalizenames=true, 
            missingstring=ms, types=Float64, kwargs...) |> DataFrame
        dropmissing!(df,1)
        dt2 = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
        select!(df, Not(1:4))
        df.date = dt2
        metadata!(df, "filename", x, style=:note)
        return df
    end


    """
    waread2 on regex
    kwargs... passed to CSV.File
    Read the text file, preserve line 1 as header column
    Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
    This avoids reading the entire file into memory at once, 
    which can be more memory-efficient for large datasets.
    """
    function waread2(x::Regex;kwargs...)
        inF = first(filter(x->!occursin(r"yrly|nc|png|svg|grd",x),Grep.grep(x,readdir())))

        ms = ["-9999", "lin", "log", "--"]
        df = CSV.File(inF; delim="\t", header=1, 
            normalizenames=true, missingstring=ms,
            #maxwarnings=2, 
            silencewarnings=true, 
            types=Float64, kwargs...) |> DataFrame
        dropmissing!(df,1)        
        df.date = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
        select!(df, Not(1:4)) #since date is at last position
        metadata!(df, "filename", x, style=:note)
        return df
    end

    """
    Fastest Reader. is also dfr.
    Read the text file, preserve line 1 as header column
    """
    function waread(x::String)
        ms = ["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame; 
            delim="\t", header=1, missingstring=ms, 
            maxwarnings = 1, #silencewarnings = true,
            normalizenames=true, types=Float64)
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

    """
    Read the text file, preserve line 1 as header column
    """
    function waread(x::Regex)
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
        dfs = filter(x->size(x,1)>0,dfs)
        return(dfs)
    end

    function loadalldfs(path::Vector{Any})
        files = path
        dfs = DataFrame[]
        for file in files
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
            file_path = file
        println("reading ",file_path,"...")
        p1 = waread(file_path)
        push!(dfs, p1)
            end
        end
        dfs = filter(x->size(x,1)>0,dfs)
        return(dfs)
    end

    function loadalldfs(files::Vector{String})
        dfs = DataFrame[]
        for file in files
            if isfile(file) && (!occursin(r"tex|xml|fzt|ftz|log|wq|txt|yrly|nc|png|svg",file))
            file_path = file
            println("reading ",file_path,"...")
            p1 = waread(file_path)
            push!(dfs, p1)
            end
        end
        dfs = filter(x->size(x,1)>0,dfs)
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
        p1 = waread(file_path)
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
        p1 = waread(file_path)
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
        p1 = waread(file_path)
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
        #p1 = waread(file_path)
        #push!(dfs, p1)
        push!(nms, file)
            end
        end
        #return(dfs)
        return(nms)
    end

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
    #     df=waread(df)
    #     nm = propertynames(df)[1:end-1];
    #     o = collect(DataFrames.metadata(df))[1][2] |>basename
    #     @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=o)     
    # end

    # function lplotf(df::String)
    #     df=waread(df)
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
        p1 = waread(file_path)
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
        #p1 = waread(file_path)
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
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
    end



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


    function rglob(prefix::String)
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
        return string.(results)
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
        return string.(results)
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



    function vjl(regex)
        # greps jl from julia folder
        pt = src_path;
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
    
"""
    reads, reduces + merges by date
    ds = innerjoin(unique.(dfs, xcol)..., on = xcol, makeunique=true)
    example:
    df = mall(glob(r"qges|qbas|qd"))
    """
    function mall(files::Vector{String};xcol=:date)
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
        # df = reduce((left, right) -> 
        # innerjoin(left, right, on = xcol,makeunique=true), 
        # dfs)
        df = innerjoin(unique.(dfs, xcol)..., on = xcol, makeunique=true)
        #xcol = string(xcol)
        #df = hcat(df[!,Not(Cols(Regex(xcol)))],df[:,Cols(Regex(xcol))])
        df = hcat(df[!,Not(Cols(xcol))],df[:,Cols(xcol)])
        return(df)
    end

    """
    reduces + merges by date
    files::Vector{DataFrame};xcol=:date
    """
    function mall(files::Vector{DataFrame};xcol=:date)
        # df = reduce((left, right) -> 
        # innerjoin(left, right, on = xcol,
        #     makeunique=true), files)
        df = innerjoin(unique.(files, xcol)..., 
            on = xcol, makeunique=true)
        df = hcat(df[!,Not(Cols(xcol))],df[:,Cols(xcol)])
        return(df)
    end
    
    """
    reduces + merges by date
    files::Vector{Any};xcol=:date
    """
    function mall(files::Vector{Any};xcol=:date)
        # df = reduce((left, right) -> 
        # innerjoin(left, right, on = xcol,makeunique=true), 
        # files)
        df = innerjoin(unique.(files, xcol)..., 
            on = xcol, makeunique=true)
        df = hcat(df[!,Not(Cols(xcol))],df[:,Cols(xcol)])
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
        first |>waread
            
        return(df)
    end
 
    function yrsum(x::String)
        df = waread(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(DataFrames.groupby(df, :year), y .=> sum .=> y);
        return(df_yearsum)
    end

    function yrsum(x::Regex)
        df = waread(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(DataFrames.groupby(df, :year), y .=> sum .=> y);
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
        df_yearsum = DataFrames.combine(DataFrames.groupby(df, :year), y .=> sum .=> y);
        return(df_yearsum)
    end

    function yrmean(x::String)
        df = waread(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(DataFrames.groupby(df, :year), y .=> mean .=> y);
        return(df_yearsum)
    end

    function yrmean(x::Regex)
        df = waread(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(DataFrames.groupby(df, :year), y .=> mean .=> y);
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
        df_yearsum = DataFrames.combine(DataFrames.groupby(df, :year), y .=> mean .=> y);
        return(df_yearsum)
    end

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
        """
        df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),
            dfs))] |> first
        return(df)
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
        df_monthsum = DataFrames.combine(DataFrames.groupby(df, :month), y .=> sum .=> y);
        return(df_monthsum)
    end

    function monsum(x::DataFrame)
        df = x
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(DataFrames.groupby(df, :month), y .=> sum .=> y);
        return(df_monthsum)
    end

    function monmean(x::String)
        df = waread(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(DataFrames.groupby(df, :month), y .=> mean .=> y);
        return(df_monthsum)
    end

    function monmean(x::DataFrame)
        df = x
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :month] = month.(df[!,:date]);
        df_monthsum = DataFrames.combine(DataFrames.groupby(df, :month), y .=> mean .=> y);
        return(df_monthsum)
    end


    function vg2(regex::AbstractString, ending::AbstractString)
        cmd = `grep --color=always -C2 -rIHn -E "$regex" --include="*.$ending"`
        run(cmd)
    end

    "--- reader with drop exept of first col ---"
    function so_read(x::AbstractString)
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


    "--- main reader ---"
    function readmhm(x::AbstractString)
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

    """
    greps from repl_history
    """
    function vgrepl(snippet::AbstractString)
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

    """
    greps from current dir iRegex
    """
    function glob(x::AbstractString)
        filter(file -> occursin(Regex(x,"i"),file), readdir())
    end

    """
    greps from current dir Regex
    """
    function glob(x::Regex)
        filter(file -> occursin(x,file), readdir())
    end

    """
    date to last position
    """
    function reorder_df(df::DataFrame)
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
                dout.nm .= file
                push!(outdf,dout)
            catch
                @warn("error! files that can't be loaded as a DataFrame")
                # Skip files that can't be loaded as a DataFrame
                continue
            end
        end
        return(outdf)
    end

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

    """
    --- reader with fewer constrains ---
    no |> dropmissing 
    df[!,Cols(r"^Col|date")] |>dfp  
    """
    function odfr(x::AbstractString)
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

    """
    --- reader with fewer constrains ---
    with |> dropmissing on 2nd col
    df[!,Cols(r"^Col|date")] |>dfp  
    """
    function odfr(x::Regex)
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

    
    """
    removes " inside textfile
    """
    function process_file2(input_file::String)
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
        dy = DataFrames.combine(DataFrames.groupby(df, :year), y .=> sum .=> y);
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

    """
    only for wasim output files
    """
    function fread(z::Union{Regex,String})
        if isa(z,Regex)
            v = filter(file -> occursin(z,file), readdir());
            z = v[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt|log"i),v)]|>first
        end 
        println("loading $z ...")   
        m = map(x->(x=>Int64),[:YY,:MM,:DD,:HH])
        ms = ["-9999","-9999.0","lin", "log", "--"]
        df = CSV.read(z,DataFrame;
            delim="\t", header=1, 
            normalizenames=true, 
            missingstring=ms,
            types=Dict(m),
            stripwhitespace=true)|>f->dropmissing(f,:YY)
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
        select!(df,Not(1:4))
        DataFrames.metadata!(df, "filename", z, style=:note)
    end

    # readall = loadalldfs
    # readmeteo = waread
    # loaddf = waread

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
        include("C:/Users/Public/Documents/Python_Scripts/julia/smallfuncs.jl")
    end
 
    function kge_fread()
        kge_read(pwd(),"outjl");
    end

 
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

    function op()
        #pwrs""" explorer . """
        open(`powershell -noprofile explorer . `,"w",stdout)
    end 

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

    """
    map(typeof, eachcol(df)) #check types of cols
    msk = broadcast(x->typeof(x)==Vector{Float64},df)
    """
    function subset_dataframe_by_mask(df::DataFrame, msk::DataFrame)
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

    # function pydf_to_julia(py_df::PyObject)
    #     """
    #     no transposing
    #     """
    #     # Convert each column of the Python DataFrame to a Julia array
    #     col_names = py_df.columns  # Get the column names from the Python DataFrame
    #     col_arrays = [convert(Array, py_df[col]) for col in col_names]
    #     # Create a Julia DataFrame using the converted arrays and column names
    #     julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, col_arrays))    
    #     return julia_df
    # end 
    
    # function pydf(py_df::PyObject)
    #     """
    #     Convert each column of the Python DataFrame to a Julia array
    #     """
    #     # fn = filename
    #     # pyo = py"""waread3($fn).reset_index(drop=False)"""
    #     # pdf = wa.pydf(pyo)
    #     # names(pdf)
    #     # dfp(pdf)
    #     col_names = py_df.columns  # Get the column names from the Python DataFrame
    #     col_names = convert(Array, col_names)
    #     col_arrays = [convert(Array, py_df[col]) for col in col_names]
    #     jdf = DataFrame(col_arrays, :auto)
    #     #size(jdf)
    #     fn = try
    #         py_df.filename
    #     catch
    #         @info "no filename present"
    #     end

    #     metadata!(jdf, "filename", fn, style=:note);
    #     rename!(jdf, col_names);
    #     return jdf
    # end

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

    """
    removes empty TS recursively; 
    use with caution!
    """
    function rmeq_rec(; rootdir = ".")
        
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

    """
    runs inside REPL
    """
    function runwasim(ctlfile)
        try
            exec = normpath("C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05/wasimvzo64.exe")
            run(`$exec $ctlfile`)
        catch e
            println("no valid input!")
        end
    end

    """
    usage: 
    
    nd = ctsum("thickness",infile)
    nd = ctsum("ksat",infile)
    nd = ctsum("Par_n",infile)
    nd = ctsum("theta_sat",infile)

    """
    function ctsum(xx, filename)
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

    """
    returns col_names, col_vectors from a DataFrame
    col_names, col_vectors = ctov(df)
    a,b = cmk.ctov(df[!,Not(:date)]|>dropmissing)
    typeof(a) : Vector{Any}
    typeof(b) : Vector{Float64}
    rainclouds(a,b)
    """
    function ctov(df::DataFrame)
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


    """
    cloudplot(df::DataFrame)
    #https://docs.makie.org/stable/examples/plotting_functions/rainclouds/
    """
    function cloudplot(df::DataFrame)

        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
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
        #category_labels = names(df[!,Not(:date)]),
        fig = Figure() 
        ax3 = Axis(fig[1, 1], 
            title = replace(ti,r"_|so"=>" "),
            xlabel="", ylabel = "") 

        rainclouds!(a,b;            
            orientation = :horizontal,
            plot_boxplots = true, 
            cloud_width=2.5, 
            side = :right, 
            violin_limits = extrema,
            #clouds=hist,
            color = selected_colors)
        hideydecorations!(ax3, grid=true, ticks=true)
        hidexdecorations!(ax3, grid=true, ticks=true,  ticklabels = false)
        ax3.xticks = (1:ncol(df),names(df))
        #ax3.yticks = #repeat("x",3)
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)

        return fig

        #cbs = names(df[!,Not(:date)])
        #CairoMakie.yticks!(;ytickrange=1:length(cbs),yticklabels=cbs)

    end

    """
    cloudplot2(df::DataFrame)
    with less options
    """
    function cloudplot2(df::DataFrame)
    
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
        end
    
        a,b = ctov(df)
        colors = Makie.wong_colors()
        unique_values = unique(a)
        value_to_color = Dict(unique_values[i] => colors[i % length(colors) + 1] for i in 1:length(unique_values))
        selected_colors = [value_to_color[value] for value in a]

        fig = Figure( #resolution=(600, 400), 
            #fonts=(;regular = "consolas")            
        )
        
        ax3 = Axis(fig[1, 1], 
            title = replace(ti,r"_|so"=>" "),
                xlabel="", 
                ylabel = "") 
            #ylabel=first(names(df))           )
    
        rainclouds!(a,b;
            #xlabel = "",
            #ylabel = " ", 
            #title = ti,
            gap=0.0002,
            #cloud_width = maximum(b),
            #clouds = violin,
            #violin_limits = extrema(b).*.8,
            #xticklabelalign = (:right, :center),
            color = selected_colors)
            #clouds=hist,
            #xticklabelrotation = π / 4,
        ax3.xticks = (1:ncol(df),names(df))
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        return fig
    end
    
    """
    plots timeseries
    """
    function tspblack(df::DataFrame)
        dt = df.date
        fig = Figure( #resolution=(600, 400), 
            fonts=(;regular = "consolas")            
            )
        
        tempo = string.(dt)
        #lentime = nrow(df)
        lentime = size(df,1)
        slice_dates = range(1, lentime, step=lentime ÷ 8)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        
        ax3 = Axis(fig[1, 1], 
            title = replace(tit,r"_|so"=>" "),
                xlabel="Date", 
                #yscale = :log10,
                ylabel=first(names(df)))
    
        #cols = propertynames(df)  
        cols = names(df)  
        filter!(x->!occursin(r"date|year",x),cols)

        # colors = Makie.wong_colors()
        # # Map unique values to colors
        # value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # # Map values in a to colors
        # selected_colors = [value_to_color[value] for value in cols]
        # colors = Makie.wong_colors()
        #v = 1:length(cols)
        #for v, col in enumerate(cols)
        for col in cols
            #vals = select(df,col)|>Matrix|>vec   
            vals = Vector(df[:, col])
            lines!(ax3, 1:lentime, vals; 
            color=:black,         #     color=colors[v], 
             linewidth=0.85)
        end
        
        #vals2 = select(df,2)|>Matrix|>vec
        #line1 = lines!(ax3, 1:lentime, vals; color=:black, linewidth=0.85)
        #line2 = lines!(ax3, 1:lentime, vals2; color=:red, linewidth=0.85)
        #
        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        return fig
    end

    """
    plots timeseries
    """
    function tsp(df::DataFrame)
        dt = df.date
        fig = Figure( #resolution=(600, 400), 
                fonts=(;regular = "consolas")            
                )
    
        
        tempo = string.(dt)
        lentime = size(df, 1)
        slice_dates = range(1, lentime, step=lentime ÷ 8)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        
        ax3 = Axis(fig[1, 1], 
                   title = replace(tit, r"_|so" => " "),
                   xlabel = "Date",
                   ylabel = first(names(df)))
        
        cols = names(df)
        filter!(x -> !occursin(r"date|year", x), cols)
        
        colors = Makie.wong_colors()
        # Map unique values to colors
        value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # Map values in a to colors
        selected_colors = [value_to_color[value] for value in cols]
        #colors = Makie.wong_colors(length(cols))
        
        for (i, col) in enumerate(cols)
            #vals = select(df, col) |> Matrix |> vec
            vals = Vector(df[:, col])
            lines!(ax3, 1:lentime, vals; color = selected_colors[i], linewidth = 0.85)
        end
        
        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        
        return fig
    end

    """
    plots timeseries for each column separately in one figure
    """
    function tsp2(df::DataFrame)
        dt = df.date
        tempo = string.(dt)
        lentime = size(df, 1)
        slice_dates = range(1, lentime, step=lentime ÷ 8)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        cols = filter(x -> !occursin(r"date|year", x), names(df))
    
        fig = Figure( #resolution=(800, 400), 
        fonts=(;regular = "consolas")            
        )
                     
        
        for (i, col) in enumerate(cols)
            ax = Axis(fig[1, i], 
                       title = replace(tit, r"_|so" => " "),
                       xlabel = raw"", #"Date",
                       ylabel = col)
            
            #vals = select(df, col) |> Matrix |> vec
            vals = Vector(df[:, col])
            
            lines!(ax, 1:lentime, vals; linewidth = 0.85)
            
            ax.xticks = (slice_dates, tempo[slice_dates])
            ax.xticklabelrotation = π / 4
            ax.xticklabelalign = (:right, :center)
        end
        
        return fig
    end

    """
    plots timeseries for each column separately in one figure
    """
    function tsp3(df::DataFrame)
        dt = df.date
        tempo = string.(dt)
        lentime = size(df, 1)
        slice_dates = range(1, lentime, step=lentime ÷ 8)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        cols = filter(x -> !occursin(r"date|year", x), names(df))
    
        fig = Figure( #resolution=(800, 800), 
        fonts=(;regular = "consolas")            
        )
                     
        
        for (i, col) in enumerate(cols)
            row, cl = fldmod1(i, 2)
            ax = Axis( fig[row, cl],
                       title = replace(tit, r"_|so" => " "),
                       xlabel = "Date",
                       ylabel = col)
            
            #vals = select(df, col) |> Matrix |> vec
            vals = Vector(df[:, col])
            lines!(ax, 1:lentime, vals; linewidth = 0.85)
            
            ax.xticks = (slice_dates, tempo[slice_dates])
            ax.xticklabelrotation = π / 4
            ax.xticklabelalign = (:right, :center)
        end
        
        return fig
    end

    """
    barplots yearsum timeseries
    baryrsum
    """
    function tsbar(df::DataFrame;kw...)
        df = yrsum(df)
        #dt = df.date
        dt = df.year
        fig = Figure(kw...) 
                #resolution=(600, 400), 
                #fonts=(;regular = "consolas")            
                #)
    
        
        tempo = string.(dt)
        lentime = size(df, 1)
        #slice_dates = range(1, lentime, step=lentime ÷ 8)
        slice_dates = range(1, lentime, step=1)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
                
        
        ax3 = Axis(fig[1, 1], 
                   title = replace(tit, r"_|so" => " "),
                   xlabel = raw" ",
                   ylabel = first(names(df)))
  
        cols = names(df)
        filter!(x -> !occursin(r"date|year", x), cols)
        
        colors = Makie.wong_colors()
        # Map unique values to colors
        value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # Map values in a to colors
        selected_colors = [value_to_color[value] for value in cols]
        
        # dfm = DataFrames.stack(df)
        # dfm.variable = string.(dfm.variable)
        # fig = Figure()
        # ax = Axis(fig[1,1])
        # barplot!(ax,
        #     dfm.year,
        #     dfm.value,
        #     color=dfm.year,
        #     stack=dfm.year )
      


        for (i, col) in enumerate(cols)
            #vals = select(df, col) |> Matrix |> vec
            vals = Vector(df[:, col]) #more efficient
            #lines!(ax3, 1:lentime, vals; color = selected_colors[i], linewidth = 0.85)
            # barplot!(ax3, 1:lentime, vals; 
            # stack = df.year,
            # color = selected_colors[i], alpha = 0.85)
            barplot!(ax3, 1:lentime, vals,
            stack = df.year,
            color = selected_colors[i])
        end
        
        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        
        return fig
    end

    function toMain()
        fnames = names(Main.cmk, all=true)
        for submodule in fnames
            @eval import Main.cmk.$submodule
        end
    end

    # using DataFrames, CSV, Statistics, Dates, Distributions

    """
    find LOG. R-SQUARE > .4 recursivley
    """
    function findlog(;lb=.4)
        v = (;recursive=true)
        #map(x->DataFrames.subset(x,3 => ByRow(<(1))),v)
        k = try 
            map(x->DataFrames.subset(x,3 => ByRow(>(lb))),v)
            catch
                @error "no df on lowerbound! "
                return
        end
        df = try 
            reduce(vcat,k) 
            catch
                @error "vcat failed! "
                return
        end
        
        df = try 
            DataFrames.subset(df,3 => ByRow(<(1.0)))
            catch
                println(first(k))
                @error "no df on upperbound! "
                return
        end
        
        
        return df
    end

    """
    Extract the coordinates from the polygon and plot them
    """
    function mkpoly(geom::Vector{ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}})
        ks = geom|>first
        x = [first(pt) for pt in GeoInterface.coordinates(ks)|>first]
        y = [last(pt) for pt in GeoInterface.coordinates(ks)|>first]
        # fig = Figure()
        # axs = Axis(fig, xlabel="", ylabel="")
        # Makie.lines!(x, y)    
        # return fig
        return Makie.lines(x, y)    
    end

    """
    Contour over heatmap
    """                     
    function mkecont(r::Raster;wasim=true,missval=0)
        if wasim
            z = r.data[:,:,1]
            reverse!(z,dims=1)
            z = z'
            replace!(z, missval=>missing)
        else
            z = r.data
            reverse!(z,dims=1)
            reverse!(z,dims=2)
            replace!(z, missval=>missing)
        end
        #ex = extrema(z|>skipmissing)
        #lscale = ex[1]:10:ex[2] # Adjusted levels
        fig = Figure(
            #size=(800, 600), 
            fontsize=22);
        axs = Axis(fig[1,1], aspect=1, xlabel="x", ylabel="y")
        p1 = heatmap!(axs, z, colormap=(:plasma, 0.87))
        contour!(axs, z; color=:black) #, levels=lscale
        Colorbar(fig[1, 2], p1, width=20, ticksize=20, tickalign=1)
        xmax = findmax(size(z))[1] #gets the index of the max value
        limits!(axs, 1, xmax, 1, xmax)
        # hideydecorations!(axs, grid=true, ticks=false)
        # hidexdecorations!(axs, grid=true, ticks=false)
        return fig
    end


    """
    contourplots from rasters 3 plots
    r::Raster;wasim=true,missval=0
    """
    function mkecont2(r::Raster;wasim=true,missval=0)
        if wasim
            z = r.data[:,:,1]
            reverse!(z,dims=1)
            z = z'
            try
                replace!(z, missval=>missing)    
            catch
                @error "MethodError: Cannot `convert` an object of type Missing to an object of type Float64"
                @warn "replace missing failed!"
            end
            
        else
            z = r.data
            reverse!(z,dims=1)
            reverse!(z,dims=2)
            try
                replace!(z, missval=>missing)    
            catch
                @error "MethodError: Cannot `convert` an object of type Missing to an object of type Float64"
                @warn "replace missing failed!"
            end
        end
        #ex = extrema(z|>skipmissing) #can be errornous due to NaN32
        #lscale = ex[1]:10:ex[2] # Adjusted levels
        fig = Figure(size=(1200, 400), fontsize=22);
        axs = [Axis(fig[1, j], aspect=1, xlabel="x", ylabel=j == 1 ? "y" : "")
                for j in 1:3]
        p1 = heatmap!(axs[1], z, colormap=:plasma)
            contour!(axs[2], z; color=:black) 
            heatmap!(axs[3], z; colormap=(:plasma, 0.85))
        contour!(axs[3], z; color=:white)
        Colorbar(fig[1, 4], p1, width=20, ticksize=20, tickalign=1)
        xmax = findmax(size(z))[1] #gets the index of the max value
        [limits!(axs[i], 1, xmax, 1, xmax) for i in 1:3]
        [hideydecorations!(axs[i], grid=false, ticks=false) for i in 2:3]
        return fig
    end

    function plot_raster_and_points(raster, points)    
        fig = Figure(fontsize=22); #size=(1200, 400),
        #axs = Axis(fig, xlabel="x", ylabel="y"
        fig = Figure();
        axs = Axis(fig, xlabel="", ylabel="",aspect=1)
        Makie.heatmap!(axs,raster)
        for pt in ptc_tuples
            Makie.scatter!(axs,points, markersize = 5, color = :red)
        end
        # Set explicit limits if needed
        #GeoInterface.bbox(ptc_tuples)
        # xmax,ymax = findmax(ptc_tuples)[1]
        # xmin,ymin = findmin(ptc_tuples)[1]
        # limits!(axs, xmin, xmax, ymin, ymax)
        #limits!(axs, 1, size(r)[2], 1,size(r)[1])
        hideydecorations!(axs, grid=false, ticks=true)
        return fig
    end

    """
    lat lon reverse
    """
    function reverse_coords(polygon)
        # Extract the coordinates from the polygon
        # crds = GeoInterface.coordinates.(polygon)
        # # Reverse each pair of coordinates
        # reversed_crds = [[[(lon, lat) for (lat, lon) in ring] for ring in polygon] for polygon in crds]
        
        # coords = GeoInterface.coordinates.(gd.geometry)
        coords = GeoInterface.coordinates(polygon)
        # Flatten the list of points
        #points = reduce(vcat, coords[1])
        ptc = coords[1]
        #ptc = GeoInterface.coordinates.(points)
        ptc_tuples = [tuple(pt...) for pt in ptc]
        prev = [(Y,X) for (X,Y) in ptc_tuples]
        #
        #rpt = ArchGDAL.createpoint.(prev)
        #reversed_polygon = ArchGDAL.createpolygon(coords)
        reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
        
        # Create a new polygon with the reversed coordinates
        #reversed_polygon = ArchGDAL.createpolygon.(reversed_crds)
        #reversed_polygon = ArchGDAL.createpolygon.(prev)

        return reversed_polygon
    end

    """
    plots timeseries like cmk, but with a legend
    ```
    dfpl(x::Union{Regex,String,DataFrame};leg::Bool=true,kw...)

    kw.. passed to Axis:
    cmk.dfpl(r"wolf"i;yscale=log10)
    using CairoMakie
    f = Figure()
    for (i, scale) in enumerate([identity, log10, log2, log, sqrt, Makie.logit])
        row, col = fldmod1(i, 2)
        Axis(f[row, col], xscale = scale, title = string(scale),
            xminorticksvisible = true, xminorgridvisible = true,
            xminorticks = IntervalsBetween(5))

        lines!(range(0.01, 0.99, length = 200), 1:200)
    end
    f
    ```
    
    https://docs.makie.org/stable/reference/blocks/legend/
    """
    function dfpl(x::Union{Regex,String,DataFrame};leg::Bool=true,kw...)
        if isa(x,DataFrame)
            df = (x)
        else
            df = waread(x)
        end

        dt = vec(Matrix(
            select(df,
                first(filter(x->occursin(r"date|year",x),
                names(df))))))
        #dt = df.date

        fig = Figure()
                # fonts=(;regular = "consolas")            
                # #size=(600, 400), 
                # )
        
        tempo = string.(dt)
        lentime = size(df, 1)
        #slice_dates = range(1, lentime, step=lentime ÷ 8)
        slice_dates = range(1, lentime, step=lentime ÷ 10)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        
        ax3 = Axis(fig[1, 1];kw...)

        cols = names(df)
        filter!(x -> !occursin(r"date|year", x), cols)

        colors = Makie.wong_colors()
        # Map unique values to colors
        value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # Map values in a to colors
        selected_colors = [value_to_color[value] for value in cols]

        for (i, col) in enumerate(cols)
            vals = Vector(df[:, col])
            lines!(ax3, 1:lentime, vals; color = selected_colors[i], linewidth = 0.85, label=col)
        end
        
        if leg
            lt = replace(tit, r"_|so" => " ")
            legend = Legend(fig, ax3, lt, 
                framevisible = true,
                nbanks = 2)

            fig[1, 2] = legend
            # axislegend(ax3, fig, 
            # cols, lt, 
            # position = :rb,
            # orientation = :horizontal)

        end

        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        
        return fig
    end

    tsdf = dfpl

    """
    density plot of dataframe monthly values
    kw passed to Axis
    msk -> only values > 0 
    ```
    mkemon(x::Union{String,Regex,DataFrame};
        col::Any="tot_average",
        tit::nothing,
        msk=false,
        mskval=0,
        offset=4, #offset denominator. increase on higher density
        kw...)
    ```
    """
    function mkemon(x::Union{String,Regex,DataFrame};col::Any="tot_average",
            ti::Union{String,LaTeXString} = "",
            msk=false,
            mskval=0,
            offset=4,
            kw...)
        if x isa String
            printstyled("reading $x\n",color=:light_red)
            df = dfr(x)
            dropmissing!(df)    
        elseif x isa Regex
            x = first(dfonly(x))
            printstyled("reading $x\n",color=:light_red)
            df = dfr(x)
            dropmissing!(df)
        else
            df = copy(x)
            dropmissing!(df)
        end

        # ti = try
        #     DataFrames.metadata(df) |> only |> last |> basename
        # catch
        #     @warn "No basename in metadata!"
        #     ti = raw""
        # end
        #check if "tot_average" in df
        if col isa Number && names(df)[col] == "date"
            if col + 1 >= ncol(df)
                col = col - 1
                @info "date column is at $col position, skipping to col-1..."
            elseif names(df)[col - 1] == "date"
                col = ncol(df)
                @info "date column is at $col-1 position, skipping to last col..."
            else
                col = col + 1
                @info "date column is at position $col, skipping to col + 1..."
            end
        end
        
        try 
            select!(df, Cols(:date,col))
        catch
            @warn "No :date & $col column in dataframe!"
            return
        end
        #only values > 0 
        if msk
            df = try 
                DataFrames.subset(df,col => ByRow(>(float(mskval))))
            catch
                cn=propertynames(df[!,Not(Cols(r"date|month"i))])|>first
                DataFrames.subset(df,cn => ByRow(>(float(mskval))))
            end
        end

        df[!, :month] = month.(df[!,:date]);
        grp = DataFrames.groupby(df, :month)

        months = ["Januar", "Februar", "März", "April",
            "Mai", "Juni", "Juli", "August", "September",
            "Oktober", "November", "Dezember"]
        
        #if isnothing(tit)
        if length(ti) == 0
            plot_title = ti*" Basin:"*string(col)
        else
            plot_title = ti
        end
        
        f = Figure()
        Axis(f[1, 1], title = plot_title,
            yticks = ((1:12) ./ 4,  reverse(months));kw...)
        for i in 12:-1:1
            values = select(grp[i], Not(Cols(r"date|month"i)))|>Matrix|>vec
            d = Makie.density!(values, offset = i / offset,
                        color = :x, colormap = :thermal, 
                        strokewidth = 1, strokecolor = :black)
                        #colorrange = (-5, 5),
            #Apply an absolute translation to the Scene
            translate!(d, 0, 0, -0.1i) 
            #translate!(d, 0, 0.5i, -0.1i) 
            #translate!(d, 0.5i, 0, -0.1i) 
            
        end
        return f
    end

    """
    ```
    mkrheat(x::Union{String,Regex,Raster};hide=false,msk=true,missval=0.0,umask=10^6,mskval=0.001,layer=1,crs=EPSG(25832),mappedcrs=EPSG(25832),turn=true,kw...) \n
    uses Rasters to read netcdf and makie to plot.
    Regex only applies to netcdf files. \n
    hide = hidedecorations \n
    turn = true flips the raster  \n
    agr = Raster(x)|>v->Rasters.aggregate(mean,v,2)
    i,m = extrema(agr|>skipmissing)
    cmk.mkrheat(agr;mskval=i,umask=m-0.0001)
    ```
    """
    function mkrheat(x::Union{String,Regex,Raster};
        hide=false,
        msk=true,
        layer=1,
        umask::Number=10^6,
        mskval::Number=0.001,
        missval::Number=0.0,
        crs=EPSG(25832),mappedcrs=EPSG(25832),turn=true,
        c=(:plasma,0.7),kw...)
        #(Reverse(:plasma),0.9) -9999.000000
        if isa(x,Regex)
            #fn = filter(z->endswith(z,"nc"),glob(x))[1]
            #ne = join([x,"nc","\$"],"+.")
            ds = try
                #fn = filter(z->endswith(z,"nc"),Grep.grep(x,readdir()))[1]
                fn = first(filter(z -> endswith(z, "nc") && occursin(x, z), readdir()))
                println("found $fn")
                if endswith(fn,"nc")
                    Raster(fn)
                else
                    Raster(fn;missingval=missval) #missingval,
                end
            catch
                @error "no file found!"
                return
            end
            
        elseif isa(x,Raster)
            ds = x #Rasters.rebuild(x;missingval=missval)
            fn = string.(name(ds))
        else
            fn = x
            ds = try
                if endswith(fn,"nc")
                    Raster(fn)
                else
                    Raster(fn;missingval=missval) #missingval,
                end
            catch
                @error "no file found!"
                return
            end
        end
        
        #A = ds[Ti=layer].data
        #A = ds[:, :, 1].data
        #mx = missingmask(ds)
            
        ds = ds[:, :, layer]
        if msk
            if umask==10^6
                umask = maximum(ds|>skipmissing)
            end
            #A = ds.data .> float(mskval)
            #ds = Rasters.rebuild(ds;missingval=0.0)
            #ds = replace_missing(ds,missingval)
            #bitmat = (dx .> msk) .& (dx .<= umask)
            zm = (ds .> float(mskval)) .& (ds .<= umask)
            #zm = ds .> float(mskval)
            #zm = boolmask(zm;missingval=0.0)
            zm = boolmask(zm;missingval=missval)
            
            A = Rasters.mask(ds; with=zm)
            min_val, max_val = extrema(skipmissing(A))
            @info join(["min: ", string(min_val), ", max: ", string(max_val)], " ")
            #
            if (typeof(A.metadata) == DimensionalData.Dimensions.LookupArrays.NoMetadata || isnothing(A.metadata))
                @error "no metadata in raster!"
            else
                try
                    @info A.metadata.val    #|>DataFame
                catch
                    @warn "no metadata in raster!"
                end
            end
            
            #A = Rasters.mask(ds; with=boolmask(ds;missingval=mskval))
            A = A.data
            #A = Rasters.rebuild(A;kw...)
            #A = ds.data
            #A = replace(A, missing => mskval)
            #A = A[A .> mskval]
        else
            A = ds.data
        end
        #see also Rasters.trim
        #extrema(A|>skipmissing)
        # A = replace(A, -9999.0 => missing)
        # A = replace(A, 0.0 => missing)
        #b = replace(A, [0 .. 70] => missing)
        #b = replace(A, float(-10..70) => missing)
        #b = map(x -> replace(x,float(-10..70) => missing),eachrow(A))     

        mrows = map(x -> !all(ismissing, x),eachrow(A))
        mcols = map(x -> !all(ismissing, x),eachcol(A))
        A = A[mrows, mcols]
        
        if turn
            A = transpose(A)
            A = reverse(A, dims=2) #very important!
        else
            #A = A #reverse(A, dims=1) #replace(A, missing => 0.0)
            #A = replace(A, missing => mskval)
            A = A
        end
        
        ti = replace(basename(fn),".nc"=>"","_"=>" ")
        #ti = replace(name(ds),".nc"=>"","_"=>" ")
        #https://docs.makie.org/stable/explanations/colors/
        fig = Figure(
            #color = Makie.wong_colors(),
            #color = :thermal, #see colormap below
            #size=(800, 600), 
            fontsize=22);
        axs = Axis(fig[1,1],
            title=ti, 
            aspect = DataAspect(), 
            xlabel="x", 
            ylabel="y")
        
        #p1 = heatmap!(A,colormap=(:turbo, 0.9))
        p1 = heatmap!(A,
        #colormap=(:thermal, 0.9)
        colormap = c,kw...
        )
        #p1 = heatmap!(A,colormap=Makie.wong_colors())
        
            contour!(axs, A; color=:black) #, levels=lscale
            Colorbar(fig[1, 2], p1, 
                width=20, ticksize=20, 
                tickalign=1)
            #xmax,ymax = size(A)
            #limits!(axs, 1, xmax, 1, ymax)
        # hideydecorations!(axs, grid=true, ticks=false)
        if hide 
            hidedecorations!(axs, grid=true, ticks=true)
        end
        return fig
    end

    """
    reads ascii grid and replaces -9999 \n returns Raster
    """
    function agread(x::String)
        dataset = ArchGDAL.read(x)
        # Get the first band
        band = ArchGDAL.getband(dataset, 1)
        # A = permutedims(band, (2, 1))
        # A = A.a
        # A = reverse(A, dims=1)
        A = reverse(band, dims=2)
        B = replace(A, -9999.0 => missing)
        xm = .!all(ismissing, B, dims=2)
        ym = .!all(ismissing, B, dims=1)
        A = B[xm[:], ym[:]]
        k = Raster(A,(X,Y))
        println(extrema(k))
        return k
    end  

    """
    streamplot
    streamplot(f::function, xinterval, yinterval; color = norm, kwargs...)
    """                     
    function mkestream(r::Raster;wasim=true,missval=0)
        if wasim
            z = replace_missing(r[t=1], missval)
            reverse!(z,dims=1)
            # z = z'       
            # replace!(z, missval=>missing)
        else
            z = r.data
            # reverse!(z,dims=1)
            # reverse!(z,dims=2)
            # replace!(z, missval=>missing)
        end
        #ex = extrema(z|>skipmissing)
        
        # Define a function to calculate the gradient
        function gradient2d(z)
            dy, dx = size(z)
            gx = [z[i, j+1] - z[i, j] for i in 1:dy, j in 1:dx-1]
            gy = [z[i+1, j] - z[i, j] for i in 1:dy-1, j in 1:dx]
            return gx, gy
        end

        gx, gy = gradient2d(z)
        f = (x, y) -> Point2(gx[round(Int, y), round(Int, x)], gy[round(Int, y), round(Int, x)])
        xinterval = 1:size(z, 2)
        yinterval = 1:size(z, 1)       
        p1 = CairoMakie.streamplot(f, xinterval, yinterval)
        return p1
    end


    """
    new mkrheat
    ```
    mkh(x::Union{String,Regex,Raster};
        hide=false,
        msk=true,
        proj=false,
        layer=1,
        mskval::Number=0.001,
        umask::Number=10^6,
        missval::Number=0.0,
        turn=true,
        psize = (900, 700),
        c=(:plasma,0.7),kw...)
    ```
    """
    function mkh(x::Union{String,Regex,Raster};
        hide=false,
        msk=true,
        proj=false,
        layer=1,
        mskval::Number=0.001,
        umask::Number=10^6,
        missval::Number=0.0,
        assign_missing::Bool=true,
        psize = (900, 700),
        turn=true,
        c=(:plasma,0.7),
        kw...)
        #(Reverse(:plasma),0.9) -9999.000000
        if isa(x,Regex)
            ds = try
                fn = first(filter(z -> endswith(z, "nc") && occursin(x, z), readdir()))
                println("found $fn")
                if endswith(fn,"nc")
                    Raster(fn;
                    crs=EPSG(25832),
                    mappedcrs=EPSG(25832))
                else
                    Raster(fn;missingval=missval)
                end
            catch
                @error "no file found!"
                return
            end
        elseif isa(x,Raster) && assign_missing
            ds = Rasters.rebuild(x;missingval=missval)
            fn = string.(name(ds))
        elseif isa(x,Raster)
            ds = x
            fn = string.(name(ds))
        else
            fn = x
            ds = try
                if endswith(fn,"nc")
                    Raster(fn;
                    crs=EPSG(25832),
                    mappedcrs=EPSG(25832))
                else
                    Raster(fn;missingval=missval)
                end
            catch
                @error "no file found!"
                return
            end
        end       
            
        ds = ds[:, :, layer]
        if proj
            #For converting between projections that are rotated, skewed or warped in any way, use resample.
            #ds = Rasters.reproject(ds, EPSG(4326))
            #If projections can be converted for each axis independently, 
            #it may be faster and more accurate to use reproject.
            ds = Rasters.resample(ds;
                crs=ProjString("+proj=longlat +datum=WGS84"),
                method=:bilinear)
        end


        if msk
            if umask==10^6
                umask = maximum(ds|>skipmissing)
            end
            zm = (ds .> float(mskval)) .& (ds .<= umask)
            zm = boolmask(zm;missingval=missval)
            
            A = Rasters.mask(ds; with=zm)
            min_val, max_val = extrema(skipmissing(A))
            @info join(["min: ", string(min_val), ", max: ", string(max_val)], " ")
            if (typeof(A.metadata) == DimensionalData.Dimensions.LookupArrays.NoMetadata || isnothing(A.metadata))
                @error "no metadata in raster!"
            else
                try
                    @info A.metadata.val    #|>DataFame
                catch
                    @warn "no metadata in raster!"
                end
            end
            
            A = Rasters.trim(A)
            #extrema(A|>skipmissing)
            #A = Rasters.rebuild(A;kw...)
        else
            A = Rasters.trim(ds)
        end
        if turn
            A = transpose(A)
            A = reverse(A, dims=2) #very important!
        end
        ti = replace(basename(fn),".nc"=>"","_"=>" ")
        fig = GeoMakie.Figure(
            #color = :thermal, #see colormap below
            size = psize, 
            fontsize=20);
        # axs = Makie.Axis(fig[1,1],
        #     title=ti, 
        #     aspect = DataAspect(), 
        #     xlabel="x", 
        #     ylabel="y")
        
        # Set custom x-axis tick positions (e.g., every 2 units)
        #Axis(xticks = 2:2:10)
        xl = lookup(A,:X)|>collect|>extrema
        yl = lookup(A,:Y)|>collect|>extrema
        xt = round.(range(xl[1], stop=xl[2], length=4),digits=2)
        #x_values = range(525600, stop=621000, length=5)
        
        axs =  Makie.Axis(fig[1,1],
            title   =   ti,
            #xticks = xl[1]:250*10^5:xl[2],
            xticks = xt,
            xticklabelrotation = 45,
            xminorticksvisible = true,
            xminorgridvisible = true,
            yminorticksvisible = true, 
            yminorgridvisible = true,
            xminorticks = IntervalsBetween(5),
            yminorticks = IntervalsBetween(5))
        p1 = GeoMakie.heatmap!(axs,A;
                colormap = c,
                kw...)
            #CairoMakie.
            Makie.contour!(axs, A; color=:black,
                linewidth = 0.5) 
            #, levels=lscale
            Colorbar(fig[1, 2], p1, 
                width=20, ticksize=20, 
                tickalign=1)


            #add 100 to each
            #xl = (xl[1]-100, xl[2]+100)
            #yl = (yl[1]-100, yl[2]+100)
            if GeoInterface.crs(ds)  != EPSG(25832)
                inc = 0.01
            else
                inc = 500
            end
            GeoMakie.limits!(
                xl[1]-inc, 
                xl[2]+inc,
                yl[1]-inc, 
                yl[2]+inc)
                    
        if hide 
            hidedecorations!(axs, grid=true, ticks=true)
        end
        return fig
        #GeoMakie.save("x.png",fig, px_per_unit = 10/1)
    end

    function create_raster(fn; missval::AbstractFloat = -9999.0)
        if endswith(fn,"nc")
            tmp=Raster(fn;crs=EPSG(25832),mappedcrs=EPSG(25832))
            return tmp[t=1]
        else
            return Raster(fn;crs=EPSG(25832),mappedcrs=EPSG(25832),missingval=missval)
        end
    end
    
    """
    ```
    mkfzs(x::Union{String,Regex,Raster};
        ti::String="Fliesszeitsummen [h]",
        tolonlat=false,hide=false,
        missval::AbstractFloat = -9999.0,
        c=(:greys,0.7), kw... )
    ```
    kw... passed to heatmap!
    """
    function mkfzs(x::Union{String,Regex,Raster};
        ti::String="Fliesszeitsummen [h]",
        tolonlat=false,hide=false,
        missval::AbstractFloat = -9999.0,
        c=(:greys,0.7), kw... )
    
        if isa(x,Regex)
            #fn = first(filter(z -> !endswith(z, "nc") && occursin(x, z), readdir()))
            fn = first(filter(z -> occursin(x, z), readdir()))
            if isfile(fn)
                ds = create_raster(fn; missval=missval)
            else
                @error "no file found!"
                return
            end
        elseif isa(x,Raster) && hasdim(x, :t)
            ds = x[t=1]
            fn = string.(name(ds))
        elseif isa(x,Raster)
            ds = x              # x[Band(1)]
            fn = string.(name(ds))
        else
            fn = x
            if isfile(fn)
                ds = create_raster(fn; missval=missval)
            else
                @error "no file found!"
                return
            end
        end
        if tolonlat
            ds = resample(ds;crs=ProjString("+proj=longlat +datum=WGS84"))
            lonlat = true
        end
        z = replace_missing(ds,0)
        z.data .= z.data ./ 3600
        xl = lookup(z,:X)|>collect|>extrema
        yl = lookup(z,:Y)|>collect|>extrema
        xt = round.(range(xl[1], stop=xl[2], length=4),digits=2)

        fig = CairoMakie.Figure(fontsize=15);
        axs = CairoMakie.Axis(fig[1,1],
            title = ti,
            xticks = xt,
            xticklabelrotation = 45,
            xminorticksvisible = true,
            xminorgridvisible = true,
            yminorticksvisible = true, 
            yminorgridvisible = true,
            xminorticks = IntervalsBetween(5),
            yminorticks = IntervalsBetween(5))
        p1 = CairoMakie.heatmap!(axs,z;
                colormap = c,
                kw...)
            # CairoMakie.contour!(axs, z; color=:black,
            #     linewidth = 0.5) 
            CairoMakie.Colorbar(fig[1, 2], p1, 
                width=20, ticksize=20, 
                tickalign=1)
            

            if lonlat
                inc = 0.01
            else
                inc = 500
            end
            CairoMakie.limits!(
                xl[1]-inc, 
                xl[2]+inc,
                yl[1]-inc, 
                yl[2]+inc)
                    
        if hide 
            hidedecorations!(axs, grid=true, ticks=true)
        end
        return fig
    end

    """
    wasim plot
    Union{String,Raster}, fn::Union{String,GeoMakie.GeoJSON.FeatureCollection
    lyr::Int = 1,
    lims::Tuple{Tuple{Int,Int},Tuple{Int,Int}}=((8, 13), (48, 52)),c=Reverse(:plasma))
    """
    function gmwas(rmo::Union{String,Raster}, fn::Union{String,GeoMakie.GeoJSON.FeatureCollection};
            lyr::Int = 1,
            #lims::Tuple{Tuple{Int,Int},Tuple{Int,Int}}=((9, 12), (49, 51)),
            #lims::Tuple{Tuple{Int,Int},Tuple{Int,Int}}=((9.5, 12), (48.5, 50.7)),
            c=Reverse(:plasma))
        if fn isa GeoMakie.GeoJSON.FeatureCollection
            pol = fn
        else
            pol = GeoMakie.GeoJSON.read(read(fn, String))
        end
        if rmo isa Raster
            ra = rmo[t=lyr]
            ra = replace_missing(ra,0)
            #ra = rst.project(ra;dst="+proj=longlat +datum=WGS84")
            if crs(ra) == EPSG(25832)
                ra = resample(ra;crs=ProjString("+proj=longlat +datum=WGS84"))
            end
        else
            r = Raster(rmo;crs=EPSG(25832),mappedcrs=EPSG(25832))
            ra = r[t=lyr]
            ra = replace_missing(ra,0)
            #ra = rst.project(ra;dst="+proj=longlat +datum=WGS84")
            ra = resample(ra;crs=ProjString("+proj=longlat +datum=WGS84"))
        end
        
        fig = Figure()
        #ax = GeoMakie.GeoAxis(fig[1,1]; limits=lims)
        #ax = Axis(fig[1,1]; limits=lims)
        ax = Axis(fig[1,1];)
        sf = GeoMakie.heatmap!(ax,ra;
            #colormap = Reverse( :plasma, 0.5 ),     
            #colormap = ["grey", "red", "orange", "yellow", "blue"],
            colormap = c)
            #shading = NoShading)
        Colorbar(fig[1,2], sf; 
            label = "",
            height = Relative(0.65))
        #n = length(pol)
        #hm = GeoMakie.poly!(ax, pol; 
        hm = poly!(ax, pol; 
            color = :transparent,    
            strokecolor = :black, 
            strokewidth = 1.5,
            )
        #GeoMakie.translate!(hm, 0, 0, 100) # move above surface plot
        translate!(hm, 0, 0, 100) # move above surface plot
        return fig
    end

    """
    ``` 
    function mkflow(TS, q=[0.9, 0.1]; text::String="", ylab::LaTeXString = <liter per second*km²>,leg=true,kw...)
    usage: cmk.mkflow(df;text="Flow Series \n") #Time of df is in Legend Title
    ```
    """
    function mkflow(TS, q=[0.9, 0.1]; text::String="", ylab::LaTeXString = L"\left[\frac{l}{s\cdot km^2}\right]",leg=true,kw...)
        TS = copy(TS)
        if !hasproperty(TS, :date)
            @error "no date column in dataframe!"
            return
        end
        if first(propertynames(TS)) == :date
            #date to last position
            TS = select(TS, Not(:date), :)
        end
        timestring = string.(extrema(TS.date))
        text = text*" "*timestring[1]*" - "*timestring[2]

        mn = ["Okt", "Nov", "Dez", "Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep"]
        #TS.date = Date.(TS.date) # Ensure Date column is of type Date
        #TS.HydroYear = year.(TS.date)
        TS.Month = month.(TS.date)
        TS.DOY = Dates.dayofyear.(TS.date)
        TS.hdoy = mod.(TS.DOY .- 275, 365) .+ 1
        #TS.hmon = mod.(TS.Month .- 10, 12) .+ 1
        
        # Map the months to their names
        TS.mn = map(x -> mn[Int(x)], TS.Month)

        rename!(TS, 1 => :Flow)
        doy = (TS.DOY)
        Qdoy = fill(NaN, maximum(map(Int, doy)), 7)      
        Qdoy[:, 1] .= [maximum(TS.Flow[doy .== d]) for d in levels(doy)]
        Qdoy[:, 2] .= [minimum(TS.Flow[doy .== d]) for d in levels(doy)]
        Qdoy[:, 3] .= [mean(TS.Flow[doy .== d]) for d in levels(doy)]
        Qdoy[:, 4] .= [quantile(TS.Flow[doy .== d], q[1]) for d in levels(doy)]
        Qdoy[:, 5] .= [quantile(TS.Flow[doy .== d], q[2]) for d in levels(doy)]
        Qdoy[:, 6] .= [median(TS.Flow[doy .== d]) for d in levels(doy)]
        Qdoy[:, 7] .= [first(TS.Month[doy .== d]) for d in levels(doy)]

        nms = ["MaxQ", "MinQ", "MeanQ", "Q90", "Q10", "Median", "Month"]
        df = DataFrame(Qdoy, nms)
        df[!, "SortOrder"] = ifelse.(df.Month .>= 10, df.Month .- 12, df.Month)
        sort!(df, :SortOrder)

        # Create the plot

        fig = Figure(size = (800, 600))
        ax3 = Axis(fig[1, 1], 
            #title = text, #see legend
            #palette = (patchcolor = [:red, :green, :blue, :grey, :orange],),
            #xlabel = "DOY", 
            ylabel = ylab) 
        cols = names(df)
        filter!(x -> !occursin(r"date|year|sortorder|month"i, x), cols)
        #colors = Makie.wong_colors()
        #colormap = Reverse(:bluesreds),
        colors = [:grey, :red, :green, :blue, :grey, :orange]
        value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        selected_colors = [value_to_color[value] for value in cols]
        lentime = size(df, 1)
        #with_theme(
            # Theme(
            # palette = (color = [:grey, :blue], linestyle = [:dash, :dot]),
            # Lines = (cycle = Cycle([:color, :linestyle], covary = true),)
            # ) #) do
        for (i, col) in enumerate(cols)
            vals = Vector(df[:, col])
            lines!(ax3, 1:lentime, vals; 
                #linestyle = :dash,    
                color = selected_colors[i], 
                linewidth = 0.95, 
                label=col)
        end
        if leg
            legend = Legend(fig, ax3, text, 
                framevisible = true,
                nbanks = 2)

            fig[1, 2] = legend
        end
        
        #tempo = string.(df.Month)
        tempo = TS.mn #month strings mapped on dates.
        slice_dates = range(1, lentime, step=lentime ÷ 10)
        ax3.xticks = (slice_dates, tempo[slice_dates])
        
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        
        return fig
    end

    """
    plots timeseries, see dfpl
    """
    function fts(df::DataFrame)
        
        dt = vec(Matrix(
            select(df,
                first(filter(x->occursin(r"date|year|time"i,x),
                names(df))))))
        fig = Figure( #resolution=(600, 400), 
            #fonts=(;regular = "cmr10")
            #fonts=(;regular = "consolas")
            )
        
        tempo = string.(dt)
        lentime = size(df,1)
        slice_dates = range(1, lentime, step=lentime ÷ 8)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        #xline =  first(DelimitedFiles.readdlm(x,'_',String)[2,:])
        ax3 = Axis(
            fig[1, 1], 
            title = replace(tit,r"_|so"=>" ")
                #xlabel="Date", 
                #yscale = :log10,
                #ylabel=first(names(df))
                )
        #cols = propertynames(df)  
        cols = names(df)  
        filter!(x->!occursin(r"date|year",x),cols)
        # colors = Makie.wong_colors()
        # # Map unique values to colors
        # value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # # Map values in a to colors
        # selected_colors = [value_to_color[value] for value in cols]
        # colors = Makie.wong_colors()
        #v = 1:length(cols)
        #for v, col in enumerate(cols)
        for col in cols
            vals = Vector(df[:, col])
            lines!(ax3, 1:lentime, vals; 
            color=:black,         
            # color=colors[v], 
             linewidth=0.75)
        end
        
        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        return fig
    end

    
    """
    density plot of dataframe monthly values
    kw passed to Figure
    msk -> only values > 0 
    ```
    mkemon2(x::Union{String,Regex,DataFrame},y::Union{String,Regex,DataFrame};
        col::Any="tot_average",
        col2::Any="tot_average",
            tit::String = "",
            msk=false,kw...)
    ```
    """
    function mkemon2(x::Union{String,Regex,DataFrame},y::Union{String,Regex,DataFrame};
        col::Any="tot_average",
        col2::Any="tot_average",
            tit::String = "",
            msk=false,kw...)
        if x isa String
            printstyled("reading $x\n",color=:light_red)
            df1 = dfr(x)
            df2 = dfr(y)
            dropmissing!(df1)
            dropmissing!(df2)
        elseif x isa Regex
            x = first(dfonly(x))
            x = first(dfonly(y))
            printstyled("reading $x\n",color=:light_red)
            printstyled("reading $y\n",color=:light_red)
            df1 = dfr(x)
            df2 = dfr(y)
            dropmissing!(df1)
            dropmissing!(df2)
        else
            df1 = copy(x)
            df2 = copy(y)
            dropmissing!(df1)
            dropmissing!(df2)
        end

        ti = try
            DataFrames.metadata(df1) |> only |> last |> basename
        catch
            @warn "No basename in df1!"
            ti = raw""
        end
        #check if "tot_average" in df1
        if col isa Number && names(df1)[col] == "date"
            if col + 1 >= ncol(df1)
                col = col - 1
                @info "date column is at $col position, skipping to col-1..."
            elseif names(df1)[col - 1] == "date"
                col = ncol(df1)
                @info "date column is at $col-1 position, skipping to last col..."
            else
                col = col + 1
                @info "date column is at position $col, skipping to col + 1..."
            end
        end
        #check if "tot_average" in df2
        if col2 isa Number && names(df2)[col2] == "date"
            if col2 + 1 >= ncol(df2)
                col2 = col2 - 1
                @info "date column is at $col2 position, skipping to col-1..."
            elseif names(df2)[col2 - 1] == "date"
                col2 = ncol(df2)
                @info "date column is at $col2-1 position, skipping to last col..."
            else
                col2 = col2 + 1
                @info "date column is at position $col2, skipping to col + 1..."
            end
        end
        
        
        try 
            select!(df1, Cols(:date,col))
            select!(df2, Cols(:date,col2))
        catch
            @warn "subset dataframe failed!"
            return
        end

        #only values > 0 
        if msk
            df1 = try 
                DataFrames.subset(df1,col => ByRow(>(0)))
            catch
                cn=propertynames(df1[!,Not(Cols(r"date|month"i))])|>first
                DataFrames.subset(df1,cn => ByRow(>(0)))
            end
            df2 = try 
                DataFrames.subset(df2,col => ByRow(>(0)))
            catch
                cn2=propertynames(df2[!,Not(Cols(r"date|month"i))])|>first
                DataFrames.subset(df2,cn2 => ByRow(>(0)))
            end
        end


        df1[!, :month] = month.(df1[!,:date]);
        grp1 = DataFrames.groupby(df1, :month)
        df2[!, :month] = month.(df2[!,:date]);
        grp2 = DataFrames.groupby(df2, :month)

        months = ["Januar", "Februar", "März", "April",
            "Mai", "Juni", "Juli", "August", "September",
            "Oktober", "November", "Dezember"]
        
        #if isnothing(tit)
        if length(tit) == 0
            plot_title = ti*" Basin:"*string(col)
        else
            plot_title = tit
        end
        
        f = Figure(;kw...)
        Axis(f[1, 1], title = plot_title,
            yticks = ((1:12) ./ 4,  reverse(months)))
        for i in 12:-1:1
            values = select(grp1[i], Not(Cols(r"date|month"i)))|>Matrix|>vec
            d = density!(values, offset = i / 4,
                        color = :x, colormap = :thermal, 
                        strokewidth = 1, strokecolor = :black)
                        values = select(grp2[i], Not(Cols(r"date|month"i)))|>Matrix|>vec
                        
                        density!(values, offset = i / 4,
                                    color = :x, 
                                    colormap = [:blue, :yellow, :red], 
                                    strokewidth = 1, 
                                    strokecolor = :black)
            #Apply an absolute translation to the Scene
            translate!(d, 0, 0, -0.1i) 
        end
        # for i in 12:-1:1
        #     values = select(grp2[i], Not(Cols(r"date|month"i)))|>Matrix|>vec
        #     density!(values, offset = i / 4,
        #                 color = :x, colormap = :thermal, 
        #                 #colorrange = (-5, 5),
        #                 strokewidth = 1, strokecolor = :black)
        #     #Apply an absolute translation to the Scene
        #     translate!(d, 0, 0, -0.1i) 
        # end
        return f
    end

    """
    density plot of dataframe on monthly values on annual basis
    kw passed to Axis
    ```
    density_plot(df::Union{DataFrame,Regex,String}; date_col::Symbol=:date, value_col::Symbol=:value,agg=mean, kwargs...)
    ```
    """
    function density_plot(df::Union{DataFrame,Regex,String}; date_col::Symbol=:date, 
            value_col::Union{Symbol,String,Int}="",agg=mean, kwargs...)
        if isa(df, Regex) || isa(df, String)
            df = dfr(df)
        end
        if value_col isa String
            value_col = Symbol(value_col)
        end
        if value_col isa Int
            value_col = propertynames(df)[value_col]
            @info "value_col is $value_col"
        end
        df = copy(df)
        if !(value_col in propertynames(df))
            value_col = first(propertynames(df[!,Not(:date)]))
            #@info "value_col not found in DataFrame! 
            @info "Using first column as value_col which is $value_col"
        end
        # Extract month and day
        df.month = month.(df[!,date_col]);
        df.day = day.(df[!,date_col]);
        # Aggregate data by month and day
        agg_df = DataFrames.combine(DataFrames.groupby(df, [:month, :day]), value_col => agg => :mean_value)
        density_data = [agg_df.mean_value[agg_df.month .== m] for m in 1:12]
        labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        reverse!(labels) # Reverse the order of the labels

        grs = ["#1C1C1C", "#2E2E2E", "#404040", "#525252", "#646464", 
        "#767676", "#888888", "#9A9A9A", "#ACACAC", "#BEBEBE", 
        "#D0D0D0", "#E2E2E2"]
        thermal_colors = [
            "#ADD8E6",  # January (Sky Blue - Cool)
            "#87CEFA",  # February (Light Sky Blue - Mild)
            "#00BFFF",  # March (Deep Sky Blue - Mild Warm)
            "#00FF7F",  # April (Spring Green - Warm)
            "#FFD700",  # May (Gold - Warm)
            "#FFA500",  # June (Orange - Hot)
            "#FF4500",  # July (Orange Red - Very Hot)
            "#FF6347",  # August (Tomato - Hot)
            "#FF7F50",  # September (Coral - Mild Hot)
            "#FFB6C1",  # October (Light Pink - Cooling)
            "#D3D3D3",   # November (Light Gray - Cool)
            "#B0E0E6"  # December (Light Blue - Cool)
        ]
        
        # Plot density
        fig = Figure()
        ax = Axis(fig[1, 1], ylabel = string(value_col);kwargs...)
        for (i, data) in enumerate(density_data)
            d = density!(ax, data, 
                label = labels[i],
                strokewidth = .66, 
                strokecolor = :black,
                color = thermal_colors[i])
            #Apply an absolute translation to the Scene
            #translate!(d, 0, 0, -0.2i) 
            translate!(d, 0, 0.1i, -0.1i) 
            #translate!(d, 0.5i, 0, -0.1i) 
            
        end
        #Label(fig[1, 1], string(value_col), halign = :right, valign = :top,
        #        justification = :center)
        #fig[1, 2] = Legend(fig, ax, string(value_col), framevisible = false)
        axislegend(; merge = true, position = :rt)
        #, orientation = :horizontal)
        # Legend(fig[1,1], ax, string(value_col), 
        #     framevisible = false,
        #     height = Relative(1 / 4),
        #     width = Relative(1 / 4),
        #     tellheight = false,
        #     tellwidth = false,
        #     #margin = (30, 10, 10, 10),
        #     halign = :right, valign = :top, 
        #     orientation = :horizontal)

        return fig
    end

    """
    Contour over heatmap
    ```
    mkc(r::Raster,wasim=false,missval=0)
    ```
    """                     
    function mkc(r::Raster,wasim=false,missval=0)
        if wasim
            z = r.data[:,:,1]
            reverse!(z,dims=1)
            z = z'
            replace!(z, missval=>missing)
        else
            z = r.data
            reverse!(z,dims=1)
            reverse!(z,dims=2)
            replace!(z, missval=>missing)
        end
        # if missingval(z) != missval
        #     z = replace_missing!(z,0)
        # end
        # ex = extrema(z|>skipmissing)
        # lscale = ex[1]:10:ex[2]   # errors by NaN32
        fig = Figure(
            #size=(800, 600), 
            fontsize=22);
        axs = Axis(fig[1,1], aspect=1, xlabel="x", ylabel="y")
        p1 = heatmap!(axs, z, colormap=(:plasma, 0.87))
        contour!(axs, z; color=:black) # levels=lscale)
        Colorbar(fig[1, 2], p1, width=20, ticksize=20, tickalign=1)
        xmax = findmax(size(z))[1] #gets the index of the max value
        limits!(axs, 1, xmax, 1, xmax)
        return fig
    end
    
    """
    ``` performance_metrics(df::DataFrame) ```
    """
    function performance_metrics(df::DataFrame)
        # Create a figure
        fig = Figure()
        
        # Create an axis
        ax = Axis(fig[1, 1], 
            title = "Performance Metrics Over Years", 
            xlabel = "Year", 
            ylabel = "Metric Value"
        )
        
        # Plot the lines
        lines!(ax, df.year, df.nse, color = :blue, linewidth = 2, label = "NSE")
        lines!(ax, df.year, df.kge, color = :red, linewidth = 2, label = "KGE")
        lines!(ax, df.year, df.ve, color = :green, linewidth = 2, label = "VE")
        
        # Add a legend
        axislegend(ax, position = :rt)
        
        return fig
    end

    """
    ``` cumplot(df::DataFrame; leg::Bool=true, scale = identity, oneyear::Bool=false) ```
    plots cumulated values of each column of a dataframe
    """
    function cumplot(x::DataFrame; leg::Bool=true, scale = identity, oneyear::Bool=false)
        if isa(x,DataFrame)
            df = (x)
        else
            df = waread(x)
        end        
        #enumerate([identity, log10, log2, log, sqrt, Makie.logit])
        dropmissing!(df)
        if oneyear
            lastfullyear = year(last(df.date))
            df = filter(row -> year(row.date) == lastfullyear, df)
            #mean_values = DataFrames.combine(df, names(df, Not(:date)) .=> mean)
            #println(mean_values)
        end
        
        fig = Figure()
        #make shure to start x-axias label at first date value
        #a,b = extrema(df[!, :date])
        #a,b = (1, size(df, 1))
        ax = Axis(fig[1, 1], yscale = scale,
            #limits = ((a, b), nothing),
            #xlabel = "Date", 
            ylabel = "", xlabel = "")
            #ylabel = "Cumulated Rainfall")
        #labs = []
        for col in filter(c -> c != :date, propertynames(df)) # exclude the "date" column
            if eltype(df[!, col]) <: Number # check if the column is numeric
                lines!(ax, df[!, "date"], 
                    cumsum(df[!, col]), 
                    label = string.(col))
                #push!(labs, col)
            end
        end
        ax.xticks = (4:size(df, 1)+4, df[!, "date"]) # show dates on x-axis
        ax.xticklabelrotation = π / 4 # rotate x-axis labels for better readability
        ax.xticklabelalign = (:right, :center)
        #hidedecorations!(ax)
        if leg
            tit = try
                values(DataFrames.metadata(df)) |> only |> basename
                catch
                    @warn "no metadata in df"
                    raw""
                end
            lt = replace(tit, r"_|so" => " ")
            legend = Legend(fig, ax, lt, 
                framevisible = true, nbanks = 2)
            fig[1, 2] = legend
        end
        fig
    end


end   ##endof end of module

# #https://docs.makie.org/stable/explanations/fonts/
# Makie.to_font("Computer Modern")
# Makie.theme(:fonts)
# #https://docs.makie.org/stable/explanations/latex/
# using CairoMakie
# with_theme(theme_latexfonts()) do
#     fig = Figure()
#     Label(fig[1, 1], "A standard Label", tellwidth = false)
#     Label(fig[2, 1], L"A LaTeXString with a small formula $x^2$", tellwidth = false)
#     Axis(fig[3, 1], title = "An axis with matching font for the tick labels")
#     fig
# end


printstyled("
    cairomakie.jl as module cmk loaded \n
    using .cmk to access functions from Main\n
            ",
    color=:green)
