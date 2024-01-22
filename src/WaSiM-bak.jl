#julia --startup-file=no -q --color=yes --project=. #-e 'using Pkg; Pkg.instantiate(); Pkg.API.precompile()'
#include("src/wa.jl")
# cd(raw"C:\Users\chs72fw\.julia\dev\WaSiM")
# pt="/mnt/c/Users/chs72fw/.julia/dev/WaSiM";cd(pt)
# using Pkg;Pkg.activate(".") 

###this is wa.jl

#using WaSiM

#module Startup
#test WaSiM

# using PrecompileTools    # this is a small dependency
    # @compile_workload begin
        # using DataFrames, DataFramesMeta, CSV, Statistics, Dates, StatsPlots, Distributions
        # using Plots.PlotMeasures
        # using DelimitedFiles, Grep, Printf
        # using Rasters, ArchGDAL
        # import NCDatasets
        # using PyCall
        # using KernelDensity
        # using SHA
    # end

module WaSiM

    using DataFrames, DataFramesMeta, CSV, Statistics, Dates, StatsPlots, Distributions
    using Plots.PlotMeasures
    #using PrettyTables #err

    using DelimitedFiles, Grep, Printf
    using Rasters, ArchGDAL
    import NCDatasets
    using PyCall
    # using PyPlot            #for pyplot_df, moved to pyjl
    using KernelDensity
    using SHA

    default(show = true)
    using PrecompileTools    # this is a small dependency
#    @compile_workload begin

    if Sys.isapple()
        platform = "osx"
        const homejl = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
        const mybash = "/Users/apfel/.bash_aliases"
        src_path = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
    elseif Sys.iswindows()
        platform = "windows"
        src_path = "C:\\Users\\Public\\Documents\\Python_Scripts\\julia"
        macro wasim() pt="C:\\Users\\chs72fw\\.julia\\dev\\WaSiM\\src\\wa.jl";include(pt);end
    else
        platform = "unix"
        winpt = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
        pcld = "~/pCloud Drive/Stuff/Python_Scripts/julia"
        src_path = isdir(winpt) ? winpt : pcld
        println("sourcepath is $src_path")
        if isdir(winpt)
            macro wasim() pt="/mnt/c/Users/chs72fw/.julia/dev/WaSiM/src/wa.jl";include(pt);end
        end
    end   

        function dfplotjs(df::DataFrame;logy::Bool,fact::Float64)
        nrows=size(df)[2]-1 
        #length(names(df))-1
        o = DataFrames.metadata(df)|>collect
        ti = basename(o[1][2])
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            #rows=2, cols=2
                );
            for i in 1:nrows;
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            fact = isnothing(fact) ? 1 : fact; #nice
            logy = isnothing(logy)==true ? logy==false : logy==true;
            if logy == true
                PlotlyJS.relayout!(fig,yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            #elseif isnothing(log) 
            else
                PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            end
            display(fig)
        end

        function dfplotjs(df::AbstractString;logy::Bool,fact::Float64)
            df=readmeteo(df)
            nrows=size(df)[2]-1
            #length(names(df))-1
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            fig = PlotlyJS.make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in 1:nrows;
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            fact = isnothing(fact) ? 1 : fact; #nice
            logy = isnothing(logy)==true ? logy==false : logy==true;
            if logy == true
                PlotlyJS.relayout!(fig,yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            else
                PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            end
            display(fig)
        end

        function dfplotjs(filepath::AbstractString)
            dfplotjs(filepath;logy=false,fact=1.0)
        end

        function dflogjs(filepath::AbstractString)
            dfplotjs(filepath;logy=true,fact=1.0)
        end

        function recursive_glob_prfx(rootdir=".", prefix="")
            results = []
            for (looproot, dirs, filenames) in walkdir(rootdir)
                for filename in filenames
                    if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
            return results
        end

        function qgk(rootdir=".", prefix="qgk")
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

        function penman_monteith(ETo, G, T, Td, u2, es, ea, Ra)
            """
            Calculates the potential evapotranspiration (PET) using the Penman-Monteith equation.

            Parameters
            ----------
            ETo : Float64
                Reference evapotranspiration (mm/day).
            G : Float64
                Soil heat flux density (mm/day).
            T : Float64
                Air temperature (°C).
            Td : Float64
                Dew point temperature (°C).
            u2 : Float64
                Wind speed at 2 m height (m/s).
            es : Float64
                Saturation vapor pressure (kPa).
            ea : Float64
                Actual vapor pressure (kPa).
            Ra : Float64
                Aerodynamic resistance (s/m).

            Returns
            -------
            PET : Float64
                Potential evapotranspiration (mm/day).
            """
            # Constants
            R = 8.314 # J/mol/K
            cp = 1.013e-3 # kJ/g/K

            # Latent heat of vaporization (MJ/kg)
            Lambda = 2.501 - 0.002361 * T

            # Psychrometric constant (kPa/°C)
            gamma = cp * P / (0.622 * Lambda)

            # Slope of the saturation vapor pressure curve (kPa/°C)
            delta = 4098 * es / (T + 237.3) ^ 2

            # Net radiation (MJ/m2/day)
            Rn = (1 - 0.23) * ETo

            # Air density (kg/m3)
            rho = P * 1000 / (R * (T + 273.15))

            # Specific heat of air (kJ/kg/K)
            cpa = 1.013 * rho ^ -0.0065 * 1000

            # Delta term (MJ/m2/day/°C)
            delta_term = (delta / (delta + gamma)) * (Rn - G)

            # Psi term (MJ/m2/day)
            psi_term = (gamma / (delta + gamma)) * rho * cp * (es - ea) / Ra * u2

            # Potential evapotranspiration (mm/day)
            PET = (delta_term + psi_term) / Lambda

            return PET
        end

        function jdd()
            cwd = pwd()
            dirs = readdir(".")
            for dir in dirs
                if isdir(dir)
                    size = 0
                    for (root, dirs, files) in walkdir(dir)
                        for file in files
                            size += stat(joinpath(root, file)).size
                        end
                    end
                @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
                end
            end
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
            # fn = fn/1024^2
            # printstyled("$fn MB",color=:blue)
        end 


        # printstyled("folders smaller 1MB will be omitted...\n",color=:red)
        # print_sorted_sizes(pwd())

        ##w endswith
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

        function ddense(path::String,skip::Int,start::Int,stop::Int)
            ms=["-999","-9999","lin","log","LIN","LOG"]
            df = CSV.read(path,DataFrame,skipto=skip,
            missingstring=ms,delim="\t",comment="-",silencewarnings=false,
            ntasks=4,downcast=true,normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
            df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
            df=df[:,Not(1:3)]
            #nrows=size(df)[2]-1
            println(propertynames(df))
            #@df df density(:_11, group = (:tot_average, :date), legend = :topleft)
            #@df df density(:tot_average, legend = :topleft)
            @df df density(cols(start:stop), legend = :topleft)
        end

        function denselog(regex::AbstractString,dfs::Vector{DataFrame})
            "selects first match and plots..."
            df = dfs[map(n->occursin(Regex(regex,"i"),n),
                map(x->basename(only(DataFrames.metadata(x))[2]),
                dfs))] |> first
                s = propertynames(df)[Not(end)];
                o = DataFrames.metadata(df)|>collect
                ti = basename(o[1][2])
                @df df density(
                    cols(s),
                    title=ti,
                    yaxis=:log,
                    legend = :topright) 
        end

        function ddense(df::DataFrame)
            s = Symbol.(filter(x->!occursin(r"date|year"i,x),names(df)))
            @df df density(cols(s), legend = :topright)
        end

        function denseplot(df::DataFrame)
            s = propertynames(df)[Not(end)] #masks last column == date     #[1:end-1]
            @df df density(cols(s), legend = :topright)
        end

        function denseplot(df::AbstractString)
            df=readdf(df)
            s = propertynames(df)[Not(end)] #masks last column == date     #[1:end-1]
            @df df density(cols(s), legend = :topright)
        end

        function rplot(reg::Regex, lyr::Int)
            """
            rplot(x::Regex)
            reads first match of regex wasim ncs
            """
            file = filter(x -> occursin(reg,x), readdir(pwd()))
            println("subsetting first nc of $file...")
            file = filter(x -> endswith(x,".nc"), file)|>first
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            Plots.plot(xr[t=lyr];c=cgrad(:thermal),size=(1200*.8, 800*.8))   
        end

        function rplot(reg::Regex)
            """
            rplot(x::Regex)
            reads first match of regex wasim ncs
            """
            file = filter(x -> occursin(reg,x), readdir(pwd()))
            println("subsetting first nc of $file...")
            file = filter(x -> endswith(x,".nc"), file)|>first
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            Plots.plot(xr[t=1];c=cgrad(:thermal),size=(1200*.8, 800*.8))   
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

        function old_waread(x::AbstractString)
            """
            skipping 6 lines - no dropmissing
                for Meteo Time Series
            """
            ms=["-9999","lin","log","-9999.0"]
            df = CSV.read(x,
            DataFrame,
            missingstring=ms,
            ntasks=8,
            skipto=6,
            limit=typemax(Int),
            delim="\t",
            silencewarnings=false,
            normalizenames=true,drop=(i, nm) -> i == 4) #|> dropmissing
            df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
            df=df[:,Not(1:3)]
            metadata!(df, "filename", x, style=:note);
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
            x = glob(x)|>first
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
                if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
                file_path = joinpath(path, file)
            println("reading ",file_path,"...")
            p1 = readdf(file_path)
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
                if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
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

        function vibx(df::String)
            df = readdf(df)
            str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
            @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
            @df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.5,marker=(:black,stroke(1)),legend=false)
        end

        function vibx(df::DataFrame)
            str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            #ln = propertynames(df[end-1])
            ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
            @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
            @df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
        end

        function vio(mm::Regex)
            """
            vioplot wasim timeseries 
            """
            df = readdf(mm)
            str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
            @warn "No basename in metadata!"
            ti = raw""
            end    
            if (any(x->occursin("year",x),names(df)))
                #s = Symbol.(filter(x->!occursin("year",x),names(df)))
                df = df[!,Not("year")]
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
                @df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
                title!(ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
            @df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
            title!(ti)
            #@df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
            end
        end

        function vio(df::DataFrame)
            str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
            @warn "No basename in metadata!"
            ti = raw""
            end    
            if (any(x->occursin("year",x),names(df)))
                #s = Symbol.(filter(x->!occursin("year",x),names(df)))
                df = df[!,Not("year")]
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
                @df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
                title!(ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
            #@df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
            title!(ti)
        #	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
        #	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
        end
        end

        function vio(df::String)
            df = readdf(df)
            str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
            @warn "No basename in metadata!"
            ti = raw""
            end    
            if (any(x->occursin("year",x),names(df)))
                #s = Symbol.(filter(x->!occursin("year",x),names(df)))
                df = df[!,Not("year")]
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
                @df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
                title!(ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
            @df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
            title!(ti)
            #@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.15,legend=false)
            #marker=(:black,stroke(1)),legend=false)
        end
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

        function pline(path::AbstractString)
        #    ms=["-999","-9999","lin","log","LIN","LOG"]
            df = CSV.read(path,DataFrame,
            missingstring="-9999", #also windows
        #    missingstring=ms,
            delim="\t",comment="-",
            silencewarnings=false,
        #    ntasks=4,downcast=true, # got unsupported keyword arguments "ntasks", "downcast" @windows                                          
            normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
            df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
            df=df[:,Not(1:3)]
            # ms=["-9999","lin","log","LIN","LOG","--"] #comment="-",
            # #df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",comment="-",ignorerepeated=true,silencewarnings=true,typemap=Dict(Int64=>String))  |> @dropna() |> DataFrame
            # df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",ignorerepeated=true,silencewarnings=true,typemap=Dict(String=>Int64))
            # df = df[completecases(df), :]

            #Use findall(completecases(df)) to get the indices of the rows.

            # #df = filter( [2]=> x -> !any(f -> f(x), (ismissing)), df)
            # #df = filter( [5]=> x -> isnumeric, df)
            # #parse.(Date, df[:,1:4])
            # #parse.(Date, string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH")
            # df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
            # df=df[:,Not(1:4)]
            nrows=size(df)[2]-1
            st=[]
            for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
            p = make_subplots(rows=nrows, cols=1, 
            shared_xaxes=true, 
            shared_yaxes=false,
            vertical_spacing=0.05,
            #subplot_titles= st;
            )
            for i in 1:nrows;
                    add_trace!(p, 
                    PlotlyJS.scatter(x=df.date, y=df[:,i],
                    name=st[i]),   row=i,     col=1);
            end
            #relayout!(p,height=600*2,width=900*2,title_text="Series of "*basename(path))
            PlotlyJS.relayout!(p,height=600*1.5,width=900*1.5,title_text="Series of "*basename(path))
            display(p)
        end


        function dfplot(df::DataFrame)
            nrows=size(df)[2]-1
            st=[]
            for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
            p = make_subplots(rows=nrows, cols=1, 
            shared_xaxes=true, 
            shared_yaxes=false,
            vertical_spacing=0.05,
            )
            for i in 1:nrows;
                    add_trace!(p, 
                    scatter(x=df.date, y=df[:,i],
                    name=st[i]),   row=i,     col=1);
            end
            relayout!(p,height=600*1.5,width=900*1.5)
            return(p)
        end

        function dfplot(df::AbstractString)
            df=readmeteo(df)
            nrows=size(df)[2]-1
            st=[]
            for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
            p = make_subplots(rows=nrows, cols=1, 
            shared_xaxes=true, 
            shared_yaxes=false,
            vertical_spacing=0.05,
            )
            for i in 1:nrows;
                    add_trace!(p, 
                    PlotlyJS.scatter(x=df.date, y=df[:,i],
                    name=st[i]),   row=i,     col=1);
            end
            PlotlyJS.relayout!(p,height=600,width=900)
            display(p)
        end

        plotdf=dfplot


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
            return (1 - (sum((simulated .- observed).^2) / sum((simulated .- mean(observed)).^2)))
        end

        function kge(df::DataFrame)
            simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
            r = cor(observed, simulated)
            α = std(simulated) / std(observed)
            β = mean(simulated) / mean(observed)
            return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
        end

        function xx(ext::AbstractString)
            path = pwd()
            v = []
            files = readdir(path)
            for file in files
                i = joinpath(path, file)
                if isfile(i) && occursin(Regex(ext),file) && endswith(file, ".nc")
            #println(i);
            i=basename(i)
            osize = stat(i).size
            @printf("%-40s %15.2f MB\n","$(i):",osize/1024^2);
            #outname=replace(i,"nc"=>"jl.png");
                #println(outname," saved!");
            push!(v,i)
                end
            return(v)
        end
        end	

        function cattoyrsum(df::DataFrame)
        ln = Symbol.(filter(x->!occursin("date",x),names(df)))
        df[!, :year] = year.(df[!,:date]);
        #ot=DataFrame[]
        it=[]
        od=(by(df,:year,ln[1]=>sum));
        for i in ln[2:end];
        x=(by(df[Not(:date)],:year,i=>sum));
        push!(it,x[end])      
        end ;
        ot = hcat(od,DataFrame(it))  
        #ot=join(od,x[end],:year,makeuniuqe=true)
        #ot = hcat(od,it)
        return(ot)
        end

        function cattoyrmean(df::DataFrame)
        ln = Symbol.(filter(x->!occursin("date",x),names(df)))
        df[!, :year] = year.(df[!,:date]);
        it=[]
        od=(by(df[Not(:date)],:year,ln[1]=>mean));
        for i in ln[2:end];
        x=(by(df,:year,i=>mean));
        push!(it,x[end])      
        end ;
        ot = hcat(od,DataFrame(it))  
        return(ot)
        end

        function lplot(regex::AbstractString,dfs::Vector{DataFrame})
            "selects first match and plots..."
            df = dfs[map(n->occursin(Regex(regex,"i"),n),
                map(x->basename(only(DataFrames.metadata(x))[2]),
                dfs))] |> first
                ln = Symbol.(filter(x->!occursin("date",x),names(df)))
                nm = propertynames(df)[1:end-1];
                o = DataFrames.metadata(df)|>collect
                ti = basename(o[1][2])
                @df df Plots.plot(:date,cols(ln),yaxis=:log,title=ti)  
        end

        function lplot(df::DataFrame)
            nm = propertynames(df)[1:end-1];
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
            @warn "No basename in metadata!"
            ti = raw""
        end
            ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot(:date,cols(ln),yaxis=:log,title=ti)     
        end


        function lplot(df::String)
            df=readdf(df)
            nm = propertynames(df)[1:end-1];
            o = collect(DataFrames.metadata(df))[1][2] |>basename
            ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot(:date,cols(ln),yaxis=:log,title=o)     
        end

        function lplot(x::Regex)
            df=readdf(x)
            nm = propertynames(df)[1:end-1];
            o = collect(DataFrames.metadata(df))[1][2] |>basename
            ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot(:date,cols(ln),yaxis=:log,title=o)     
        end


        function dfl(regex::Regex)
            "selects first match and plots in log y-axis..."
            df=readf(regex)
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
            ti = raw"" 
            end
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
                @df df Plots.plot(:year,cols(s),yaxis=:log,legend = :topright, title=ti)
            else   
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df Plots.plot(:date,cols(s),yaxis=:log, legend = :topright, title=ti)
                end
        end

        function dfl!(regex::Regex)
            "adds first match and plots in log y-axis..."
            df=readf(regex)
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
            ti = raw"" 
            end
            println("adding $ti")
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
                @df df Plots.plot!(:year,cols(s),yaxis=:log,legend = :topright)
                Plots.annotate!([(20,5,text(ti, 12, :left, :top, :green))])
            else   
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df Plots.plot!(:date,cols(s),yaxis=:log, legend = :topright)
                Plots.annotate!([(20,5,text(ti, 12, :left, :top, :green))])
                end
        end

        function aplot(df::DataFrame)
            df = copy(df)
            df[!,:year]=year.(df[!,:date]) ;
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            o = DataFrames.metadata(df)|>collect
            ti = "AndrewsPlot of "*basename(o[1][2])
            @df df andrewsplot(:year, cols(s), legend = :topleft,title=ti)
        end

        function aplot(df::Regex)
            df = readf(df)
            df = copy(df)
            df[!,:year]=year.(df[!,:date]) ;
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            o = DataFrames.metadata(df)|>collect
            ti = "AndrewsPlot of "*basename(o[1][2])
            @df df andrewsplot(:year, cols(s), legend = :topleft,title=ti)
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
        #    cwd==nothing
        #    dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
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

        function lplot(df::DataFrame)
            nm = propertynames(df)[1:end-1];
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=ti)     
        end

        function lplot(df::String)
            df=readdf(df)
            nm = propertynames(df)[1:end-1];
            o = collect(DataFrames.metadata(df))[1][2] |>basename
            @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=o)     
        end

        function lplotf(df::String)
            df=readdf(df)
            nm = propertynames(df)[1:end-1];
            o = collect(DataFrames.metadata(df))[1][2] |>basename
            @df df plot(:date,cols(nm[1:end-1]),yaxis=:log,title=o)     
        end

        function getnames(ncs::Vector{Raster})
            x=map(x->name(x),ncs)
            return(x)
        end

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
            z = v[broadcast(x->!endswith(x,"nc"),v)];
            return(z)
        end

        function dfon(x1::AbstractString)
            z = filter(file -> occursin(Regex(x1,"i"),file), 
            readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
            return(z)
        end

        function dfonly(x1::Regex)
            z = filter(file -> occursin(x1,file), 
            readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
            return(z)
        end

        function nconly(Any)
            v = readdir();
            z = v[broadcast(x->endswith(x,"nc"),v)];
            return(z)
        end

        function readallras(path::AbstractString)
            v = readdir(path);
            v = v[broadcast(x->endswith(x,"nc"),v)];
            z::Vector{Raster}=[];
            for s in v; 
            #if contains(x1,s) & occursin(r"nc$",s)
            ts=read(Raster(s,missingval=0))
            push!(z,ts);
            end
            return(z)
        end

        function readallras(path::AbstractString, ex::AbstractString)
            v = readdir(path);
            v = v[broadcast(x->endswith(x,"nc") & occursin(ex,x),v)];
            z::Vector{Raster}=[];
            for s in v; 
            #if contains(x1,s) & occursin(r"nc$",s)
            ts=read(Raster(s,missingval=0))
            push!(z,ts);
            end
            return(z)
        end

        function readallras(ex::Regex)
            v = readdir(".");
            v = v[broadcast(x->endswith(x,"nc") & occursin(ex,x),v)];
            z::Vector{Raster}=[];
            for s in v; 
            #if contains(x1,s) & occursin(r"nc$",s)
            ts=read(Raster(s,missingval=0))
            push!(z,ts);
            end
            return(z)
        end

        function denseplot(df::String)
            df=readdf(df)
            s = propertynames(df)[Not(end)]
            @df df density(cols(s), legend = :topright)
        end


        function dfp(mm::Regex)
            """
            plots wasim timeseries
            """
            df=readdf(mm)
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin("year",x),names(df)))
                @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
            end
        end

        function dfp(df::DataFrame)
            #ti = DataFrames.metadata(df)|>only|>last|>basename 
            ti = try
                    DataFrames.metadata(df)|>only|>last|>basename
                catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
                @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
            @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
            end
        end


        function dfp(df::String)
            df=readdf(df)
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin("year",x),names(df)))
                @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
            end
        end

        function dfp(regex::AbstractString,dfs::Vector{DataFrame})
            "selects first match and plots..."
            df = dfs[map(n->occursin(Regex(regex,"i"),n),
                map(x->basename(only(DataFrames.metadata(x))[2]),
                dfs))] |> first
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin("year",x),names(df)))
                @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
            end
        end


        function dfp!(df::DataFrame)
            ti = try
                    DataFrames.metadata(df)|>only|>last|>basename
                catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
                @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
            @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
            end
        end

        function dfp!(df::String)
            df=readdf(df)
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin("year",x),names(df)))
                @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
            end
        end

        function dfp!(mm::Regex)
            """
            plots wasim timeseries
            """
            df=readdf(mm)
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin("year",x),names(df)))
                @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
            end
        end

        function dfp!(regex::AbstractString,dfs::Vector{DataFrame})
            "selects first match and plots..."
            df = dfs[map(n->occursin(Regex(regex,"i"),n),
                map(x->basename(only(DataFrames.metadata(x))[2]),
                dfs))] |> first
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin("year",x),names(df)))
                @df df Plots.plot!(:year,cols(s),legend = :topright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df Plots.plot!(:date,cols(s),legend = :topright, title=ti)
            end
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

        # getf(".*(^th)+.*(nc)+.*")  
        # #SAME
        # getf("^th+.*nc")
        # ###lookbehind	
        # #getf("stack?+.*nc") 
        # #getf("!stack?+.*nc") 

        function plotf(ext::AbstractString)
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
            return(
            dfp(readdf(m[1])))
        end 

        function plotf(ext::String)
            dfp(readdf(ext))
            plot!(title=basename(ext))
        end 

        function plotf(ext::DataFrame)
        dfp(ext)
        end 


        function homg()
            pt="D:/Wasim/Goldbach/revision/"
            cd(pt)
            println("you are here: ",pwd())
        end


        function hombr()
            pt="D:/Wasim/Tanalys/DEM/brend_fab/out/m4/"
            cd(pt)
            println("you are here: ",pwd())
        end

        function hometeo()
            cd("D:/Wasim/Tanalys/DEM/Input_V2/meteo/")
            println("you are here: ",pwd())
        end

        function homreg()
            cd("D:/Wasim/regio/out/");
            println("you are here: ",pwd())
            fd()
        end

        function homg()
            cd("D:/Wasim/Goldbach/");
            println("you are here: ",pwd())
            fd()
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
        function wcl(file::AbstractString)
            open(file) do f
                println(count(_ -> true, eachline(f)))
            end
        end

        function wcl(file::AbstractString,Bool)
            open(file) do f
                ct=(count(_ -> true, eachline(f)))
                #println(file,ct)
                println("$file:\t $ct")
            end
        end

        #wcl(file,true)
        #faster and pure julia:
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
            cd(src_path)
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
            cd(src_path*"/win")
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


        function vgjlrec(snippet::AbstractString)
            """
            recursive grep
            """
            owd=src_path
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

        function vgpyo(snippet::AbstractString)
            owd=pwd()
            cd(dirname(src_path))
            files = filter(file -> endswith(file, ".py"), readdir())
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

        function vgpy(snippet::AbstractString)
            owd=dirname(src_path)
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

        function fdf(df::DataFrame)
            nrows=size(df)[2]-1 
            o = DataFrames.metadata(df)|>collect
            ti = basename(o[1][2])
            fig = PlotlyJS.make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in 1:nrows;
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            fact = 0.7
            logy = true;
            if logy == true
                PlotlyJS.relayout!(fig,
                template="seaborn",
                yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            else
                PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            end
            println(DataFrames.describe(df))
            println("showing plot...")
            display(fig)
        end

        function xdf(df::DataFrame)
            try
                nrows=size(df)[2]-1 
                fig = make_subplots(
                    shared_xaxes=true, 
                    shared_yaxes=true    
                    );
                for i in 1:nrows;
                    add_trace!(fig, 
                    PlotlyJS.scatter(x=df.date, y=df[:,i],
                    name=names(df)[i]));
                end
                PlotlyJS.relayout!(fig,
                template="plotly_dark",
                yaxis_type="log")
                display(fig)
            catch e
                println("An error occurred: ", e)
            finally
                println("showing plot...")
                println(describe(df))
        end
        end

        #go dir up
        function cdb()
            dirname(pwd())|>cd
            pwd()|>println
        end

        #like jdd to vector of strings.
        function fdd(;cwd=pwd())
            dirs = readdir(cwd)
            
            if length(dirs) == 0 
                println("$cwd is empty!")
                return
            end
            
            if filter(x -> (isdir(x)),dirs) == []
                bn = basename(cwd)
                @info "no dirs in $bn !"
                dd()
                return
            end
            
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
            #files = filter(file -> endswith(file, file_ending), readdir())
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
            #pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia";
            pt=src_path;
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

        function median_filter(ras::Raster)
            # Get the array and dimensions of the raster
            Z=Band(1)
            #arr = ras[:Z]
            arr = ras[Z]
            nx, ny = size(arr)
            # Create an output array with the same size and type
            out = similar(arr)
            # Loop over the pixels, excluding the borders
            for i in 2:nx-1, j in 2:ny-1
            # Get the values in the 3x3 window
            window = arr[i-1:i+1, j-1:j+1]
            # Calculate the median of the window
            out[i,j] = median(window)
            end
            # Return a new raster with the filtered array
            return rebuild(ras,out)
        end

        function dfilter(regex::AbstractString,dfs::Vector{DataFrame})
            filter(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
            )
        end

        function filterplot(regex::AbstractString,dfs::Vector{DataFrame})
            "selects first match and plots..."
            df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
            )] |> first
            dfp(df)
        end

        function filterplot!(regex::AbstractString,dfs::Vector{DataFrame})
            "selects first match and add to plot..."
            df = dfs[map(n->occursin(Regex(regex,"i"),n),
            map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
            )] |> first
            y = filter(x->!occursin("date",x),names(df))
            s = map(y -> Symbol(y),y)
            #Symbol(names(df))
            #s = propertynames(df)[Not(end)] #geht auch, aber positionsabhängig
            @df df Plots.plot!(:date,cols(s),legend = :topright)
        end

        function freadfs(ext::AbstractString)
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

        function dfpall(files::Vector{Any})
            "reads, reduces + merges by date and plots"
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
            innerjoin(left, right, on = :date,makeunique=true), 
            dfs)
            y = filter(x->!occursin("date",x), names(df))
            s = map(y -> Symbol(y),y)
            @df df Plots.plot(:date,
                    cols(s),
                    #yaxis = :log,
                    #legend = :bottom)
                    legend = false)
        end


        function dfpall(dfs::Vector{DataFrame})
            "reduces + merges by date and plots"
            df = reduce((left, right) -> 
            innerjoin(left, right, on = :date,makeunique=true), 
            dfs)
            y = filter(x->!occursin("date",x), names(df))
            s = map(y -> Symbol(y),y)
            @df df Plots.plot(:date,
                    cols(s),
                    #yaxis = :log,
                    legend = :bottom)
        end


        function mall(files::Vector{Any})
            "reads, reduces + merges by date"
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
            innerjoin(left, right, on = :date,makeunique=true), 
            dfs)
            return(df)
        end

        function mall(files::Vector{DataFrame})
            "reduces + merges by date"
            df = reduce((left, right) -> 
            innerjoin(left, right, on = :date,makeunique=true), 
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

        cs= [:default		,
            :blues		,
            :bluesreds		,
            :darkrainbow		,
            :darktest		,
            :grays		,
            :greens		,
            :heat		,
            :lightrainbow		,
            :lighttest];



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

        function bardf(x::String)
            "with String"
            df = readdf(x)
            y = filter(x->!occursin("date",x),names(df))
            s = map(y -> Symbol(y),y)
            #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
                ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            df[!, :year] = year.(df[!,:date]);
            df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
            @df df_yearsum Plots.plot(:year,
                cols(s),
                legend = :topright, 
                title=ti,
                seriestype=:bar)
        end

        function bardf(x::Regex)
            "with regex, and new metadata extraction"
            df = readdf(x)
            y = filter(x->!occursin("date",x),names(df))
            s = map(y -> Symbol(y),y)
            #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
                ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            df[!, :year] = year.(df[!,:date]);
            df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
            @df df_yearsum Plots.plot(:year,
                cols(s),
                legend = :topright, 
                title=ti,
                seriestype=:bar)
        end

        function bardf(x::DataFrame)
            "with DataFrame input"
                df = x
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
                @df df_yearsum Plots.plot(:year,
                    cols(s),
                    legend = :topright, 
                    title=ti,
                    seriestype=:bar)
        end

        function bardfm(x::String)
            "with String"
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
            @df df_yearsum Plots.plot(:year,
                cols(s),
                legend = :topright, 
                title=ti,
                seriestype=:bar)
        end

        function bardfm(x::Regex)
            "with regex, and new metadata extraction"
            df = readdf(x)
            y = filter(x->!occursin("date",x),names(df))
            s = map(y -> Symbol(y),y)
            #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
                ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            df[!, :year] = year.(df[!,:date]);
            df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
            @df df_yearsum Plots.plot(:year,
                cols(s),
                legend = :topright, 
                title=ti,
                seriestype=:bar)
        end

        function bardfm(x::DataFrame)
            "with DataFrame input"
                df = x
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
                @df df_yearsum Plots.plot(:year,
                    cols(s),
                    legend = :topright, 
                    title=ti,
                    seriestype=:bar)
        end

        function bardfm!(x::DataFrame)
            "with DataFrame input"
                df = x
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
                @df df_yearsum Plots.plot!(:year,
                    cols(s),
                    legend = :topright, 
                    title=ti,
                    seriestype=:bar)
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

        function readras(file::AbstractString)
            x=read(Raster(file,missingval=0)) #read all in RAM
            #describe(x)
            return(x)
        end

        function readras(mm::Regex)
            "reads first match"
            v::Vector{String} = readdir();
            v = v[broadcast(x->endswith(x,"nc"),v)];
            #file = v[(broadcast(x->occursin(path,x),v))] |>first;
            file = v[(broadcast(x->occursin(mm,x),v))];
            println("reading first match of: ",file)
            f1 = first(file)
            x::Raster = read(Raster(f1,missingval=0))
            return(x)
        end

        function readrasrec(prefix::Regex)
            """
            readras(prefix::Regex)
            reads first match of regex raster
            """
            rootdir="."
            results = []
            for (looproot, dirs, filenames) in walkdir(rootdir)
                for filename in filenames
                    #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                    if (occursin(prefix,filename)) && (endswith(filename,"nc"))
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
            println(results)
            file = first(results)
            x = read(Raster(file,missingval=0)) #read all in RAM
            describe(x)
            return(x)
        end

        function median_filter(ras::Raster)
            # Get the array and dimensions of the raster
            Z=Band(1)
            #arr = ras[:Z]
            arr = ras[Z]
            nx, ny = size(arr)
            # Create an output array with the same size and type
            out = similar(arr)
            # Loop over the pixels, excluding the borders
            for i in 2:nx-1, j in 2:ny-1
            # Get the values in the 3x3 window
            window = arr[i-1:i+1, j-1:j+1]
            # Calculate the median of the window
            out[i,j] = median(window)
            end
            # Return a new raster with the filtered array
            return rebuild(ras,out)
        end

        # function kge_df(ext::String;recursive::Bool)
        #     path = pwd()
        #     files = readdir(path)
        #     v = []
        #     for file in files
        #         file_path = joinpath(path, file)
        #         if isfile(file_path) && endswith(file, ext)
        #             dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
        #             observed  = dd[:,5]
        #             simulated = dd[:,6]
        #             kge_value = kge2(observed, simulated)
        #             nse_value = nse(observed, simulated)
        #             nm = basename(file_path)
        #             println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
        #             printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
        #             push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
        #             v = DataFrame(v)
        #         elseif isdir(file_path) & (recursive==true)
        #             dfs_in_subdir = kge_df(file_path, ext)
        #          end
        #     end
        #     return(v)
        # end

        # function kge_df(path::String,ext::String;recursive=true)
        #     path = pwd()
        #     files = readdir(path)
        #     v = []
        #     for file in files
        #         file_path = joinpath(path, file)
        #         if isfile(file_path) && endswith(file, ext)
        #             dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
        #             observed  = dd[:,5]
        #             simulated = dd[:,6]
        #             kge_value = kge2(observed, simulated)
        #             nse_value = nse(observed, simulated)
        #             nm = basename(file_path)
        #             println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
        #             printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
        #             push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
        #             v = DataFrame(v)
        #         elseif isdir(file_path) & (recursive==true)
        #             dfs_in_subdir = kge_df(file_path, ext)
        #          end
        #     end
        #     return(v)
        # end

        function kge_df(ext::String)
            """
            should be non-recursive
            """
            path = pwd()
            files = readdir(path)
            v = []
            for file in files
                file_path = joinpath(path, file)
                if isfile(file_path) && endswith(file, ext)
                    dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
                    observed  = dd[:,5]
                    simulated = dd[:,6]
                    kge_value = kge2(observed, simulated)
                    nse_value = nse(observed, simulated)
                    nm = basename(file_path)
                    println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
                    printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
                    push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
                    v = DataFrame(v)
                end
            end
            return(v)
        end

        function kge1(simulations, evaluation)
            r = cor(simulations, evaluation)
            α = std(simulations) / std(evaluation)
            β = mean(simulations) / mean(evaluation)
            return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
        end

        function kgerec(ext::String)
            """
            should be recursive
            """
            path::String = pwd()
            results = []
            v = []
            for (looproot, dirs, filenames) in walkdir(path)
                for filename in filenames
                    #if (endswith(filename, ext)) 
                    if (occursin(Regex(ext,"i"),filename)) && 
                        (!occursin(r"output|yrly|nc|png|svg|jpg|xml|ctl",filename))
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
                    sz=length(results)
                    println("found $sz files...\n$results")
                    for file in results
                        dd = CSV.read(file,DataFrame,missingstring="-9999",delim="\t")
                        observed  = dd[:,5]
                        simulated = dd[:,6]
                        kge_value = kge2(observed, simulated)
                        nse_value = nse(observed, simulated)
                        nm = basename(file)
                        printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
                        println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
                        push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm,:path=>file))
                        v = DataFrame(v)
                    end
                    return(v)
                end


        #geht.
        # for (looproot, dirs, filenames) in walkdir(path)
        #     for filename in filenames
        #         if (occursin(Regex(ext,"i"),filename)) && (!occursin(r"yrly|nc|png|svg|jpg|xml|ctl",filename))
        #             push!(results, joinpath(looproot, filename)) 
        #         end
        #     end
        # end

        #map(x->unique())
        #unique(results)
        # results=results[5:8]
        # v = []
        # for file in results
        #     file_path = file
        #         #dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
        #         nm = basename(file_path)
        #         println(nm)
        #         push!(v,Dict(:name=>nm,:dfa=>file))
        #         v = DataFrame(v) # 
        # end


        #ds = kge_df("qout")
        #https://stackoverflow.com/questions/42499528/julia-convention-for-optional-arguments
        #. Mostly you'll define two methods for your function:

        function rp(x::Raster)
            contourf(x; c=cgrad(:matter),size=(1200, 800),
            xlabel="",ylabel="")
        end

        mask_trim(raster, poly) = Rasters.trim(Rasters.mask(raster; with=poly); pad=10)


        function rp(file::AbstractString, lyr::Int)
            #xlyr = length(lyr)!=1 ? 1 : lyr
            ts=read(Raster(file,missingval=0))
            x = ts[t=lyr]
            contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        function rp(file::AbstractString)
            ts=read(Raster(file,missingval=0))
            #ts=Raster(file,missingval=0)
            #contourf(ts; c=cgrad(:thermal),size=(1200, 800))
            contourf(ts; c=cgrad(:matter))
        end

        function rp(x::Regex)
            v = rglob(x)
            file = v[broadcast(x->endswith(x,"nc"),v)]|>first;
            ts=read(Raster(file,missingval=0))
            contourf(ts; c=cgrad(:thermal))
            #contourf(ts; c=cgrad(:thermal),size=(1200, 800))
        end

        function rp(file::AbstractString, lyr::Int)
            ts=read(Raster(file,missingval=0))
            x = ts[t=lyr]
            contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        function rplot(file::AbstractString)
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            Plots.plot(xr;c=cgrad(:thermal),
            xlabel="",
            ylabel="",
            size=(1200*.66, 800*.66))
        end

        function rplot(file::AbstractString, lyr::Int)
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            Plots.plot(xr[t=lyr];c=cgrad(:thermal),size=(1200*.8, 800*.8))   
        end
            #xr[t=20,cname="RdBl"]|>plot
            # c=:thermal]|>Plots.plot

        function rplot(file::Raster, lyr::Int)
            xr = file
            Plots.plot(xr[t=lyr];c=cgrad(:thermal),size=(1200*.8, 800*.8))
        end

        function rplot(x::Raster;ex::Int)
            """
            subset raster by last dimension.
            excludelayer and plot the rest.
            rplot(tras;ex=2)
            rplot(tras,2)
            keyword arguments by using a semicolon (;) in the parameter list. 
            """
            xr = x[Dim{Rasters.name(x.dims)[end]}(Rasters.Where(x -> x >= ex))]
            Plots.plot(xr;
                c=cgrad(:thermal),
                size=(1200*.66, 800*.66),
                xlabel="",
                ylabel="",
                title=Rasters.name(xr))
        end

        #tras[Dim{Rasters.name(tras.dims)[end]}(Rasters.Where(x -> x >= 6))] |>Plots.plot

        function rp(x::Regex;msk::Float64,gt=false)
            """
            subset raster by mask and last dimension.
            excludelayer and plot the rest.
            keyword arguments by using a semicolon (;) in the parameter list. 
            """
            #msk=0 #msk::Float64
            v = rglob(x)
            x = v[broadcast(x->endswith(x,"nc"),v)]|>first
            x = read(Raster(x,missingval=0;lazy=true))
            #lastdim = x.data|>size|>last
            #lastdim = Int(x.dims[3][end])
            #rn = r[t=2:ee];    #subs
            #dimname = (Rasters.name(x.dims)[end]) #gehtnet
            #x = x[t=1]
            #x = x[Rasters.Where(k)=lastdim]
            #xs = x[lastdim=Rasters.Where(Symbol.(Rasters.name(x.dims)[end]))]
            xs = x[t=Int(x.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            Plots.contourf(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:thermal),
                xlabel="",
                ylabel="",
                size=(1200, 800))
        end



        function rpall(file::AbstractString)
            xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            #xr = read(Raster(file;crs=EPSG(25832),missingval=-9999))
            Plots.plot(xr;c=cgrad(:thermal),
            xlabel="",
            ylabel="",
            size=(1200*.8, 800*.8))
        end



        function cntplt(file::AbstractString)
            x=read(Raster(file,missingval=0))
            Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        #describe(ar[2])  

        function cpl(file::Raster)
            #x=Rasters.rebuild(file;missingval=-9999)
            x=Rasters.rebuild(file;missingval=0)
            #x=x[t=1]
            Plots.contourf(x; 
            c=cgrad(:matter),
            xlabel="",
            ylabel="")
        end

        function cntplt(file::Union{Missing, Float64})
            x=file
            Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        function cplt(file::AbstractString)
            x=read(Raster(file))
            Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
        end

        function cpal(ext::AbstractString)
            path = pwd()
            files = readdir(path)
            for file in files
                i = joinpath(path, file)
                if isfile(i) && occursin(Regex(ext),file) && (!occursin("stack",file)) && endswith(file, ".nc")
            outname=replace(i,"nc"=>"jl.png");
                #println(outname)
                r=read(Raster(i,missingval=0));
                p=Plots.contourf(r;
                title=replace(basename(i),".nc"=>""), #split(outname,"/")[end], #basename(i)
                c=cgrad(:thermal),
                size=(1200, 800));
                savefig(p,outname)
                println(basename(outname)," saved!");
                end
        end
        end

        function stackplot(ext::AbstractString)
            path = pwd()
            files = readdir(path)
            for file in files
                i = joinpath(path, file)
                if isfile(i) && occursin(Regex(ext),file) && (occursin("stack",file)) && endswith(file, ".nc")
            outname=replace(i,"nc"=>"jl.png");
                #println(outname)
                r=read(Raster(i,missingval=0,mappedcrs=EPSG(25832)));
            #(i,missingval=-9999,mappedcrs=EPSG(25832))
            ee = Int(r.dims[3][end])
            rn = r[t=2:ee];    #subset till end
                p=Plots.plot(rn;
        #		title=replace(basename(i),".nc"=>""), #no title cause problems
                c=cgrad(:thermal),
                size=(1200, 800));
                savefig(p,outname)
                println(basename(outname)," saved!");
                end
        end
        end

        #########tdiff 
        function tdifnc()
            """
            path="."
            read non-recursivley and plots tdiff
            glob<->rglob
            stores tdiff.nc

            """
            pot = filter(x->occursin(regand("Layer","nc"),x),glob(r"etp"))|>last|>readras
            real= filter(x->occursin(regand("Layer","nc"),x),glob(r"etr"))|>last|>readras
            td = pot-real
            nm = filter(x->occursin(regand("Layer","nc"),x),glob(r"etp"))|>last
            ti = split(nm,".")|>x->x[end-1]
            p = Plots.plot(td;
                title = "Tdiff[mm] of "*ti,
                c=cgrad(:matter))
            display(p)
            ##for overwriting force
            #write("tdiff.nc",td;force=true)
            @warn raw"for writing: write(tdiffjl.nc,RasterObject;force=true) "
            return(td)
        end

        function ncmean(x::Regex)
            """
            path="."
            read non-recursivley and plots tdiff
            glob<->rglob
            stores to mean.nc

            """
            #x=r"win"
            file = filter(file -> occursin(x,file) && endswith(file,".nc"), 
                readdir())|>first
            td::Raster = read(Raster(file,missingval=0))

        #    nm = filter(x->occursin(regand("Layer","nc"),x),glob(r"etp"))|>last
        #    ti = split(nm,".")|>x->x[end-1]
            # plot(td;
            # title = "Tdiff[mm] of "*ti,
            # c=cgrad(:matter))
            ##for overwriting force
            td = td/365
            outname = basename(file)
            m = match(r".*[.]",basename(file))
            outfile = contains(basename(file),".") ? string(m.match,"nc") : basename(file)*".html"
            # replace!(outname,".nc","mean.nc")
            # replace!(outname,"sum","")    
            write("mean.nc",td;force=true)
        end


        function facets_loop(ext::AbstractString)
            """
            like stackplot, but for interactive view
            """
            path = pwd()
            files = readdir(path)
            for file in files
                i = joinpath(path, file)
                if isfile(i) && occursin(Regex(ext,"i"),file) && (occursin("stack",file)) && endswith(file, ".nc")
                    @warn("subsetting first layer...")
                    r=read(Raster(i,missingval=0,mappedcrs=EPSG(25832)));
                    ee = Int(r.dims[3][end])
                    rn = r[t=2:ee];    #subset till end
                    p=Plots.plot(rn;
                    xlabel="",
                    ylabel="",
                    c=cgrad(:thermal),
                    size=(1200, 800));
                    display(p)
                end
        end
        end

        function facets(ext::AbstractString)
            """
            like stackplot, but for interactive view
            """
            #grids = filter(x -> occursin(regand(ext,"nc"),x), readdir())
            grids = filter(x -> occursin(Regex(ext,"i"),x), readdir())
            grids = filter(x->endswith(x,"nc"),grids)
            if (length(grids)>0)
                file = first(grids)
                @warn("taking first match of $grids\n -> $file")
            else
                dir=pwd()
                @warn "No netcdf match found in $dir !"
                return
            end
                if isfile(file)
                    r=read(Raster(file,missingval=0,mappedcrs=EPSG(25832)));
                    if (r.dims[end]|>length == 1)
                        @warn("only one layer available...")
                        p=Plots.plot(r;
                        # xlabel="",
                        # ylabel="",
                        #c=cgrad(:thermal),
                        c=cgrad(:matter),
                        size=(1200, 800));
                        display(p)
                    else
                        @warn("subsetting first layer...")
                        ee = Int(r.dims[3][end])
                        rn = r[t=2:ee];    #subset till end
                        p=Plots.plot(rn;
                        xlabel="",
                        ylabel="",
                        #title=replace(basename(i),".nc"=>""), #no title cause problems
                        c=cgrad(:thermal),
                        size=(1200, 800));
                        display(p)
                #savefig(p,outname)
                #println(basename(outname)," saved!");
                end
        end
        end

        function dfyrs(df::DataFrame;)
            ti = DataFrames.metadata(df)|>only|>last|>basename
            fact,logy = 1,0
            y = filter(x->!occursin("date",x),names(df))
            s = map(y -> Symbol(y),y)
            df[!, :year] = year.(df[!,:date]);
            df = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
            #df = df[!,Not("date")]
            ln = Symbol.(filter(x->!occursin("year",x),names(df)))
            nrows=size(df)[2]-1
            if nrows == 1
                ln = only(ln)
                fig = 
                PlotlyJS.plot(
                PlotlyJS.scatter(x=df.year, y=df[!,ln],
                name=ln,type="bar")
                );
                PlotlyJS.relayout!(fig,
                    height=600*fact,width=900*fact,
                    title_text="Series of "*ti)
            else
                fig = PlotlyJS.make_subplots(
                    shared_xaxes=true, 
                    shared_yaxes=true    
                    );
                for i in ln
                    PlotlyJS.add_trace!(fig, 
                    PlotlyJS.scatter(x=df.year, y=df[:,i],
                    name=i));
                end
                if logy == true
                    PlotlyJS.relayout!(fig,yaxis_type="log",
                    height=600*fact,width=900*fact,
                    title_text="Series of "*ti)
                else
                    PlotlyJS.relayout!(fig,
                    height=600*fact,width=900*fact,
                    title_text="Series of "*ti)
                end
            end
            display(fig)
        end

        function dfpjs(df::DataFrame;)
            nrows=size(df)[2]-1 
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            fig = PlotlyJS.make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in 1:nrows;
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            fact,logy = 1,0
            if logy == true
                PlotlyJS.relayout!(fig,yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti,
                xaxis=PlotlyJS.attr(    rangeslider_visible=true,rangeselector=PlotlyJS.attr(
                buttons=[
                    PlotlyJS.attr(count=1, label="1m", step="month", stepmode="backward"),
                    PlotlyJS.attr(count=6, label="6m", step="month", stepmode="backward"),
                    PlotlyJS.attr(count=1, label="YTD", step="year", stepmode="todate"),
                    PlotlyJS.attr(count=1, label="1y", step="year", stepmode="backward"),
                    PlotlyJS.attr(step="all")
                ]    )))
            else
                PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti,
                xaxis=PlotlyJS.attr(rangeslider_visible=true,
            rangeselector=PlotlyJS.attr(        buttons=[
                PlotlyJS.attr(count=1, label="1m", step="month", stepmode="backward"),
                    PlotlyJS.attr(count=6, label="6m", step="month", stepmode="backward"),
                    PlotlyJS.attr(count=1, label="YTD", step="year", stepmode="todate"),
                    PlotlyJS.attr(count=1, label="1y", step="year", stepmode="backward"),
                    PlotlyJS.attr(step="all")
                ]    )))
            end
            display(fig)
        end

        function dfpjs(df::String;)
            df = readf(df)
            nrows=size(df)[2]-1 
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            fig = PlotlyJS.make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in 1:nrows;
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            fact,logy = 1,0
            if logy == true
                PlotlyJS.relayout!(fig,yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            else
                PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            end
            display(fig)
        end

        function dfpjs(df::Regex;)
            df = readdf(df)
            nrows=size(df)[2]-1 
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            fig = PlotlyJS.make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in 1:nrows;
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            fact,logy = 1,0
            if logy == true
                PlotlyJS.relayout!(fig,yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            else
                PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            end
            display(fig)
        end

        function dfbarjs(df::Regex;)
            df = readdf(df)
            nrows=size(df)[2]-1
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            fig = PlotlyJS.make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in 1:nrows;
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.bar(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            fact,logy = 0.66,0
            if logy == true
                PlotlyJS.relayout!(fig,
                template="seaborn",
                yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            else
                PlotlyJS.relayout!(fig,
                template="seaborn",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            end
            display(fig)
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

        function barp(x::DataFrame)
            "with DataFrame input"
                df = x
                    ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
                if any(x->occursin("year",x),names(df))
                    ln = Symbol.(filter(x->!occursin("year",x),names(df)))
                    @df df Plots.plot(:year,
                        cols(ln),
                        legend = :topright, 
                        title=ti,
                        seriestype=:bar) #color=:lightrainbow
                elseif any(x->occursin("month",x),names(df))
                    ln = Symbol.(filter(x->!occursin("month",x),names(df)))
                    @df df Plots.plot(:month,
                        cols(ln),
                        legend = :topright, 
                        title=ti,
                        seriestype=:bar)
                elseif (
                    any(x->occursin("month",x),names(df)) & 
                    any(x->occursin("year",x),names(df))            
                    )
                    ln = (filter(x->!occursin("month",x),names(df)))
                    ln = Symbol.(filter(x->!occursin("year",x),ln))
                    @df df Plots.plot(:month,
                        cols(ln),
                        legend = :topright, 
                        title=ti,
                        seriestype=:bar)
                else
                    dfp(df)        
                end
        end


        #"D:\Wasim\Tanalys\DEM\brend_fab\out\m6\hhydfab.m6.2012"

        function vg2(regex::AbstractString, ending::AbstractString)
            cmd = `grep --color=always -C2 -rIHn -E "$regex" --include="*.$ending"`
            run(cmd)
        end

        # cd("/mnt/c/Users/Public/Documents/Python_Scripts/julia")
        # vg2("readallras","jl")

        function wslp(winpt::AbstractString)
            printstyled("needRAW like raw",color=:red)
            #return(raw$winpt)
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

        # function wslp(winpt::AbstractString)
        #     cmd = `$file = $winpt.Path; echo $file`
        #     run(cmd)
        # end

        # function wslp(winpt::AbstractString)
        #     cmd = `$file = $pwd.Path|sed -e 's,.*:,/mnt/d,' -e 's/[\\]/\\/\/g';`
        #     run(cmd)
        # end

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

        function dfsplog(dfs::Vector{DataFrame};save="")
            "plots and adds"
            df = dfs[1]
            s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
            p = @df df Plots.plot(:date,
                    cols(s),
                    yaxis = :log,
                    legend = false)
                    #legend = :bottom)
            for i in 2:length(dfs)
                nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
                println("adding $nm")
                s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
                @df dfs[i] Plots.plot!(:date,cols(s),
                label="$nm") #geht, wenn oben legend true ist.
                # label="$nm",
                # legend = false)
                # Plots.annotate!(0.5, 0.5, text(nm, 14))
            end
            display(p)
            if !isempty(save) 
                Plots.savefig(p,save*".png")
                printstyled("$save saved as $save*.png! \n",color=:green)
            end
        end

        function dfsp(dfs::Vector{DataFrame};save="")
            "plots and adds"
            df = dfs[1]
            s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
            p = @df df Plots.plot(:date,
                    cols(s),
                    legend = false)    #legend = :bottom)
            for i in 2:length(dfs)
                nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
                println("adding $nm")
                s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
                @df dfs[i] Plots.plot!(:date,cols(s),
                label="$nm") #geht, wenn oben legend true ist.
            end
            display(p)
            if !isempty(save) 
                Plots.savefig(p,save*".png")
                printstyled("$save saved as $save*.png! \n",color=:green)
            end
        end

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

        # function zzg(snippet::AbstractString)
        #     """
        #     greps from repl_history
        #     """
        #     #file = raw"/home/ubu/.julia/logs/repl_history.jl"
        #     file = raw"C:\Users\chs72fw\.julia\logs\repl_history.jl"
        #         open(file) do f
        #             counter = 0 # Zähler initialisieren
        #             before  = counter - 1
        #             after   = counter + 2
        #             for line in eachline(f)
        #                 counter += 1 # Zähler erhöhen
        #                 if (counter==before) || (counter==after)
        #                     #printstyled("$counter:\t",color=:light_red) 
        #                     printstyled("$file:\t",color=:light_magenta,bold=false) 
        #                     printstyled("$line\n",color=:yellow,bold=true) 
        #                 end
        #                 if contains(line,snippet)
        #                     printstyled("$counter:\t",color=:light_red) 
        #                     printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
        #                     printstyled("$line\n",color=:green,bold=true)

        #                 end
        #             end
        #         end
        # end


        function ll(; reg::Bool=false,
            x::AbstractString="nc")
            """
            lists from current dir and optionally greps iRegex
            """
            #readdir(join ? pwd() : ".", join=join, sort=sort)
            reg ? filter(file -> occursin(Regex(x,"i"),file), readdir()) : readdir()
        end
        # ll(;reg=true,x="temp")
        # ll(;reg=false)
        # ll(;reg=true)    #<- here only nc files
        # ll()
        # @edit readdir()


        ##
        second(x) = x[2]
        third(x) = x[3]
        lastbefore(x) = x[end-1]
        # geht auch:
        # v |> x -> x[2]
        # rglob(r"pr")|>lastbefore     
        #rglob(r"pr")|>lastbefore|>dfread|>describe   
        # rglob(r"pr")|>lastbefore|>yrmean

        #p()
        #cd("D:/Wasim/regio/out/c8")

        function reorder_df(df::DataFrame)
            """
            date to last position
            """
            df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
            return(df)
        end


        function tree(cwd::AbstractString, prefix=" ")
            paths::Vector{String} = []
            #cwd = abspath(cwd)
            cwd = relpath(cwd)
            for (looproot, dir, filenames) in walkdir(cwd)
                for relpath in dir
                    push!(paths, joinpath(looproot,relpath))
                end
            end
            if length(paths)==0
                printstyled("no subfolders present in \n"*pwd(),color=:red)
                return
            end
            ap=abspath(cwd)
            printstyled("$ap\n",color=:red)
            for relpath in paths
                prefix = " "
                        relpath = replace(relpath, r"[^\\]"=> ":",count=1)
                        ###helper to fast cdinto....
                        relpath = replace(relpath, r"[^\\]*."=>"-- ",count=1)
                        relpath = replace(relpath, "\\"=> "/")
                println(relpath)
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
                    @warn("error! ")
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
                    @warn("error! ")
                    # Skip files that can't be loaded as a DataFrame
                    continue
                end
            end
            return(outdf)
        end

        #ds=qall()

        function theplot(x::AbstractString)
            df = DataFrame(CSV.File(x))
            nm=names(df)[end-1] #lastbefore column (qobs)
            ##subset DF by value (all positive vals..)
            df = filter(nm => x -> x > 0, df)
            df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
            ndf = df[!,Not(1:4)]
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
            r2 = overall_pearson_r^2
            #nse(simulations, evaluation)
            nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
            kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
            #ti = "Time Series of $(uppercase(first(split(basename(x), '-'))))"
            ti = first(split(basename(x),"_"))
            #subs = "Pearson r: $(round(overall_pearson_r, digits=2))\nPearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
            subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
            #p = plot(title=[ti, subs], ylabel="[unit/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
            p = Plots.plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
            Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
            Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
            Plots.annotate!(
            #nrow(ndf), 0.95*maximum(ndf.Observed),
            :bottomright,
            text("$subs", 10, :black, :right))
            return p
        end

        function theplot(x::DataFrame)
            ndf = copy(x)
            @warn "renaming to :Simulated,:Observed,:Date !"
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
            r2 = overall_pearson_r^2
            nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
            kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
            ti = try
                basename(last(collect(DataFrames.metadata(ndf)))[2])
            catch
                @warn "No basename in metadata!"
                raw""
            end 
            subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
            p = Plots.plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
            Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
            Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
            Plots.annotate!(
            :bottomright,
            text("$subs", 10, :black, :right))
            return p
        end

        ftp(z::AbstractString) = theplot(first(
            filter(x->occursin(Regex(z,"i"),x),
            filter(x->endswith(x,"qoutjl"),readdir())))
            )

        #filter(x->endswith(x,"qout"),readdir()))))

        function ftp(df::DataFrame)
            if size(df)[2]!=3
                #throw(@warn "wrong number of columns - using dfp!")
                @warn "wrong number of columns - using dfp!"
                @warn "need :Simulated,:Observed,:Date !"
                display(dfp(df))
                return
            end
            ndf = copy(df)
            # rename!(ndf,3=>"date")
            # @warn "last col renamed! !"
            #reorder
            ndf = hcat(ndf[!,Not(Cols(r"date"i))],ndf[:,Cols(r"date"i)])
            #nm=names(ndf)[end-2] #lastbefore date column (qobs)
            ##subset DF by value (all positive vals..)
            #ndf = filter(nm => x -> x > 0, ndf)
            #filter([1]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), ndf)
            ndf = filter(:date=> x -> !any(f -> f(x), (ismissing, isnothing)), ndf)
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            dropmissing!(ndf) ##hmm sketchy..
            overall_pearson_r = cor(ndf[!, :Simulated],ndf[!, :Observed])
            r2 = overall_pearson_r^2
            #nse(simulations, evaluation)
            nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
            kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
            ti = try
                basename(last(collect(DataFrames.metadata(ndf)))[2])
            catch
            @warn "No basename in metadata!"
                raw""
            end 
            subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
            p = plot(ndf[!, :Date], ndf[!, :Simulated], 
            title=ti, 
            line=:dash, color=:blue, label="Modeled",
            ylabel="[mm/day]", xlabel="modeled time", 
            yscale=:log2, legend=:topleft)
            plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
            annotate!(
            :bottomright,
            text("$subs", 10, :black, :right))
            return p
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

        #for i in range(1,last(size(r)));println(describe(r[t=i]));end 
        #for i in range(1,(size(r)[end]));println(describe(r[t=i]));end
        function descr(r::Raster)
            nm=name(r);
        for i in 1:last(size(r));
            printstyled("$nm t $i \n",color=:green)
            describe(r[:,:,i])
        end
        end

        function readalloutput()
            """
            reads timeseries and stores to Vector{DataFrame}
            reads NetCDFs and stores to Vector{Any}
            dfs,ncs = readalloutput()
            filterplot("qg",dfs)
            filterplot("rad",ncs)
            """
            cwd = "."
            dfs=loadalldfs(cwd)
            ncs=readallras(cwd)
            return(dfs,ncs)
        end

        function stats(r::Raster)
            m = mean(r) # get the mean for each band
            n = minimum(r) # get the minimum for each band
            x = maximum(r) # get the maximum for each band
            d = median(r) # get the median for each band
            s = std(r) # get the standard deviation for each band
            # get the number of missing values for each band
            
            c = try
                parse(Float64, Rasters.missingval(r)) 
                catch
                @warn "No missval in metadat! -set to 0.0"
                c = 0.0
                end
            
            
            arr=[m,n,x,d,s,c]'
            #println("$nm\n",arr)
            df = DataFrame(arr,:auto)
            nm=["mean", "min", "max", "median", "sd", "missval"]
            rename!(df,nm)

            # Matrix(arr) # convert the adjoint to a matrix
            # m = collect(arr) # convert the adjoint to a matrix
            # xc=[
            #     "mean", "min", "max", "median", "sd", "missval",
            # m,n,x,d,s,c]
            
            return(df)
        end

        function rp3(x::String)
            """
            3D plot with geoarrays
            """
            ga = GeoArrays.read(x)
            values = ga.A # a 3D array of raster values
            #GeoArrays.coords(ga) # a tuple of x, y and band coordinates
            #crs = ga.crs # a string of CRS definition
            t = GeoArrays.coords(ga)|>size
            coords = (1:t[1], 1:t[2]) # a Tuple{UnitRange{Int64}, UnitRange{Int64}}
            # Plots.surface(coords[1], coords[2], values[:, :, 1]) # plot the first band
            # xlabel!("x")
            # ylabel!("y")
            # zlabel!("value")
            # ti=basename(x) 
            # title!(ti)
            ti=basename(x)     #title!("3D Raster Plot")
            #p1=
            Plots.surface(coords[1], coords[2], 
            values[:, :, 1]    ,
            xlabel="x",ylabel="y",zlabel="value",title=ti)
            #display(p1)
            
            # Plots.surface(coords[1], coords[2], values[:, :, 1]    ,
            # xlabel="x",ylabel="y",zlabel="value",
            # xaxis=:log,
            # title=ti)
        end

        function gofbatch()
            println("batch R Script for GOF")
            arr = filter(x -> isfile(x) && endswith(x, "_qout") && !occursin(r"\.(png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly)$", x), readdir())
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

        function globdf(x::Regex)
            """
            greps ts from current dir Regex
            """
            filter(file -> (occursin(x,file) & 
            (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
            ), readdir())
        end

        #globdf(r"g")
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

        #"wa"|>globdf

        function ldfpall(x::Regex)
            "reads, reduces + merges by date and plots log y-axis"
            files = rglob(x)
            dfs = DataFrame[]
            for file in files
                if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
                file_path = file
            println("reading ",file_path,"...")
            p1 = readdf(file_path)
            ##renamer
            for x in 1:size(p1,2)-1
                rename!(p1,x=>basename(file_path)*names(p1)[x])
            end
            push!(dfs, p1)
                end
            end
            df = reduce((left, right) -> 
            innerjoin(left, right, on = :date,makeunique=true),dfs)
            ##to preserve column order and names (date at last position)
            #df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
            y = filter(x->!occursin("date",x), names(df))
            s = map(y -> Symbol(y),y)
            @df df Plots.plot(:date,
                    cols(s), yaxis = :log,
                    legend = :outertopright)
        end

        function dfl(x::DataFrame)
            "selects first match and plots in log y-axis..."
            df=x
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
            @warn "No basename in metadata!"
            ti = raw"" 
            end
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
                @df df Plots.plot(:year,cols(s),yaxis=:log,legend = :topright, title=ti)
            else   
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df Plots.plot(:date,cols(s),yaxis=:log10, 
                legend = :topright, title=ti)
                end
        end

        function rpr(x::String)
            """
            3D plot with Rasters
            """
            #x="D:/temp/saale/output/thulba/v2/All_HydrologicResponseUnits.nc"
            ga = Rasters.read(Rasters.Raster(x))
            values = ga.data[:,:,1]
            #coords = (1:length(lookup(ga, X)), 1:length(lookup(ga, Y)))
            s=size(ga)
            coords = (1:s[2], 1:s[1]) #swiched on netdf!
            ti=basename(x)
            #p1=
            Plots.surface(        coords[1], coords[2], 
                values    ,   xlabel="x",
                ylabel="y",zlabel="rastervalue",title=ti)
        end

        function rpr(x::Raster)
            """
            3D plot with Rasters
            x .> 0 
            """
            msk = float(0.001)
            xs = x[t=Int(x.dims[3][end])+1]
            zm = (xs .> msk)
            # #zm = (gt) ? (xs .> msk) : (xs .< msk)
            # x = Rasters.mask(xs; with=zm); 
            # values = x.data[:,:,1] .> 0
            # s=size(x)
            # coords = (1:s[2], 1:s[1]) #swiched on netdf!
            ti=Rasters.name(x)
            M=Rasters.mask(xs; with=zm);
            #M=Rasters.rebuild(M;missingval=msk)
            #M=Rasters.rebuild(M;missingval=msk)|>skipmissing
            Plots.surface(
                M,
                #legend=:outerbottomleft,
                legend=:bottom,
                #legend= false,
                xlabel=" ",
                ylabel=" ",
                zlabel=" ",
                title=ti,
                camera = (-20, 75)        
                )
            # Plots.surface(coords[1], coords[2], 
            #     values    ,   
            #     xlabel="x",
            #     ylabel="y",zlabel=ti,
            #     #title=ti,
            #     camera = (-20, 75))
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
            arr = filter(x -> isfile(x) && endswith(x, "qoutjl") && !occursin(r"\.(png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly)$", x), readdir())
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

        # ot = tdiff()
        # wawrite("tdiff-jl.txt",ot)

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

        # function lat(directory::String)
        #     """
        #     list_files_sorted_by_last_change
        #     """
        #     files = readdir(directory)
        #     file_times = Dict{String, Dates.DateTime}()
        #     for file in files
        #         file_path = joinpath(directory, file)
        #         stat_info = stat(file_path)
        #         file_times[file] = Dates.unix2datetime(stat_info.mtime)
        #     end
        #     sorted_files = sort(files, by=file -> file_times[file], rev=true)
        #     return DataFrame(name=sorted_files,last_modified=file_times[file])
        # end

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

        #globf(r"rout")|>npp

        loadall = readalloutput

        function readfall(x::AbstractString)
            """
            all ts to dfs
            """
            path = pwd()
            files = glob(x)
            dfs = DataFrame[]
            for file in files
                if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
                file_path = joinpath(path, file)
            println("reading ",file_path,"...")
            #p1 = readdf(file_path)
            p1 = readf(file_path) #new
            push!(dfs, p1)
                end
            end
            return(dfs)
        end

        # dfs = readfall("v4")

        # filterplot("qbas",dfs)
        # filterplot("wi",dfs)
        # filterplot("pre",dfs)
        # filterplot!("qges",dfs)


        function tpjs(x::AbstractString)
            """
            theplot optimized to PlotlyJS
            """
            df = DataFrame(CSV.File(x))
            nm=names(df)[end-1] #lastbefore column (qobs)
            ##subset DF by value (all positive vals..)
            df = filter(nm => x -> x > 0, df)
            df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
            ndf = df[!,Not(1:4)]
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
            r2 = overall_pearson_r^2
            nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
            kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
            ti = first(split(basename(x),"_")) 
            #subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
            subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"

            p = Plots.plot(
            ndf[!, :Date], ndf[!, :Simulated], color=:red, label="Modeled",
            title=ti, ylabel="[mm/day]", xlabel="modeled time", 
            yscale=:log10, 
            legend=:outerbottomleft)

            #plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
            plot!(p, ndf[!, :Date], ndf[!, :Observed],line=:dash, color=:blue, 
            label="Observed")
            annotate!(
                #:bottomright,
                last(ndf.Date), mean(ndf[!,1]),
                #last(ndf.Date), 0.95*maximum(ndf[!,1]),
                text("$subs", 10, :black, :right))
                #nrow(ndf), 0.95*maximum(ndf.Observed),
                #nrow(ndf),minimum(ndf.Observed),
            return p #or!
            # fact=.6
            # PlotlyJS.relayout!(p,
            # template="seaborn",
            # height=600*fact,width=900*fact,
            # title_text="Series of "*ti)
            # display(p)
        end

        #tpjs("Wechterswinkel-qoutjl")

        function tpjs(x::DataFrame)
            """
            theplot optimized to PlotlyJS
            """
            ndf = x
            nm=names(ndf)[2]     #obscol
            ##subset DF by value (all positive vals..)
            ndf = filter(nm => x -> x > 0, ndf)
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            dropmissing!(df)
            overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
            r2 = overall_pearson_r^2
            #qpl(ndf)
            nse_score = nse(ndf)
            kge_score = kge(ndf)
            ti = try
                basename(last(collect(DataFrames.metadata(ndf)))[2])
            catch
            @warn "No basename in metadata!"
                raw""
            end 
            subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"
            p = Plots.plot(
            ndf[!, :Date], ndf[!, :Simulated], color=:red, 
            label="Modeled",
            title=ti, ylabel="[mm/day]", xlabel="modeled time", 
            yscale=:log10, 
            legend=:outerbottomleft);
            #legend=:outertopleft);
            plot!(p, ndf[!, :Date], ndf[!, :Observed],
            line=:dash, color=:blue, 
            label="Observed")
            annotate!(
                #:bottomright, #nope
                last(ndf.Date), mean(ndf[!,2]),
                #last(ndf.Date), 0.95*maximum(ndf[!,1]),
                text("$subs", 10, :black, :right))
            return p
        end

        function tpjs(x::Regex)
            """
            theplot optimized to PlotlyJS
            """
            ndf = waread(x)
            dropmissing!(ndf)
            #nm=names(ndf)[2]     #obscol
            #ndf = filter(nm => x -> x > 0, ndf)
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            
            overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
            r2 = overall_pearson_r^2
            nse_score = nse(ndf)
            kge_score = kge(ndf)
            ti = try
                basename(last(collect(DataFrames.metadata(ndf)))[2])
            catch
            @warn "No basename in metadata!"
                raw""
            end 
            subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"
            p = Plots.plot(
            ndf[!, :Date], ndf[!, :Simulated], color=:red, 
            label="Modeled",
            title=ti, ylabel="[mm/day]", xlabel="modeled time", 
            yscale=:log10, 
            legend=:outerbottomleft);
            plot!(p, ndf[!, :Date], ndf[!, :Observed],
            line=:dash, color=:blue, 
            label="Observed")
            annotate!(
                last(ndf.Date), mean(ndf[!,2]),
                text("$subs", 10, :black, :right))
            return p
        end

        readalldfs = loadalldfs
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
        #name1="Brend_Brend_24431002_Unterweissenbrunn_"
        #replace(name1, r"(\w+)[ _]\1" => s"\1")
        #replace(name1, r"(\w+)[ _]\1" => s"\1" , "_"=>"-")

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

        #  # f = glob("so")
        # tff2(f[1:3])
        # tff2(f[5:10])

        function ssup()
            include("C:/Users/Public/Documents/Python_Scripts/julia/win/smallfuncs.jl")
        end


        function filterplot(regex::AbstractString,ncs::Vector{Raster})
            "selects first match and plots..."
            r = ncs[map(n->occursin(Regex(regex,"i"),n),
            map(x->string.(name(x)),ncs)
            )] |> first
            plot(trim(r),xlab="",ylab="",title=name(r))
            cz=r.metadata.val|>collect|>last|>last
            Plots.annotate!(:bottomright,
            text("cellsize: $cz", 8, :black, :right))
        end

        function filterplot!(regex::AbstractString,ncs::Vector{Raster})
            "selects first match and add to plot..."
            r = ncs[map(n->occursin(Regex(regex,"i"),n),
            map(x->string.(name(x)),ncs)
            )] |> first
            plot!(trim(r),xlab="",ylab="",title=name(r))
        end

        function facets(r::Raster)
            """
            like stackplot, but for interactive view
            """
            # if (length(grids)>0)
            #     file = first(grids)
            #     @warn("taking first match of $grids\n -> $file")
            # else
            #     dir=pwd()
            #     @warn "No netcdf match found in $dir !"
            #     return
            # end
            #        if isfile(file)
            #         r=read(Raster(file,missingval=0,mappedcrs=EPSG(25832)));
                    if (r.dims[end]|>length == 1)
                        @warn("only one layer available...")
                        p=Plots.plot(r;
                        # xlabel="",
                        # ylabel="",
                        #c=cgrad(:thermal),
                        c=cgrad(:matter),
                        size=(1200, 800));
                        display(p)
                    else
                        @warn("subsetting first layer...")
                        ee = Int(r.dims[3][end])
                        rn = r[t=2:ee];    #subset till end
                        p=Plots.plot(rn;
                        xlabel="",
                        ylabel="",
                        #title=replace(basename(i),".nc"=>""), #no title cause problems
                        c=cgrad(:thermal),
                        size=(1200, 800));
                        display(p)
                end
        #end
        end

        function rhist(r::Raster)
            M=r.data[:,:,1]
            Plots.histogram(skipmissing(M),title=r.name,legend=false)
        end

        function gofbatch_nosvg()
            println("batch R Script for GOF")
            arr = filter(x -> isfile(x) && endswith(x, "qoutjl") && !occursin(r"\.(png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly)$", x), readdir())
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

        function ftplin(df::DataFrame)
            if size(df)[2]!=3
                #throw(@warn "wrong number of columns - using dfp!")
                @warn "wrong number of columns - using dfp!"
                @warn "need :Simulated,:Observed,:Date !"
                display(dfp(df))
                return
            end
            ndf = copy(df)
            # rename!(ndf,3=>"date")
            # @warn "last col renamed! !"
            #reorder
            ndf = hcat(ndf[!,Not(Cols(r"date"i))],ndf[:,Cols(r"date"i)])
            nm=names(ndf)[end-2] #lastbefore date column (qobs)
            ##subset DF by value (all positive vals..)
            ndf = filter(nm => x -> x > 0, ndf)
            #filter([1]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), ndf)
            ndf = filter(:date=> x -> !any(f -> f(x), (ismissing, isnothing)), ndf)   
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            dropmissing!(ndf) ##hmm sketchy..
            overall_pearson_r = cor(ndf[!, :Simulated],ndf[!, :Observed])
            r2 = overall_pearson_r^2
            #nse(simulations, evaluation)
            nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
            kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
            ti = try
                basename(last(collect(DataFrames.metadata(ndf)))[2])
            catch
            @warn "No basename in metadata!"
                raw""
            end 
            subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
            p = plot(ndf[!, :Date], ndf[!, :Simulated], 
            title=ti, 
            line=:dash, color=:blue, label="Modeled",
            ylabel="[mm/day]", xlabel="", 
            legend=:topleft)
            plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
            annotate!(
            :bottomright,
            text("$subs", 10, :black, :right))
            return p
        end

        function kge_fread()
            kge_read(pwd(),"outjl");
        end

        function plotlybaryr(df::DataFrame)
            p = PlotlyJS.plot(df, kind = "bar");
            #rename(df,replace(names(df),"_1"=>""))
            s = (filter(x->!occursin(r"year|date",x),names(df)))
            #renamer - remove chars   
            for x in s
                newname=replace(x,"_1"=>"")
                #println(newname)
                #rename!(df, Dict(:i => "A", :x => "X"))
                rename!(df,Dict(x=>newname))
            end
            
            s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
            
            for i in s;
                PlotlyJS.add_trace!(p, 
                PlotlyJS.bar(x=df.year, y=df[:,i],
                name=i)       );
            end
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            PlotlyJS.relayout!(p,
            template="seaborn",
            #height=600*1.5,width=900*1.5, #dirname(pwd()
            #title_text="Series of "*basename(path)
            title_text=ti)
            display(p)
        end

        baryrjs = plotlybaryr

        function baryrsum(df::DataFrame)
            """
            automatically sums if only datecolumn is available
            """
            v = map(
                (x->occursin(r"date", x) & !occursin(r"year", x)),
                (names(df))
                )

            #v = map(x->occursin(r"date", x),(names(df)))
            
            if any(v)
                df = yrsum(df) 
                # inplace yrsum:
                # y = filter(x->!occursin("date",x),names(df))
                # df[!, :year] = year.(df[!,:date])
                # df = DataFrames.combine(groupby(df, :year), 
                #     y .=> sum .=> y)
            end
            
            s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
            ti = try
                z=DataFrames.metadata(df)|>only|>last|>basename
                basename(pwd())*" $z"
            catch
                @warn "No basename in metadata!"
                ti = "Series of "*basename(pwd())
            end
            @df df groupedbar(df.year,cols(s), 
            legend = :outertopright,
            xticks = df.year,
            xrotation = 45,
            xlabel = "", ylabel = "[mm]", title = ti)
        end

        function baryrmean(df::DataFrame)
            """
            automatically sums if only datecolumn is available
            """
            v = map(
                (x->occursin(r"date", x) & !occursin(r"year", x)),
                (names(df))
                )

            if any(v)
                # #df = yrsum(df), but inplace
                # y = filter(x->!occursin("date",x),names(df))
                # df[!, :year] = year.(df[!,:date])
                # df = DataFrames.combine(groupby(df, :year), 
                #     y .=> sum .=> y)
                df = yrmean(df)
            end
            
            s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
            ti = try
                z=DataFrames.metadata(df)|>only|>last|>basename
                basename(pwd())*" $z"
            catch
                @warn "No basename in metadata!"
                ti = "Series of "*basename(pwd())
            end
            @df df groupedbar(df.year,cols(s), 
            legend = :outertopright,
            xticks = df.year,
            xrotation = 45,
            xlabel = "", ylabel = "[mm]", title = ti)
        end

        function baryr(df::DataFrame)
            s = filter(x -> !(occursin(r"year|date", x)), names(df))
            for x in s
                newname = replace(x, "_1" => "")
                rename!(df, Dict(x => newname))
            end
            s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
            # Create a grouped bar plot for each column    
            # StatsPlots.groupedbar(df.year,[df[!, col] for col in s], 
            # #group = s,#repeat(s, outer = size(df, 1)),
            # xlabel = "Year", ylabel = "Value", title = "Grouped Bar Plot")
            ti = try
                z=DataFrames.metadata(df)|>only|>last|>basename
                basename(pwd())*" $z"
            catch
                @warn "No basename in metadata!"
                ti = "Series of "*basename(pwd())
            end
            @df df groupedbar(df.year,cols(s), 
            #group = s,#repeat(s, outer = size(df, 1)),
            legend = :outertopright,
            xticks = df.year,
            xrotation = 45,
            xlabel = "", ylabel = "[mm]", title = ti)
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

        function kge_df3()
            """
            should be non-recursive
            """
            x1=r"qoutjl"
            files = filter(file -> occursin(x1,file),
                readdir()[broadcast(x->!endswith(x,r"nc|txt"),
                readdir())])
            v = []
            for file_path in files
                if isfile(file_path)
                    df = waread(file_path);
                    dropmissing!(df)
                    observed  = df[:,1]
                    simulated = df[:,2]
                    kge_value = kge1(simulated,observed)
                    nse_value = nse(simulated,observed)
                    nm = basename(file_path)
                    println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
                    printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
                    push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
                    v = DataFrame(v)
                end
            end
            return(v)
        end


        function dsbar(ds::DataFrame)
            ds.name=map(x->replace(x,r"-qoutjl*" => ""),ds.name)
            ann = map(x->string.(round(x;sigdigits=1)),ds.KGE)
            bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
                title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
                fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
                annotations = (ds.name,ds.KGE, ann, :top),
                bar_width = 0.6)
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


        function dpr(x::Regex)
            """
            correlation plots on dataframe
            """
            df = globdf(x)|>first|>waread
            df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
            v = (names(df[!,1:2]))
            a = reshape(v, 1, 2)
            Plots.plot(df.date,[df[!,1], df[!,2]], 
            label=a, 
            xlabel="Date", ylabel="[mm/day]",legend = :topleft)
            r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
            # annotate!(last(df.date), 0.85*maximum(df[!,1]),
            # text("R² = $r2", 10, :black, :right))
            kge = round(kge2(df[!,1], df[!,2]), digits=2)
            nse_value = round(nse(df[!,1], df[!,2]), digits=2)
            annotate!(
                #:topright,
                last(df.date), 0.95*maximum(df[!,1]),
            text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
            )
        end

        function dpr!(x::Regex)
            """
            correlation plots on dataframe
            """
            df = globdf(x)|>first|>waread
            df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
            v = (names(df[!,1:2]))
            a = reshape(v, 1, 2)
            Plots.plot(df.date,[df[!,1], df[!,2]], 
            label=a, 
            xlabel="Date", ylabel="[mm/day]",legend = :topleft)
            r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
            #annotate!(last(df.date), 0.85*maximum(df[!,1]),
            #text("R² = $r2", 10, :black, :right))
            kge = round(kge2(df[!,1], df[!,2]), digits=2)
            nse_value = round(nse(df[!,1], df[!,2]), digits=2)
            annotate!(
                #:topright,
                last(df.date), 0.95*maximum(df[!,1]),
            text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
            )
        end

        function dpr(x)
            df=readdf(x)
            v = (names(df[!,1:2]))
            a = reshape(v, 1, 2)
            Plots.plot(df.date,[df[!,1], df[!,2]], 
            label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
            r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
            # annotate!(last(df.date), 0.85*maximum(df[!,1]),
            # text("R² = $r2", 10, :black, :right))
            kge = round(kge2(df[!,1], df[!,2]), digits=2)
            nse_value = round(nse(df[!,1], df[!,2]), digits=2)
            annotate!(
                #:topright,
                last(df.date), 0.95*maximum(df[!,1]),
            text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
            )
        end

        function dpr!(x)
            df=readdf(x)
            v = (names(df[!,1:2]))
            a = reshape(v, 1, 2)
            Plots.plot(df.date,[df[!,1], df[!,2]], 
            label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
            r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
            # annotate!(last(df.date), 0.85*maximum(df[!,1]),
            # text("R² = $r2", 10, :black, :right))
            kge = round(kge2(df[!,1], df[!,2]), digits=2)
            nse_value = round(nse(df[!,1], df[!,2]), digits=2)
            annotate!(
                #:topright,
                last(df.date), 0.95*maximum(df[!,1]),
            text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
            )
        end


        function dpr(x::DataFrame)
            """
            correlation plots on dataframe
            """
            df = copy(x)
            df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
            v = (names(df[!,1:2]))
            a = reshape(v, 1, 2)
            Plots.plot(df.date,[df[!,1], df[!,2]], 
            label=a, 
            xlabel="Date", ylabel="[mm/day]",legend = :topleft)
            r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
            # annotate!(last(df.date), 0.85*maximum(df[!,1]),
            # text("R² = $r2", 10, :black, :right))
            kge = round(kge2(df[!,1], df[!,2]), digits=2)
            nse_value = round(nse(df[!,1], df[!,2]), digits=2)
            annotate!(
                #:topright,
                last(df.date), 0.95*maximum(df[!,1]),
            text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
            )
        end

        function dpr!(x::DataFrame)
            """
            correlation plots on dataframe
            """
            df = copy(x)
            df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
            v = (names(df[!,1:2]))
            a = reshape(v, 1, 2)
            Plots.plot(df.date,[df[!,1], df[!,2]], 
            label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
            r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
            # annotate!(last(df.date), 0.85*maximum(df[!,1]),
            # text("R² = $r2", 10, :black, :right))
            kge = round(kge2(df[!,1], df[!,2]), digits=2)
            nse_value = round(nse(df[!,1], df[!,2]), digits=2)
            annotate!(
                #:topright,
                last(df.date), 0.95*maximum(df[!,1]),
            text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
            )
        end


        # function rhist(r::Raster{Union{Missing, Float64}, 3})
        #     M=Array(r.data|>collect)
        #     Plots.histogram(skipmissing(M),title=r.name,legend=false)
        # end

        function ftsp(x::AbstractString)
            nc = NCDataset(x);
            #nc.attrib
            dict = nc|>Dict   
            mykeys = keys(dict)
            println(string.(mykeys))
            time = nc["time"][:]
            v = filter(x->!occursin(r"time|lon|lat|x|y|spatial_ref",x),string.(mykeys))|>first
            xm = nc[v]|>size|>first
            xm = Int(round(median(1:xm);digits=0))
            ym = nc[v]|>size|>second
            ym = Int(round(median(1:ym);digits=0))
            plot(time, nc[v][xm,ym,:],
            label=nc[v].attrib["units"],
            title=nc[v].attrib["long_name"])        
        end

        function nctodf(x::AbstractString)
            nc = NCDataset(x);
            #nc.attrib
            dict = nc|>Dict   
            mykeys = keys(dict)
            #println(string.(mykeys))
            v = filter(x->!occursin(r"time|lon|lat|x|y",x),string.(mykeys))|>first
            time = nc["time"][:]
            datetime_vector = coalesce.(time, missing)
            #df = hcat(nc[v][end,end,:],datetime_vector)
            xm = nc[v]|>size|>first
            xm = Int(round(median(1:xm);digits=0))
            ym = nc[v]|>size|>second
            ym = Int(round(median(1:ym);digits=0))
            df = DataFrame(
                    v=>nc[v][xm,ym,:],      #x, y, indices
                    "date"=>datetime_vector)
            # df = DataFrame(
            #         v=>nc[v][end,end,:],      #x, y, indices
            #         "date"=>datetime_vector)
            DataFrames.metadata!(df, "filename", x, style=:note);        
            #df.date = Date.(string.(df.x2),"yyyy-mm-ddTHH:MM:SS") #not needed
                # plot(time, nc[v][end,end,:],
                    # label=nc[v].attrib["units"],
            # title=nc[v].attrib["long_name"])        
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

        #dfs = jldfnm(globdf(r"so"))
        #dfs = jldf(globdf(r"so"))

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





        # function umlauts(input_file::AbstractString, output_file::AbstractString)
        #     # Read input file
        #     data = readdlm(input_file, String)

        #     # Apply string replacements
        #     for i in 1:size(data, 1)
        #         data[i] = replace(data[i], r"ß" => "ss")
        #         # data[i] = replace(data[i], r"\/" => "_")
        #         # data[i] = replace(data[i], r"_" => "-")
        #         # data[i] = replace(data[i], r",," => "")
        #         # data[i] = replace(data[i], r"\xc4" => "Ae")
        #         # data[i] = replace(data[i], r"\xd6" => "Oe")
        #         # data[i] = replace(data[i], r"\xdc" => "Ue")
        #         # data[i] = replace(data[i], r"\xe4" => "ae")
        #         # data[i] = replace(data[i], r"\xf6" => "oe")
        #         # data[i] = replace(data[i], r"\xfc" => "ue")
        #         # data[i] = replace(data[i], r"\xdf" => "ss")
        #     end

        #     # Write modified data back to the file
        #     output_file = open(output_file, "w")
        #     for i in 1:size(data, 1)
        #         println(output_file, data[i])
        #     end
        #     close(output_file)
        # end
        
        #umlauts(input_file, output_file)


        macro bash_str(s) open(`bash`,"w",stdout) do io; print(io, s); end;end
        #bash""" which python """

        macro pwrs_str(s) open(`powershell -noprofile`,"w",stdout) do io; print(io, s); end;end
        #pwrs""" which python """
        #pwrs""" pwd """
        # function op()
        #     pwrs""" explorer . """
        # end 

        macro pwp_str(s) open(`powershell`,"w",stdout) do io; print(io, s); end;end
        #pwp""" fdm """
        #cdb()

        macro cmd_str(s) open(`cmd \c`,"w",stdout) do io; print(io, s); end;end
        #cmd""" pwd """
        ##nope, bad idea
        # run(`cmd \c pwd";" exit`)

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

        function rp3(x::Raster)
            """
            3D plot with Rasters from Raster.
            """
            #t = x|>size
            #coords = (1:t[1], 1:t[2])
            ti= try 
                basename(x);
            catch 
                @warn "no basename available, trying to parse name..."; 
                ti=Rasters.name(x); 
            end
            x.dims    
            Plots.surface(x[:, :, 1],
            xlabel="x",ylabel="y",zlabel="value",title=ti)
        end

        function rpmcf(x::Regex;msk::Float64,gt=false)
            """
            subset raster by mask and last dimension.
            excludelayer and plot the rest.
            keyword arguments by using a semicolon (;) in the parameter list. 
            rpm(r"rad";msk=0.0,gt=true) ->plots
            """
            #msk=0 #msk::Float64
            v = glob(x)
            r = v[broadcast(x->endswith(x,"nc"),v)]|>first
            #r = nconly(x)|>first

            #first("∀ϵ≠0: ϵ²>0", 1)

            r = read(Raster(r,missingval=0;lazy=true))
            xs = r[t=Int(r.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            fact=0.5
            Plots.contourf(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:thermal),
                xlabel="",
                ylabel="",
                size=(1200*fact, 800*fact))
        end

        function rpm(x::Regex;msk::Float64,gt=false)
            """
            subset raster by mask and last dimension.
            excludelayer and plot the rest.
            keyword arguments by using a semicolon (;) in the parameter list. 
            rpm(r"rad";msk=0.0,gt=true) ->plots
            """
            #msk=0 #msk::Float64
            v = glob(x)
            r = v[broadcast(x->endswith(x,"nc"),v)]|>first
            #r = nconly(x)|>first

            #first("∀ϵ≠0: ϵ²>0", 1)

            r = read(Raster(r,missingval=0;lazy=true))
            xs = r[t=Int(r.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            fact=0.5
            Plots.plot(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:thermal),
                xlabel="",
                ylabel="",
                size=(1200*fact, 800*fact))
        end


        function maskplot(xs::Raster;msk::Float64,gt=true)
            """
            subset raster by mask and plot.
            no subset on last dim!
            keyword arguments by using a semicolon (;) in the parameter list. 
            maskplot(r;msk=0.0) ->plots
            """
            #msk=0 #msk::Float64
            #xs = r[t=Int(r.dims[3][end])+1]
            zm = (gt) ? (xs .> msk) : (xs .< msk)
            fact=0.5
            Plots.plot(
                Rasters.mask(xs; with=zm); 
                c=cgrad(:matter),
                xlabel="",
                ylabel="",
                size=(1200*fact, 800*fact))
        end


        function baryrsum(df::Regex)
            """
            automatically sums if only datecolumn is available
            """
            df = globdf(df)|>first|>readf
            v = map(
                (x->occursin(r"date", x) & !occursin(r"year", x)),
                (names(df))
                )

            if any(v)
                df = yrsum(df) 
            end
            
            s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
            ti = try
                z=DataFrames.metadata(df)|>only|>last|>basename
                basename(pwd())*" $z"
            catch
                @warn "No basename in metadata!"
                ti = "Series of "*basename(pwd())
            end
            @df df groupedbar(df.year,cols(s), 
            legend = :outertopright,
            xticks = df.year,
            xrotation = 45,
            xlabel = "", ylabel = "[mm]", title = ti)
        end

        function baryrmean(df::Regex)
            """
            automatically sums if only datecolumn is available
            """
            df = globdf(df)|>first|>readf
            v = map(
                (x->occursin(r"date", x) & !occursin(r"year", x)),
                (names(df))
                )

            if any(v)
                df = yrmean(df)
            end
            
            s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
            ti = try
                z=DataFrames.metadata(df)|>only|>last|>basename
                basename(pwd())*" $z"
            catch
                @warn "No basename in metadata!"
                ti = "Series of "*basename(pwd())
            end
            @df df groupedbar(df.year,cols(s), 
            legend = :outertopright,
            xticks = df.year,
            xrotation = 45,
            xlabel = "", ylabel = "[mm]", title = ti)
        end

        function kgeval()
            """
            kge barplot 
            """
            ds = kge_df3()
            ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.name)
            ann = map(x->string.(round(x;sigdigits=1)),ds.KGE)
            Plots.bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
            title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
            fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
            annotations = (ds.name,ds.KGE, ann, :top),
            bar_width = 0.6)
        end

        function nseval()
            """
            nse barplot with values > 0
            """
            ds = kge_df3()
            ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.name)
            dfi = filter(row -> row.NSE .> 0, ds)
            ann = map(x->string.(round(x;sigdigits=1)),dfi.NSE)
            Plots.bar(dfi.name, dfi.NSE, xlabel = "Name", ylabel = "NSE", legend = false, 
                title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
                fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
                annotations = (dfi.name,dfi.NSE, ann, :top),
                bar_width = 0.6)        
        end

        function nsevalraw()
            """
            nse barplot with all values
            """
            ds = kge_df3()
            ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.name)
            dfi = ds
            ann = map(x->string.(round(x;sigdigits=1)),dfi.NSE)
            Plots.bar(dfi.name, dfi.NSE, xlabel = "Name", ylabel = "NSE", legend = false, 
                title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
                fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
                annotations = (dfi.name,dfi.NSE, ann, :top),
                bar_width = 0.6)        
        end
        #    using Images
        #img=load("waba-jl.png")
        function waba()
            wpth="C:/Users/Public/Documents/Python_Scripts/julia/water-balance.jl"
            include(wpth)
            #img=load("waba-jl.png")
            yd=waread("waba-input.wa")|>yrsum
            @warn "try baryr(yd) !"
            baryr(yd)|>display
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
        macro gl(s) glob(s)|>first;end
        #fastplot
        macro fp(s) dfp(Regex(s));end
        macro flog(s) dfl(Regex(s));end
        macro ncrm() ncrem="C:/Users/Public/Documents/Python_Scripts/julia/ncremover.jl";include(ncrem);end
        macro nco(s) nconly(s);end
        macro dfo(s) first(dfonly(s));end

        function read_soildata(data::String)
            # remove { and } characters from the data
            data = broadcast(x -> replace(x,    
                r"method" => "",
                r"MultipleHorizons" => "",
                r"}" => ""), data)
            
            # split the data into lines
            lines = split(data, '\n')
            
            # initialize an array to store the results
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
                    if occursin(" = ", field)
                        # split the field into key and value
                        key, value = split(field, " = ")
                        
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
                            elseif occursin(r"^\d+ \d+", value)
                                # convert the value to an array of integers
                                value = parse.(Int, split(value))
                            elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)
                                # convert the value to an array of floats
                                value = parse.(Float64, split(value))
                            end
                        end
                        
                        # store the key-value pair in the dictionary
                        data[key] = value
                    end
                end
                
                # append the dictionary to the result array
                push!(result, data)
            end
            
            return result
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
            !occursin(r"sh|csv|html|xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file)
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
                for df in dfs
                    x=collect(DataFrames.metadata(df))[1][2]|>basename
                    println(x)
                end
            catch
                @error "no metadata in $df !"
                return
            end
        end

        # cdb()
        # dfs=qbb()
        # getnames(dfs)
        # #@vv "append"
        # #write("test",dfs)
        # #df=dfs[5]

        # for df in dfs
        #     x=collect(DataFrames.metadata(df))[1][2]|>basename
        #     println("writing $x to ",x*".csv")
        #     #hdr = "QBB of "*joinpath(pwd(),x)
        #     hdr = "QBB of $x \n"
        #     onam = x*".csv"
        #     write(onam, hdr)
        #     CSV.write(onam,df,
        #     header=false,
        #     #quotestrings=false,
        #     #openquotechar=Char(' '),
        #     #closequotechar=Char(' '),
        #     transform = (col, val) -> something(val, missing),
        #     append=true,
        #     delim=" ")
        # end

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

        function nsx(dfi::DataFrame)
            """
            nse barplot with all values
            """
            dfi.name=map(x->replace(x,r"-qoutjl.*" => ""),dfi.name)
            ann = map(x->string.(round(x;sigdigits=3)),dfi.NSE)
            Plots.bar(dfi.name, dfi.NSE, xlabel = "Name", ylabel = "NSE", legend = false, 
                title = "x", xrotation = -25, fmt = :png, size = (800, 600), 
                fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
                annotations = (dfi.name,dfi.NSE, ann, :top),
                xtickfont = font(6),    
                bar_width = 0.77)        
        end

        function ctlg(;dir_path::String, file_ending::String, match::String)
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

        function rmqout()
            map(x->rm(x),glob("qoutjl"))
        end

        #ctlg(pwd(),"txt","weiss")

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
            """
            filter(row -> any(occursin(Regex(x, "i"), 
                string(value)) for value in row), 
                    eachrow(df))
        end

        function agjson(jsonfile::AbstractString)
            """
            reads json and transforms points...
            but AG transformation is wrong :(


            ERROR: type UnionAll has no field read
            works only in WSL
            """
            #const AG = ArchGDAL

            geojson_file=jsonfile
            jsonbytes = read(geojson_file) # read the geojson file as bytes
            fc = GeoJSON.read(jsonbytes)
            pts=[]
            for geom in fc.geometry
                xc = [(x) for x in geom]|>first|>first
                yc = [(x) for x in geom]|>first|>last
                pt = AG.createpoint(xc,yc)
                pt = AG.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
            # df = [] ##geht auch
            # for x in pts
            #     x1 = AG.getpoint(x,0)
            #     tmp = DataFrame(x=[x1[1]], y=[x1[2]])
            #     push!(df,tmp)
            # end
            # df = reduce(vcat, df) 
            #Plots.plot(df.x, df.y, seriestype=:scatter)
            df = DataFrame( 
                "x" => [AG.getpoint(x,0)[1] for x in pts],
                "y" => [AG.getpoint(x,0)[2] for x in pts] )
            return(df)
        end

        function wslpath()
            # Run the `wslpath` command to convert the current directory to a WSL path
            wsl_cmd = `wsl wslpath -a $(pwd)  `
            wsl_path = readchomp(pipeline(wsl_cmd))
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

        function readras2(file::AbstractString)
            x=read(Raster(file,missingval=-9999.000000)) #read all in RAM
            return(x)
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

            #@gl "mo"
        end

        function qplot(df::DataFrame)
            """
            takes first two cols of df and plots r2 QQ
            """
            if any(map(x->contains(x,"date"),names(df)))
                df = df[!,Not(:date)]
            end
            if ncol(df)>2
                df = df[!,1:2]
            end
            r2 = round(cor(df[!,1], df[!,2])^2, digits=3)
            p = df |> x -> qqplot(x[!,1],x[!,2], 
                #title = "R² = "*string(ti),
                qqline = :fit)
                #color = :grays) # erstellt ein QQ-Diagramm <- black
            xlabel!(p,names(df)[1])
            ylabel!(p,names(df)[2])
            annotate!(p,:bottomright, text("R² = "*string(r2), :black))
            #xr = maximum(df[!,1])
            #annotate!(xr-0.5, 2.0, text("R² = $ti", 12))
        end

        function fzplot(x::AbstractString)
            #import ArchGDAL
            r = x
            r = Raster(r;missingval=-9999)
            r = r ./ 3600|>Rasters.trim
            msk = float(0.001)
            #msk = Float32(0)
            zm = (r .> msk)
            r = Rasters.mask(r; with=zm); 
            r = Rasters.rebuild(r,missingval=minimum(r));
            plot(r,c=:blues,title="flowtime in [h]",xaxis="",yaxis="")
        end

        function fzplot2(x::AbstractString)
            """
            AG only.
            import ArchGDAL as AG
            """
            r = AG.readraster(x)
            r = AG.getband(r,1)
            dx = r ./ 3600
            dx = permutedims(dx, (2, 1))
            dx = reverse(dx, dims=1)
            Plots.contour(dx,title="flowtime in [h]",xaxis="",yaxis="")
            #c=:blues,
        end

        function dfpjs(df::DataFrame;)

            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end

            s = (filter(x->!occursin(r"year|date",x),names(df)))
            #renamer - remove char _   
            for x in s
                newname=replace(x,"_"=>" ")
                rename!(df,Dict(x=>newname))
            end
            s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))

            ##make scores
            overall_pearson_r = cor(df[!,2], df[!,1])
            r2 = overall_pearson_r^2
            nse_score = nse(df)
            kge_score = kge(df)

            subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"

            fig = PlotlyJS.make_subplots(shared_xaxes=true, shared_yaxes=true)

            for i in s
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:, i], name=i)
                )
            end

            fact = .66
            PlotlyJS.relayout!(fig,
                template="seaborn",
                #template="simple_white",
                height=650*fact,
                width=1200*fact,
                title_text=ti,
                xaxis_rangeslider_visible=true,
                annotations=[attr(
                            text=subs,
                            #x=minimum(df[!,2]),
                            y=maximum(df.date),
                            xanchor="right",
                            yanchor="bottom",
                            xref="x",
                            yref="y",
                            showarrow=false,
                            bordercolor="#c7c7c7",
                            borderwidth=2,
                            borderpad=4,
                            bgcolor="#ff7f0e",
                            opacity=0.6
                        )],
                updatemenus=[
                    Dict(
                        "type" => "buttons",
                        "direction" => "left",
                        "buttons" => [
                            Dict(
                                "args" => [Dict("yaxis.type" => "linear")],
                                "label" => "Linear Scale",
                                "method" => "relayout"
                            ),
                            Dict(
                                "args" => [Dict("yaxis.type" => "log")],
                                "label" => "Log Scale",
                                "method" => "relayout"
                            )
                        ],
                        "pad" => Dict("r" => 1, "t" => 10),
                        "showactive" => true,
                        "x" => 0.11,
                        "xanchor" => "left",
                        "y" => 1.1,
                        "yanchor" => "auto"
                    ),
                ]
                )

                display(fig)

        end

        function dfbar(df::DataFrame)
            ti = try
                    DataFrames.metadata(df)|>only|>last|>basename
                catch
                @warn "No basename in metadata!"
                ti = raw""
            end
            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
                #@df df Plots.bar(:year,cols(s),legend = :topright, title=ti)
                @df df groupedbar(df.year,cols(s), legend = :outertopright, title=ti)
            else    
            s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
            #@df df Plots.bar(:date,cols(s),legend = :topright, title=ti)
            @df df groupedbar(df.date,cols(s), legend = :outertopright, title=ti)
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

        function vio(df::DataFrame)
            """
            with mon abbr. see ovio for numbered months
            """   
                ti = try
                    DataFrames.metadata(df)|>only|>last|>basename
                catch
                @warn "No basename in metadata!"
                ti = raw""
                end    
                str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
                month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
                if (any(x->occursin("year",x),names(df)))
                    df = df[!,Not("year")]
                    s = Symbol.(filter(x->!occursin("date",x),names(df)))
                    @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
                    xticks!(0.5:11.5 , month_abbr)
                    title!(ti)
                else    
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
                xticks!(0.5:11.5 , month_abbr)
                title!(ti)
            end
        end


        function mbx(df::DataFrame)
            """
            annotated boxplot with mean
            """
            str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            ln = Symbol.(filter(x->!occursin(r"date|year|month"i,x),names(df)))
            month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
            p = @df df StatsPlots.boxplot(str,cols(ln),
                fillalpha=0.75, 
                linewidth=0.25,
                notch = true,
                whisker_width = :match,
                legend=false)
            xticks!(0.5:11.5 , month_abbr)
            df.Month = month.(df.date)
            means = DataFrames.combine(groupby(df,:Month ), ln[1] => mean)
            #mean(means[!,2])
            for i in eachrow(means)
                m = i[2]
                annotate!(i.Month - 0.5, m, #+ 1 
                text(round(m; digits=2), 6, :center, :top))
            end
            display(p)
        end

        function dprbig(x::Regex)
            """
            correlation plots on dataframe
            size=(1200,800)
            """
            df = globdf(x)|>first|>waread
            df = hcat(df[!,Not(Cols(r"date"i))],df[:,Cols(r"date"i)])
            v = (names(df[!,1:2]))
            a = reshape(v, 1, 2)
            Plots.plot(df.date,[df[!,1], df[!,2]], 
            label=a, 
            xlabel="Date", ylabel="[mm/day]",
            legend = :topleft,
            size=(1200,800)
            )
            r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
            kge = round(kge2(df[!,1], df[!,2]), digits=2)
            nse_value = round(nse(df[!,1], df[!,2]), digits=2)
            annotate!(
                last(df.date), 0.95*maximum(df[!,1]),
            text("KGE = $kge\nNSE = $nse_value\nR² = $r2", 10, :black, :right)
            )
        end

        function qplot(x::Vector{Float64},y::Vector{Float64})
            # if any(map(x->contains(x,"date"),names(df)))
            #     df = df[!,Not(:date)]
            # end
            # if ncol(df)>2
            #     df = df[!,1:2]
            # end
            r2 = round(cor(x, y)^2, digits=3)
            p = qqplot(x,y,
                qqline = :fit)
            # xlabel!(p,names(df)[1])
            # ylabel!(p,names(df)[2])
            annotate!(p,:bottomright, text("R² = "*string(r2), :black))
        end

        function agheat(s::AbstractString;step=30,lyr=1,msk=0.001)
            """
            AG only.
            import ArchGDAL as AG errors...
            plots heatmap from raster with ArchGDAL
            """
            if !isfile(s)
                error("file not found!")
            end
            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF !")
            end
            r = ArchGDAL.readraster(s)
            #r = AG.readraster(x)
            r = ArchGDAL.getband(r,lyr)
            #nd = AG.getnodatavalue(r)
            dx = reverse(r, dims=2)
            dx = reverse(dx, dims=1)
            println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
            bitmat = dx .> msk
            # Filter the dx matrix using bitmat
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]

            dx = dx_output
            heatmap_plot = heatmap(dx, c=:matter,
                            title=basename(s), 
                            xaxis="", yaxis="")
            step = step
            for i in 1:step:size(dx, 1)
                for j in 1:step:size(dx, 2)
                    value = round(dx[i, j]; digits=2)
                    color = isnan(value) ? :white : :black
                    annotate!(j, i, Plots.text(string(value), 7, color, :center, 
                        halign=:center, rotation=-35.0))
                end
            end
            # Show the plot
            display(heatmap_plot)
        end

        function agcont(s::AbstractString;lyr=1,msk=0.001)
            """
            AG only.
            import ArchGDAL as AG
            plots raster with ArchGDAL
            """
            if !isfile(s)
                error("file not found!")
            end
            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF !")
            end
            
            r = ArchGDAL.readraster(s)
            println(ArchGDAL.getdriver(r))
            dx = ArchGDAL.getband(r,lyr)
            ArchGDAL.setnodatavalue!(dx, -9999.0)
            dx = reverse(dx, dims=2)
            dx = reverse(dx, dims=1)
            println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
            #println(extrema(dx))
            bitmat = dx .> msk
            # Filter the dx matrix using bitmat
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]
            dx = dx_output
            #Plots.contour(dx,c=:matter,title=basename(s),xaxis="",yaxis="")
            Plots.contourf(dx, fill=true, levels=10,
                #c=:matter,
                title=basename(s),xaxis="",yaxis="")

        end

        function agcont2(s::AbstractString;step=30,lyr=1,msk=0.001)
            """
            AG only.
            import ArchGDAL as AG
            agcont2(f;lyr=5,msk=-1)
            agcont2(f;lyr=2,msk=-1,step=50)
            """
            if !isfile(s)
                error("file not found!")
            end
            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF !")
            end
            r = ArchGDAL.readraster(s)
            #r = AG.readraster(x)
            r = ArchGDAL.getband(r,lyr)
            #nd = AG.getnodatavalue(r)
            dx = reverse(r, dims=2)
            dx = reverse(dx, dims=1)
            println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
            bitmat = dx .> msk
            # Filter the dx matrix using bitmat
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]

            dx = dx_output
            c_plot = Plots.contour(dx,
                    #c=:matter,
                        title=basename(s),xaxis="",yaxis="")
            step = step
            for i in 1:step:size(dx, 1)
                for j in 1:step:size(dx, 2)
                    value = round(dx[i, j]; digits=2)
                    color = isnan(value) ? :white : :black
                    annotate!(j, i, Plots.text(string(value), 7, color, :center, 
                        halign=:center, rotation=35.0))
                end
            end

            # for i in 1:step:size(dx, 1)
            #     for j in 1:step:size(dx, 2)
            #         value = round(dx[i, j]; digits=2)
            #         color = isnan(value) ? :white : :black
            #         rotation_char = isnan(value) ? "" : "↺"
            #         annotate!(j, i, text(string(value) * rotation_char, 5, color, :center))
            #     end
            # end

            # Show the plot
            display(c_plot)
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

        function agmask(s::AbstractString;step=50,lyr=1,lower=0,upper=1000)
            """
            AG only.
            plots heatmap from raster with ArchGDAL
            upper and lower bound
            """
            if !isfile(s)
                error("file not found!")
            end
            if !endswith(s,".nc")
                @warn("file seems not to be a NetCDF !")
            end
            r = ArchGDAL.readraster(s)
            dx = ArchGDAL.getband(r,lyr)
            #ArchGDAL.getnodatavalue(dx)
            
            if endswith(s,".nc")
                dx = reverse(dx, dims=2)
                dx = reverse(dx, dims=1)
            end
            #findmin(dx)
            if (minimum(dx) <= -9999.0) # || (maximum(dx) >= 9999.0)
                println("extrema: ",join(extrema(dx)," <-> "))
                dmin=minimum(dx)
                @info "
                $dmin set to NaN!
                lower bound: $lower
                upper bound: $upper
                "
                dx = replace(dx, minimum(dx)=>NaN)
                #replace!(dx, minimum(dx)=>missing) #err
                #replace!(dx, minimum(dx)=>NaN)
            end
            
            #heatmap(dx)
            #extrema(dx)
            
            #println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
            
            bitmat = (dx .> lower) .& (dx .< upper)
            # Filter the dx matrix using bitmat
            dx_filtered = dx .* bitmat
            dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
            dx_output .= NaN
            dx_output[bitmat] .= dx_filtered[bitmat]

            if allequal(dx_output)
                @error("all values are equal!")
                return
            end

            dx = dx_output
            heatmap_plot = heatmap(dx, c=:matter,
                            title=basename(s), 
                            xaxis="", yaxis="");
            step = step
            for i in 1:step:size(dx, 1)
                for j in 1:step:size(dx, 2)
                    value = round(dx[i, j]; digits=2)
                    color = isnan(value) ? :white : :black
                    annotate!(j, i, Plots.text(string(value), 7, color, :center, 
                        halign=:center, rotation=-35.0))
                end
            end
            # Show the plot
            display(heatmap_plot)
        end

        function all_values_equal(df::DataFrame)
            for col in eachcol(df)
                if any(col .!= col[1])
                    return false
                end
            end
            return true
        end

        function pall(files::Vector{DataFrame};toyr=false)
            """
            reduces + merges by date + plots all
            """
            #files = dfs
            bns = try
                broadcast(x->DataFrames.metadata(x)|>only|>last|>basename,files)
            catch
                @warn "No basenames in metadata!"
                raw""
                end

            try
                for i in 1:length(files)
                    s = Symbol.(filter(x->!occursin(r"year|date",x),names(files[i])))
                    nn = bns[i]|>x->split(x,".")|>first
                    for x in s
                        #println(string(x)*"-"*nn)
                        newname = string(x)*"-"*nn
                        #rename!(i,s[x]=>s[x]*"-"*bns[x])
                        rename!(files[i],x=>newname)
                    end
                end
            catch
                @warn "error in renaming!"
                @debug begin
                    nms=map(x->names(x),files)
                    "names are: $nms"
                end
                #@warn "error in renaming!"
            end

            if toyr==true
                dfs = broadcast(x->yrsum(x),files)
                df = reduce((left, right) -> 
                innerjoin(left, right, on = :year,
                makeunique=true
                ),dfs)
            else
                df = reduce((left, right) -> 
                innerjoin(left, right, on = :date,
                makeunique=true
                #renamecols = lowercase => uppercase
                ),files)
            end
            
            
            ti = try
                #DataFrames.metadata(df)|>only|>last|>basename
                z = map(x->split(x,".")|>first,bns)
                if length(z)>5
                    z = z[1:5]
                    "Merged Series"
                end
                #"Series of "*join(z," ")
                join(z," ")
            catch
                @warn "No basenames for title"
            ti = raw""
                end

            if length(ti)>30
                ti = ti[1:30]*"..."
            end

            if all_values_equal(df[!,Not(Cols(r"date|year|month|day"))])==true
                @error "all values are equal!"
                return
            end

            if (any(x->occursin("year",x),names(df)))
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
            #    df = df[!,Not(:date)]
                p = @df df Plots.plot(:year,cols(s),legend = :outerbottomright, title=ti)
            else    
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
                p = @df df Plots.plot(:date,cols(s),legend = :outerbottomright, title=ti)
            end
            return p
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

        function climateplot(temp::Regex,prec::Regex;col="tot_average")
            """
            col = subbasin of interest
            wa.climateplot(r"temp",r"pre";col="tot_average")
            ws. prc and temp tauschen un opacitiy einstellen....

            """
            col = Symbol(col)
            prec = waread(prec)
            yrs = year.(prec.date)|>unique|>length
            prec = monsum(prec)
            precvec = vec(Matrix(select(prec, col)))
            precvec = precvec ./ yrs
            
            temp = waread(temp)|>monmean
            tempvec = vec(Matrix(select(temp, col)))
            month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
            p1 = Plots.bar(prec.month, precvec, color=:cornflowerblue, 
                xlabel="", 
                #xlabel="Months", 
                xflip=false,
                ylabel="Precipitation [mm]", 
                legend=false, yflip=true);
            xticks!(1:12, month_abbr)
            for i in prec.month
                val = round(precvec[i]; digits=1)
                annotate!(i, precvec[i], text("$(val)",8, :bottom))
            end
            p2 = twinx();
            #ann2 = map(x->string.(round(x;sigdigits=0)),tempvec)
            plot!(p2, tempvec, xlabel="", 
                ylabel="Temperature [°C]", color=:coral2,
                #annotations = (temp.month,tempvec, ann2, :center),
                label=false, linewidth=3);
            # # Add annotations for temperature values
            # for i in 1:length(tempvec)
            #     val = round(tempvec[i]; sigdigits=1)
            #     annotate!(i, tempvec[i], text("$(val)",7, :center))
            # end
            return p1
        end

        function climateplot(temp::DataFrame,prec::DataFrame;col::AbstractString)
            """
            ws. prc and temp tauschen un opacitiy einstellen....
            ,col::String name of 
            twinx() dreht komplett alles.
            """
            #col = propertynames(temp)[col]
            #temp = t
            
            #col = first(Symbol.(filter(x->occursin(r"$col"i,x),names(temp))))
            #prec = pr
            yrs = year.(prec.date)|>unique|>length
            prec = monsum(prec)
            precvec = vec(Matrix(select(prec, col)))
            #precvec = vec(Matrix(select(prec, Not(:month))))
            precvec = precvec ./ yrs
            
            temp = (temp)|>monmean
            tempvec = vec(Matrix(select(temp, col)))
            #tempvec = vec(Matrix(select(temp, Not(:month))))
            month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
            p1 = Plots.plot(temp.month,tempvec, xlabel="", 
                ylabel="Temperature [°C]", color=:coral2,
                yflip = false,
                label=false, linewidth=3)
            # # Add annotations for temperature values
            for i in 1:length(tempvec)
                val = round(tempvec[i]; sigdigits=1)
                annotate!(i, tempvec[i], text("$(val)",7, :center))
            end
            
            #twinx()
            #prec.month, 
            Plots.bar!(precvec, color=:cornflowerblue, 
            xlabel="",     #xlabel="Months", 
            xflip=false,
            ylabel="Precipitation [mm]", 
            legend=false    #, yflip=true
            )
            xticks!(1:12, month_abbr)
            for i in prec.month
                val = round(precvec[i]; digits=1)
                annotate!(i, precvec[i], text("$(val)",8, :bottom))
            end
            return p1
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

        function pydf_to_julia(py_df::PyObject)
            """
            no transposing
            """
            # Convert each column of the Python DataFrame to a Julia array
            col_names = py_df.columns  # Get the column names from the Python DataFrame
            col_arrays = [convert(Array, py_df[col]) for col in col_names]
            # Create a Julia DataFrame using the converted arrays and column names
            julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, col_arrays))    
            return julia_df
        end 

        function pydf(py_df::PyObject)
            """
            Convert each column of the Python DataFrame to a Julia array
            """
            # fn = filename
            # pyo = py"""waread3($fn).reset_index(drop=False)"""
            # pdf = wa.pydf(pyo)
            # names(pdf)
            # dfp(pdf)
            col_names = py_df.columns  # Get the column names from the Python DataFrame
            col_names = convert(Array, col_names)
            col_arrays = [convert(Array, py_df[col]) for col in col_names]
            jdf = DataFrame(col_arrays, :auto)
            #size(jdf)
            fn = try
                py_df.filename
            catch
                @info "no filename present"
            end

            metadata!(jdf, "filename", fn, style=:note);
            rename!(jdf, col_names);
            return jdf
        end

        function tovec(x::DataFrame, col::Any)
            df = select(x,col)
            println(names(df))
            return vec(Matrix(df))
        end

        function zp(func::Any)
            """
            prints out the function definition
            """
            #pt = joinpath(@__DIR__,"func-win.jl")
            pt = src_path*"/func-win.jl"
            _str = "$(func)"
            #_str = r"^\s*"*"$(func)"
            #_str = r"(?:^\s*)"*"$(func)"
            #r"\n\s*\n" ##empty line regex.
            #readbetween(open(pt),string(_str), "end")
            #readbetween(open(pt),Regex(_str),r"end$")
            readbetween(open(pt),Regex(_str),r"^\s*function")
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

        function read_soildata_4(filename::String)
            """
            soil types
            returns a Vector of Strings
            """
            file = open(filename) do io
                a = readbetween(io, "soil_table", "special_output")
                return(a)
            end
            file = file[2:end-1]
            #data = Grep.grep(r"^[0-9].* {",file)
            data = Grep.grep(r"^[0-9]|^[ ;]|^[;]",file)
            filter!(x->length(x)>3,data)
            return data
        end

        function pyread_meteo(s::AbstractString;hdr=0)
            pd = pyimport("pandas")
            ddf = pd.read_csv(s, delim_whitespace=true, 
                header=hdr,
                na_values=-9999,
                low_memory=false,
                verbose=true)
            #ddf.filename=basename(s)
            ddf = wa.pydf(ddf)
            
            #dt = Date.(map(x -> join(x, "-"), eachrow(ddf[:, 1:3])))
            dt = select(ddf,1:3)
            Z = []
            for i in eachcol(dt)
                col = tryparse.(Int, i)
                push!(Z,col)
            end
            dt = DataFrame(Z,:auto)
            
            msk = tryparse.(Int, ddf.YY)
            ddf = ddf[findall(!isnothing,msk),:]
            dt = dt[findall(!isnothing, dt.x1),:]
            dt = Date.(map(k -> join(k, "-"), eachrow(dt[:, 1:3])))
            nd = hcat(ddf[!, 5:end], dt)
        
            for col in names(nd)[(eltype.(eachcol(nd)) .<: String)]
                nd[!, col] .= tryparse.(Float64, col)
            end
            rename!(nd, ncol(nd) =>"date")
            metadata!(nd, "filename", s, style=:note);
            return nd
        end

        function pyread_old(s::AbstractString;hdr=0)
            pd = pyimport("pandas")
            ddf = pd.read_csv(s, delim_whitespace=true, 
                header=hdr,
                na_values=-9999,
                low_memory=false,
                verbose=true)
            #ddf.filename=basename(s)
            ddf = wa.pydf(ddf)
            
            dt = Date.(map(x -> join(x, "-"), eachrow(ddf[:, 1:3])))
            
            nd = hcat(ddf[!, 5:end], dt)

            for col in names(nd)[(eltype.(eachcol(nd)) .<: String)]
                nd[!, col] .= tryparse.(Float64, col)
            end
            rename!(nd, ncol(nd) =>"date")
            metadata!(nd, "filename", s, style=:note);
            return nd
        end

        function waread3_py(x, flag=true)
            # use py"""...""" syntax to access the Python function
            py"""
            import pandas as pd
            import datetime

            def waread3(x, flag=True):
                if flag:
                    df = pd.read_csv(x, delim_whitespace=True, header=0,
                                    na_values=-9999, verbose=True,engine='c')
                    if 'YY' not in df.columns:
                        print("Column 'YY' not found in the CSV file.")
                        return None
                    if df.iloc[:, 0].dtype == 'object':
                        print('First column contains strings, subsetting to Int...')
                        df = df[~df.iloc[:, 0].str.contains("[A-z]|-", na=False)]
                    source_col_loc = df.columns.get_loc('YY')        
                    df['date'] = df.iloc[:, source_col_loc:source_col_loc +
                                        3].apply(lambda x: "-".join(x.astype(str)), axis=1)
                    df = df.iloc[:, 4:]
                    df['date'] = pd.to_datetime(df['date'])
                    #df.set_index('date', inplace=True)
                    df.iloc[:,0:-2] = df.iloc[:,0:-2].apply(lambda x: x.astype(float), axis=1)
                    df.filename = x
                    print(df.filename,"done!")
                    return df
                else:
                    print('Date parse failed, try reading without date transformation...')
                    df = pd.read_csv(x, delim_whitespace=True, comment="Y", skip_blank_lines=True).dropna()
                    df.filename = x
                    print(df.filename,"done!")
                    return df
            """
            # call the Python function with the arguments and return the result
            return py"waread3($x, $flag)"
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

        function nqp(a::Regex,b::Regex;)
            """
            qplot from regex
            """
            a = waread(a)
            b = waread(b)
            col = ncol(a)-1

            a = a[!,Cols(col,:date)]
            b = b[!,Cols(col,:date)]  

            df = mall(a,b)
            #qplot(x::Vector{Float64},y::Vector{Float64})
            return qplot(df)
        end

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

        function wqplot(file_path::AbstractString)
            data = read_wq(file_path)
            p = @df data plot(data[!,1],cols(propertynames(data)[2:end]))
            println(describe(data))
            return p
        end

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
                        !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly|eps)$", filename)
                        )
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
            # results = filter(x -> isfile(x) && 
            # !occursin(r"yr|mon|grid|scn"i,x) && 
            # !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly|eps)$", x), results)
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
            
        #grep_with_context("routing_model", file_paths[1], 2)

        function grep_files(pattern, file_paths, context)
            """
            pts = readdir(dirname(infile);join=true)
            filter!(x->occursin(r"ctl",x),pts)
            grep_files("routing_model", pts, 2)
            """
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

        function stplot(fn::String)
            """
            reads from stationdata, reprojects and plots
            """
            fl = CSV.read(fn,DataFrame;limit=4)
            xc = fl[2,5:end]|>collect
            yc = fl[3,5:end]|>collect
            pts = ArchGDAL.IGeometry[]
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
            od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
            p = plot(od.geometry);
            for (i, pt) in enumerate(od.geometry)
                #x = od.xc[i]
                x = ArchGDAL.getx(od.geometry[i], 0)
                #y = od.yc[i]
                y = ArchGDAL.gety(od.geometry[i], 0)
                name = od.name[i]
                annotate!(x, y, text(name, 8, :black, :bottom, :left))
            end
            return plot!(p)
        end

        function runwasim(ctlfile;exepath="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05/wasimvzo64.exe")
            """
            runs inside REPL
            """
            try
                exec = normpath(exepath)
                run(`$exec $ctlfile`)
            catch e
                println("no valid input!")
            end
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

        function corrbar(a::Regex, b::Regex)
            # Load the DataFrames
            df_A = (a |> dfonly |> first |> readdf)
            df_B = (b |> dfonly |> first |> readdf)

            nma = DataFrames.metadata(df_A)|>only|>last|>basename
            nmb = DataFrames.metadata(df_B)|>only|>last|>basename

            # Remove last column from df_B
            select!(df_B, Not(ncol(df_B)))
            select!(df_A, propertynames(df_B)) 

            # Calculate correlations and replace NaN with 0
            #correlations = cor.(eachcol(df_A_subset), eachcol(df_B)).^2
            correlations = Vector{Float64}(undef, size(df_A, 2))
            for i in 1:size(df_A, 2)
                correlations[i] = cor(df_A[:, i], df_B[:, i])^2
            end
            replace!(correlations, NaN => 0)

            p0 = bar(1:size(df_A, 2), correlations,
                legend = false,
                title = "$nma vs $nmb",
                fillcolor = ifelse.(correlations .> 0.35, 
                    "cornflowerblue", "coral2"),
                xticks = (1:size(df_A, 2), propertynames(df_A)),
                xrotation = 45,
                xlabel = "",
                ylabel = "Correlation R²",
                left_margin = 10mm,
                bottom_margin = 2mm);

            ann = map(x->string.(round(x;sigdigits=2)),correlations)

            for i in 1:size(df_A, 2)
                Plots.annotate!(i,correlations[i],
                (ann[i],9,:center,:top,:black))
                println("R² "*ann[i]*" of Basin "*names(df_A)[i]*" added")
            end

        return p0
        end

        function corrbar(a::DataFrame, b::DataFrame)
            # Load the DataFrames
            df_A = a
            df_B = b
            ti = try 
                nma = DataFrames.metadata(df_A)|>only|>last|>basename
                nmb = DataFrames.metadata(df_B)|>only|>last|>basename
                ti = "$nma vs $nmb"
                
            catch
                @info    "no metadata in $a or $b !"
                ti = "$a vs $b"
                
            end

            # Remove last column from df_B
            select!(df_B, Not(ncol(df_B)))
            select!(df_A, propertynames(df_B)) 

            # Calculate correlations and replace NaN with 0
            #correlations = cor.(eachcol(df_A_subset), eachcol(df_B)).^2
            correlations = Vector{Float64}(undef, size(df_A, 2))
            for i in 1:size(df_A, 2)
                correlations[i] = cor(df_A[:, i], df_B[:, i])^2
            end
            replace!(correlations, NaN => 0)

            p0 = bar(1:size(df_A, 2), correlations,
                legend = false,
                title = ti,
                fillcolor = ifelse.(correlations .> 0.35, 
                    "cornflowerblue", "coral2"),
                xticks = (1:size(df_A, 2), propertynames(df_A)),
                xrotation = 45,
                xlabel = "",
                ylabel = "Correlation R²",
                left_margin = 10mm,
                bottom_margin = 2mm);

            ann = map(x->string.(round(x;sigdigits=2)),correlations)

            for i in 1:size(df_A, 2)
                Plots.annotate!(i,correlations[i],
                (ann[i],9,:center,:top,:black))
                println("R² "*ann[i]*" of Basin "*names(df_A)[i]*" added")
            end

        return p0
        end

        function nse2(df::DataFrame)
            """
            takes cols from CSV.read
            """
            observed, simulated = df[:,6],df[:,5]
            return (1 - (sum((simulated .- observed).^2) / sum((observed .- mean(observed)).^2)))
        end

        function nse_rec()
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
                        !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly|eps)$", filename)
                        )
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
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
                push!(out,Dict("name"=>x,"NSE"=>nse2(dd)))
                catch
                    @warn("$x can't be loaded as a DataFrame, skipping ...")
                    continue
                end
                
                
            end
            df = DataFrame(out)
            df.NSE .= replace(df.NSE, NaN => missing)
            dropmissing!(df)
            return sort(df,:NSE;rev = true)
        end

        function nctodf(str::String)
            """
            use of Rasters.lookup
            """
            #x = @gl "tsoil"
            rr = Raster(str)
            # dims(rr)|>collect
            # rr.dims|>collect
            # rr[Dim{:t}(Rasters.Between(2,end))] |>contourf
            
            mt = split(string.(rr[:,:,1].refdims)[1],",")|>first
            #occursin("Dim{:t}",)

            #if last(dims(rr))=="Dim{:t}"
            if mt=="Dim{:t"
                ti = Rasters.lookup(rr, Dim{:t})
            else
                ti = try
                    Rasters.lookup(rr, Ti)
                catch
                    @error "no time dimension found!"
                    return
                end 
            end
            
            #map(x->mean(x),rr.dims)
            #dims(rr, (Dim{:t})) isa Vector{Float64}
            
            
            #ag = Rasters.aggregate(Rasters.Center(), rr, (Y(20), X(20));)
            #plot(ag)
            #x,y,z = map(x->round(x ./2;digits=0),size(rr))
            x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
            #x,y,z = map(x->parse.(Int,x ./2),size(rr))
            #length(rr[1,:,:])
            
            df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
            df = hcat(df, parent(ti),makeunique=true)
            #rename!(df,1=>Rasters._maybename(rr),2=>"date")
            rename!(df,1=>Rasters._maybename(rr),2=>"layer")

            DataFrames.metadata!(df, "filename", str, style=:note);        
            return df

        end

        function vef2(df::DataFrame)
            """
            observed, simulated = df[:,6],df[:,5]
            """
            obs, sim = df[:,6],df[:,5]
            return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
        end

        function rec()
            """
            reads recursively from qout and calulates KGE and NSE
            sorted by NSE
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
                        !occursin(r"\.(log|png|svg|txt|html|ftz|ftz_0|list|nc|zip|7z|xml|sh|grd|yrly|eps)$", filename)
                        )
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
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
                push!(out,Dict("name"=>x,"NSE"=>nse2(dd),"KGE"=>kge2(dd),"VE"=>vef2(dd)))
                catch
                    @warn("$x can't be loaded as a DataFrame, skipping ...")
                    continue
                end
                
            end
            df = DataFrame(out)
            df.NSE .= replace(df.NSE, NaN => missing)
            df.KGE .= replace(df.KGE, NaN => missing)
            df.VE .= replace(df.VE, NaN => missing)
            dropmissing!(df)
            df.bn .= map(x->splitdir(x)|>last,df.name)
            df.pt .= map(x->splitdir(dirname(x))|>last,df.name)
            return sort(df,:NSE;rev = true)
        end

        function readall(path::Regex)
            """
            tries to read all files in a directory as df
            """
            v::Vector{String} = readdir();
            #v = v[broadcast(x->!endswith(x,"nc"),v)];
            files = v[(broadcast(x->occursin(path,x),v))];
            dfs::Vector{DataFrame} = []
            for file in files
                if isfile(file) && (!occursin(r"xml$|fzt$|ftz$|log$|ini$|^wq|yrly$|nc$|png$|svg$",file))
                    file_path = file
                    println("reading ",file_path,"...")
                    #p1 = waread(file_path)
                    ##modified waread
                    ms = ["-9999","lin","log","--"]
                    df = CSV.read(file_path, DataFrame; 
                        delim="\t", header=1, missingstring=ms, 
                        #maxwarnings = 1, 
                        silencewarnings = true,
                        rows_to_check = 100, ignorerepeated = true,
                        normalizenames=true, types=Float64)
                    df = dropmissing(df, 1)
                    dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
                    df.date = dt2
                    df = select(df, Not(1:4))
                    DataFrames.metadata!(df, "filename", file_path, style=:note)
                    for x in names(df)
                        if startswith(x,"_")
                            newname=replace(x,"_"=>"C", count=1)
                            rename!(df,Dict(x=>newname))
                        end
                    end
                    
                    push!(dfs, df)
                end
            end
            return(dfs)
        end

        function readall(path::Vector{Any})
            """
            tries to read all files in a directory as df
            """
            files = dfonly(path)
            dfs::Vector{DataFrame} = []
            for file in files
                if isfile(file) && (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
                    file_path = file
                    println("reading ",file_path,"...")
                    #p1 = waread(file_path)
                    ##modified waread
                    ms = ["-9999","lin","log","--"]
                    df = CSV.read(file_path, DataFrame; 
                        delim="\t", header=1, missingstring=ms, 
                        #maxwarnings = 1, 
                        silencewarnings = true,
                        rows_to_check = 100, ignorerepeated = true,
                        normalizenames=true, types=Float64)
                    df = dropmissing(df, 1)
                    dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
                    df.date = dt2
                    df = select(df, Not(1:4))
                    DataFrames.metadata!(df, "filename", file_path, style=:note)
                    for x in names(df)
                        if startswith(x,"_")
                            newname=replace(x,"_"=>"C", count=1)
                            rename!(df,Dict(x=>newname))
                        end
                    end
                    
                    push!(dfs, df)
                end
            end
            return(dfs)
        end

        function readbetween(io::IO, start::Regex, stop::Regex)
            output = Vector{String}()
            while !eof(io)
                line = readline(io)
                if occursin(start,line)
                    push!(output, line)
                    while !eof(io)
                        line = readline(io)
                        if occursin(stop,line)
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

        function vef(sim, obs)
            return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
        end

        function vef(x::DataFrame)
            """
            observed, simulated = df[:,6],df[:,5]
            """
            xd = x[:,Not(Cols(r"date|year"))]
            sim = DataFrames.select(xd, map(z->occursin(r"C*[0-9]",z),
                names(xd)))|>dropmissing|>Matrix|>vec
            obs = DataFrames.select(xd, map(z->!occursin(r"C*[0-9]",z),
            names(xd)))|>dropmissing|>Matrix|>vec
            
            return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
        end

        function vef(x::String)
            """
            observed, simulated = df[:,6],df[:,5]
            """
            ms=["-9999","lin","log"]
            df::DataFrame = CSV.read(x,DataFrame,
            missingstring=ms,
            types = Float64,
            delim="\t",
            silencewarnings=true,
            normalizenames=true) |> dropmissing
            #drop=(i, nm) -> i == 4) |> dropmissing
            
            obs, sim = df[:,6],df[:,5]
            
            return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
        end

        function seldf(str::String,dfs::Vector{DataFrame})
            """
            selects columns from vector of dfs ...
            """
            str  = join([str,"date"],"|")
            #filter(x->select(occursin(str,x)),dfs)
            return map(k -> k[!,Cols(Regex(str))],dfs)
            
        end

        function monp(lk::String;doleg=true)
            df = waread(lk)
            df = monmean(df)
            #str = [ @sprintf("%02i", x) for x in (df.month) ];
            #mn = [ monthname(x) for x in (df.month) ]
            mn = [ monthabbr(x) for x in (df.month) ]

            #linestyle dashdtriangle is unsupported with Plots.GRBackend(). 
            #Choose from: [:auto, :dash, :dashdot, :dashdotdot, :dot, :solid]
            
            # if ncol(df)>length(linestyles)
            #     linestyles = vcat(linestyles,linestyles)
            # end

            #linestyles = :auto

            ln = Symbol.(filter(x->!occursin(r"date|month|year",x),names(df)))
            # @df df StatsPlots.plot(str,cols(ln),
            #     fillalpha=0.75, linewidth=1.25, legend=true)   
            Plots.plot()

            if doleg
                for (i, col) in enumerate(ln)
                    Plots.plot!(mn, df[!,col], fillalpha=0.75, 
                    #linestyle=linestyles[i], 
                    linestyle=:auto, 
                    linewidth=1.25, 
                    #label=[string(col)],
                    label=string(col),
                    legend=true)
                end

                return Plots.plot!() 
            end 

            for (i, col) in enumerate(ln)
                Plots.plot!(mn, df[!,col], fillalpha=0.75, 
                linestyle=:auto, 
                linewidth=1.25, 
                legend=false)
            end
            return Plots.plot!() 
        end

        function monp(lk::DataFrame;doleg=true)
            
            df = monmean(lk)
            #str = [ @sprintf("%02i", x) for x in (df.month) ];
            #mn = [ monthname(x) for x in (df.month) ]
            mn = [ monthabbr(x) for x in (df.month) ]
            ln = Symbol.(filter(x->!occursin(r"date|month|year",x),names(df)))
            Plots.plot()

            if doleg
                for (i, col) in enumerate(ln)
                    Plots.plot!(mn, df[!,col], fillalpha=0.75, 
                    #linestyle=linestyles[i], 
                    linestyle=:auto, 
                    linewidth=1.25, 
                    #label=[string(col)],
                    label=string(col),
                    legend=true)
                end

                return Plots.plot!() 
            end 

            for (i, col) in enumerate(ln)
                Plots.plot!(mn, df[!,col], fillalpha=0.75, 
                linestyle=:auto, 
                linewidth=1.25, 
                legend=false)
            end
            return Plots.plot!() 
        end

        function ezplot(dirpath::String)
            #cd(dirpath)
            
            files = readdir(dirpath;join=true)
            rs = try
                filter(file -> occursin(Regex(".ezg"),file),files)|>first
                
            catch
                @error "No .ezg file found in $dirpath !"
                return nothing  # 
            end 
            
            lstr = filter(file -> occursin(Regex(".lnk"),file),files)|>first
            ezg = Raster(rs;missingval=-9999)   |>Rasters.trim
            lnk = Raster(lstr;missingval=-9999)|>Rasters.trim
            p = Plots.plot(ezg,c=:lightrainbow,
                xaxis="",yaxis="");
            Plots.plot!(lnk,c=:white,
                title=dirname(dirpath),
                xaxis="",yaxis="")

            return p
            
            # r = ArchGDAL.readraster(rs)
            # dx = ArchGDAL.getband(r,1)
            # dp = AG.polygonize(dx)
        end

        function looks_like_number(str::AbstractString)
            try
                parse(Float64, str)
                return true
            catch
                return false
            end
        end

        function pyread(x)
            """
            pyreader, reads all as stings, conversion later.
            """
            pd = pyimport("pandas")
            df = pd.read_table(x, 
                engine="c",
                verbose=true,
                low_memory=false,
                header=0,
                skipinitialspace=true,
                dtype="str",              #new!
                na_values=[-9999]
                )
            col_names = df.columns  # Get the column names from the Python DataFrame
            col_names = convert(Array, col_names)
            col_arrays = [convert(Array, df[col]) for col in col_names]
            filtered_rows = broadcast(x->looks_like_number(x),col_arrays[1])
            
            df = DataFrame(col_arrays, :auto)
            
            df = df[filtered_rows, :]
            
            rename!(df, col_names)
            if "YY" ∉ names(df)
                println("Column 'YY' not found in the CSV file.")
                return nothing
            end

            function parse_date(row)
                try
                    year = parse(Int, row[1])
                    month = parse(Int, row[2])
                    day = parse(Int, row[3])
                    return Date(year, month, day)
                catch
                    return missing  # Return missing if parsing fails
                end
            end
            
            df.date = map(parse_date, eachrow(df))
            
            df = select(df, Not(1:4))
            
            for x in names(df)
                if looks_like_number(x)
                    newname=replace(x,"$x"=>"C$x", count=1)
                    rename!(df,Dict(x=>newname))
                end
            end

            #map(y->typeof(y),eachcol(df))

            # Iterate over column names
            for colname in names(df)
                # Check if the column type is not Date
                if eltype(df[!, colname]) != Date
                    df[!, colname] .= tryparse.(Float64, df[!, colname])
                end
            end

            DataFrames.metadata!(df, "filename", x, style=:note)
            return df 
        end

        function addplot(specfile::String)
            """
            adds points to an existing plot from a specfile
            """
            fl = CSV.read(specfile,DataFrame;limit=4)
            xc = fl[2,5:end]|>collect
            yc = fl[3,5:end]|>collect
            pts = ArchGDAL.IGeometry[]
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                #pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
            od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
            Plots.plot!(od.geometry;c=:black,shape = :x) #:utriangle)
            for (i,pt) in enumerate(od.geometry)
                x = ArchGDAL.getx(od.geometry[i], 0)
                y = ArchGDAL.gety(od.geometry[i], 0)
                name = od.name[i]
                #Plots.text(name,:consolas, 8, :black, :bottom, :left)
                #Plots.annotate!(x, y, Plots.text(name, 8, :black, :bottom, :left))
                Plots.annotate!(x, y, Plots.text(name, 8, :black, :bottom, :left))
            end
            Plots.plot!()
        end

        function qplot(x::Regex, y::Regex)
            
            """
            takes first two cols of df and plots r2 QQ
            """
            df1,df2 = map(z->wa.waread(z),[x,y])
            # df1 = waread(x)
            # df2 = waread(y)
            
            # if any(map(x->contains(x,"date"),names(df2)))
            #     df2 = df2[!,Not(:date)]
            # end
            @info "takes first and last column, merges, and plots QQ"
            df = mall(
                df1[!,Cols(1,ncol(df1))],
                df2[!,Cols(1,ncol(df2))])
            
            if any(map(x->contains(x,"date"),names(df)))
                df = df[!,Not(:date)]
            end       
            if ncol(df)>2
                df = df[!,1:2]
            end
            r2 = round(cor(df[!,1], df[!,2])^2, digits=3)
            p = df |> x -> qqplot(x[!,1],x[!,2], 
                #title = "R² = "*string(ti),
                qqline = :fit)
                #color = :grays) # erstellt ein QQ-Diagramm <- black
            xlabel!(p,names(df)[1])
            ylabel!(p,names(df)[2])

            #ve = round(vef(df), digits=3)
            ve = df |> x -> vef(x[!,1],x[!,2])
            ve = round(ve, digits=3)
            Plots.annotate!(p,:topleft, 
            Plots.text(
                "R² = "*string(r2)*"\nVE = "*string(ve), 
            11,:black, halign=:left)
            )
            #
            
            #annotate!(j, i, Plots.text(string(value), 7, color, :center, 
            #halign=:center, rotation=-35.0))
        end

        function jsrast(R::String)
            rpath="/mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions"
            run(`wsl -d Ubuntu-18.04 -e Rscript "$rpath/plotly_rast.R" $R`)
        end

        function oplat()
            files = filter(isfile, readdir())
            sorted_files = sort(files, by = mtime, rev = true)
            if !isempty(sorted_files)
                lat = sorted_files[1]
                println("opening: $lat ...")
                run(`cmd /c start $lat`)
            end
        end   

        function rmlat()
            files = filter(isfile, readdir(;sort=false))
            sorted_files = sort(files, by = mtime, rev = true)
            if !isempty(sorted_files)
                lat = sorted_files[1]
                println("This deletes the latest created file, i.e: ", lat)
                print("continue? (y/n): ")
                reply = readline()
                println(reply)
                if lowercase(string.(reply)) == "y"
                    rm(lat, force=true)
                    println("Deleted: ", lat)
                else
                    @error "abort...."
                end
            end
        end    

        function lat()
            files = filter(isfile, readdir(;sort=false))
            sf = sort(files, by = mtime, rev = true)
            lat = sf[1]
            if length(sf) < 6
                sf = join(sf,"\n")
            else
                sf = join(sf[1:6],"\n")
            end
            printstyled("$sf\n--> ",color=:green)
            return lat
        end

        function toMain()
            fnames = names(Main.WaSiM, all=true)
            for submodule in fnames
                @eval import Main.WaSiM.$submodule
            end
        end

        function rfread(x::String)
            dt = rimport("data.table")
            #rdf = dt.fread(x,nThread=8,skip="Y",colClasses="numeric")        #verbose=true,
            rdf = dt.fread(x,nThread=8,skip="Y",check=true,strip=true)
                #na="-9999")
            dt.setnames(rdf,1,"YY")
            df = rcopy(rdf)
            filter!(x -> x.YY != "YY",df)
            for i in 5:ncol(df)
                df[!,i] .= replace(df[!,i],-9999.0 => missing)
            end
            for i in 1:4
                if isa(df[!,i],Vector{String})
                    df[!,i] .= tryparse.(Int,df[!,i])
                end
            end
            # dropmissing!(df)
            dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
            df.date = dt2
            df = select(df, Not(1:4))
            DataFrames.metadata!(df, "filename", x, style=:note)
            for x in names(df)
                if startswith(x,"X")
                    newname=replace(x,"X"=>"C", count=1)
                    rename!(df,Dict(x=>newname))
                end
            end
        
            return df     
        end

        function ncdf(str::String)
            """
            use of Rasters.lookup
            recode
            """
            rr = Raster(str)
            @info "this extracts values from the midpoint of the raster.
            rename!(df,1=>Rasters._maybename(rr),2=>...)
            last dimension is time..."

            #nm = Rasters.name(rr)
            nm = Rasters._maybename(rr)

            if dims(rr)[1]|>name != :X
                new_dims = (X(parent(Rasters.lookup(rr,1))), Y(parent(Rasters.lookup(rr,2))), Ti(parent(Rasters.lookup(rr,3))))
                # Create a new raster with the new dimensions and the same data
                rr = try
                    Raster(rr.data, dims=new_dims)
                catch
                    @error "var conversion failed!"
                    return
                end 
            end

            ti = Rasters.lookup(rr,3)
            x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
            df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
            df = hcat(df, parent(ti), makeunique=true)
            #rename!(df,1=>Rasters._maybename(rr),2=>"date")
            rename!(df,1=>Symbol(nm),2=>:date)
            DataFrames.metadata!(df, "filename", str, style=:note);        
            return df
        end

        function ncdf(rr::Raster)
            """
            use of Rasters.lookup
            recode
            """
            @info "this extracts values from the midpoint of the raster.
            rename!(df,1=>Rasters._maybename(rr),2=>...)
            last dimension is time..."

            #nm = Rasters.name(rr)
            nm = Rasters._maybename(rr)

            if dims(rr)[1]|>name != :X
                new_dims = (X(parent(Rasters.lookup(rr,1))), Y(parent(Rasters.lookup(rr,2))), Ti(parent(Rasters.lookup(rr,3))))
                # Create a new raster with the new dimensions and the same data
                rr = try
                    Raster(rr.data, dims=new_dims)
                catch
                    @error "var conversion failed!"
                    return
                end 
            end

            ti = Rasters.lookup(rr,3)
            x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
            df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
            df = hcat(df, parent(ti), makeunique=true)
            #rename!(df,1=>Rasters._maybename(rr),2=>"date")
            rename!(df,1=>Symbol(nm),2=>:date)
            DataFrames.metadata!(df, "filename", nm, style=:note);
            
            #if df.date[1]|>typeof == CFTime.DateTimeNoLeap
            #parse(DateTime,string(df.date[1]))
            if df.date[1] != Date
                df.date = map(x->tryparse(DateTime,string(x)),df.date)
            end
            

            return df
        end

        function findlog(;lb=.4)
            """
            find LOG. R-SQUARE > .4 recursivley
            """
            v = qall(;recursive=true)
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

        function skipyr(df::DataFrame;fun=nothing)
            """
            skips first year and returns df
            optionally apply fun to df
            """
            fst = year.(df.date)[1]
            data = filter(:date => x -> Dates.year(x) > fst, df)
            #return fun ? fun(data) : data #if boolean
            if isnothing(fun)
                return data
            else
                return fun(data)
            end
        end

        function surf(x::Raster;msk=0.001,c=:jet,cam=(-20, 40))
            """
            3D plot with Rasters
            wa.surf(rr[Ti=12];msk=3,cam=(12,45),c=:lightrainbow)
            
            """
            #0.001 isa Float64
            #ti=Rasters._maybename(x)
            #ti=Rasters.name(x)

            if msk !== Float64
                msk = tryparse(Float64,string(msk))
            end       
            
            #xs = (length(x.dims) > 2) ? x[t=Int(x.dims[3][end])+1]  : x
            xs = x
            zm = (xs .>= msk)
                
            M = Rasters.mask(xs; with=zm);
            
            Plots.surface(
                M,
                c = c,
                #title = ti,
                legend= false,
                xlabel=" ",
                ylabel=" ",
                #zlabel=" ",
                camera = cam        
                )
            #     camera = (-20, 75))
        end

        
        function bardf(x::DataFrame;leg=:topright)
            "with DataFrame input"
                df = copy(x)
                y = filter(x->!occursin(r"date|year|month",x),names(df))
                s = map(y -> Symbol(y),y)
                    ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
                #df[!, :year] = year.(df[!,:date]);
                #df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
            nm = filter(x->occursin(r"date|year|month",x),
                names(df))|>first
            
            @df df Plots.plot(select(df,Symbol(nm))|>Matrix,
                    cols(s),
                    legend = leg, 
                    title = ti,
                    seriestype=:bar)
        end

        function ncdf(xs::RasterSeries)
            """
            use of Rasters.lookup
            @info "this extracts values from the midpoint of the raster.
            rename!(df,1=>Rasters._maybename(rr),2=>...)
            last dimension is time..."
            """
            #for rr in xs;println(name(rr));end
            dfs = []
            for rr in xs;
                @info "generating" (name(rr))
                nm = Rasters._maybename(rr)
                if dims(rr)[1]|>name != :X
                    new_dims = (X(parent(Rasters.lookup(rr,1))), Y(parent(Rasters.lookup(rr,2))), Ti(parent(Rasters.lookup(rr,3))))
                    # Create a new raster with the new dimensions and the same data
                    rr = try
                        Raster(rr.data, dims=new_dims)
                    catch
                        @error "var conversion failed!"
                        return
                    end 
                end
            
                ti = Rasters.lookup(rr,3)
                x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
                df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
                df = hcat(df, parent(ti), makeunique=true)
                #rename!(df,1=>Rasters._maybename(rr),2=>"date")
                rename!(df,1=>Symbol(nm),2=>:date)
                #DataFrames.metadata!(df, "filename", str, style=:note); 
                push!(dfs,df)
            
                end
                
            odf = reduce((left, right) -> 
                innerjoin(left, right, on = :date,makeunique=true), 
                dfs)

            return hcat(odf[!,Not(Cols(r"date"))],odf[:,Cols(r"date")])
        end

        function colsums(df::DataFrame)
            colsum = sum.(eachcol(df[!,Not(Cols(r"date|month|year"))]))
            return colsum
        end

        function subsum(df::DataFrame)
            """
            Subset a DataFrame to exclude only columns with zero sums.
            """
            column_sums = sum.(eachcol(df[!, Not(:date)]))
            # Find the column indices with sums equal to 0.0
            zero_sum_columns = findall(==(0.0), column_sums)
            # Subset the DataFrame to include only columns with zero sums
            dout = hcat(df.date, select(df, Not(zero_sum_columns)))
            rename!(dout, 1=>"date")
            return dout
        end


        function fsz(;rec=false)
            """
            Returns the total size of all files in the current directory non recursively.
            """
            total_size = 0
            files = readdir()  # Get a list of files in the current directory
            
            for file in files
                filepath = joinpath(pwd(), file)  # Get the full path of the file
                if isfile(filepath)
                    size = stat(filepath).size  # Get the file size
                    total_size += size
                end
            end
            
            nr=length(files)
        
            total_size_mb = total_size / (1024 * 1024)  # Convert size to megabytes
            total_size_mb = round(total_size_mb, digits=2)  # Round to 2 decimal places
            printstyled("Total size of $nr files in $(pwd()): $total_size_mb MB\n", color=:green)
            
            if rec
                dirs = readdir()
                for dir in dirs
                    if isdir(dir)
                        cd(dir)
                        fsz()
                        cd("..")
                    end
                end
            end
        end
        
        function fsize()
            """
            ALL names and size in MB via readdir
            """
            files = filter(x -> isfile(x), readdir())
            fzs = [(file, filesize(file) / 2^20) for file in files]
            tot = round(sum(map(x->x[2],fzs));digits=3)
            #printstyled("Total Size: $tot MB\n",color=:green)
            #return(DataFrame(Dict(fzs)))
            fzs = sort(fzs, by = x -> x[2], rev = true)
            #fzs = Dict(fzs)
            odf = rename(DataFrame(fzs),["file","size"])
            DataFrames.metadata!(odf, "Total Size", tot, style=:note)
            #DataFrames.metadata(odf)
            printstyled("Total Size: $tot MB\n",color=:green)
            return(odf)
        end

        function dfm(x::Union{Regex, String, DataFrame}; 
            ann = true, 
            log = false, 
            title = true,
            leg = false,  
            fun = wa.monmean, 
            mode=:line)
            """
            Plots a DataFrame and applies a function.
            
            Parameters:
            - x: Union{Regex, String, DataFrame} - Input data, can be a DataFrame or file path.
            - leg: Symbol - Legend position, default is :topright.
            - fun: Function - The function to apply to the DataFrame, default is wa.monmean.
            - mode: Symbol - Plot mode, can be :bar, :scatter, :line, :steppre, :steppost, :hist, :box; default is :line.
            - log: Bool - Logarithmic y-axis, default is false.
            - title: Bool - Title of the plot, default is true.
            
            Example Usage:
            dfm(s; fun=yrsum, mode=:scatter, leg=false)
            dfm(s; fun=wa.monmean, mode=:box, leg=false)
            """
            if isa(x, DataFrame)
                df = x
            else
                df = waread(x)
            end
            
            ti = try
                DataFrames.metadata(df)|>only|>last|>basename
            catch
            @warn "No basename in metadata!"
                ti = raw""
            end

            ln = Symbol.(filter(x -> !occursin(r"date|month|year", x), names(df)))
            
            if fun == false
                # Reorder the DataFrame
                dx = hcat(df[:, Cols(r"date")], df[!, Not(Cols(r"date"))])
                @info "No function applied!"
            else
                @info "Applying $fun ..."
                dx = fun(df)
            end

            if (ann == false || names(dx)[1] == "date")
                @info "No annotations to plot..."
                p1 = @df dx Plots.plot(
                    Vector(dx[!, 1]),
                    cols(ln),
                    legend = leg,
                    title = ti,
                    seriestype = mode)
                #return p1
            else
                # Annotate column names along the lines
                anns = map(x->string(x),ln)
                #anns, fontsize, rotation, halign, color
                #xans = map(x->Plots.text(x, 8, -20.0, :top, :red), anns) #bottom
                fnt = Plots.font(
                    family="sans-serif",
                    pointsize=8, 
                    #halign=:hcenter, 
                    valign=:bottom, 
                    #valign=:top, #:vcenter, 
                    rotation= -10.0, 
                    color=:red)
                
                xans = map(x->Plots.text(x, fnt), anns)
                # y_offset = map(x->median(dx[!,x]),ln)
                # y_offset = map(x->minimum(dx[!,x]),ln)
                #y_offset = map(x->minimum(dx[!,x]),ln) .* map(x->maximum(dx[!,x]),ln)
                #y_offset = map(x->maximum(dx[!,x]),ln) .- map(x->mean(dx[!,x]),ln)
                y_offset = map(x->minimum(dx[!,x]),ln) .+ map(x->mean(dx[!,x]),ln)
                #col_annotations = (median(Vector(dx[!, 1])), y_offset, xans) #middle of x-axis
                col_annotations = (median(Vector(dx[!, 1])), y_offset, xans) #middle of x-axis
                #col_annotations = (y_offset, y_offset, xans)
                
                # annotations = [(
                #     x_position, y_offset, 
                # Plots.text(label, fnt)) 
                # for (x_position, label) in 
                # zip(y_offset, anns)]



                p1 = @df dx StatsPlots.plot(
                    Vector(dx[!, 1]),
                    cols(ln),
                    legend = leg,
                    title = ti,
                    seriestype = mode,
                    annotations = col_annotations)  # Add annotations to the plot
                #return p1

                # for annotation in annotations
                #     Plots.annotate!(annotation)
                # end


            end

            # if names(dx)[1] == "date"
            #     col_annotations = nothing
            # end
            
            if names(dx)[1] == "month"
                Plots.xticks!(
                    1:12, monthabbr.(1:12))
            end

            if log === true
                yaxis!(:log10)
            end

            if title !== true
                title!(raw"")
            end

            return p1
        end

        function qplot(df1::DataFrame,df2::DataFrame;col1=1,col2=1)
            """
            takes two dfs and plots r2 QQ of selected cols
            """
            if any(map(x->contains(x,"date"),names(df1)))
                df1 = df1[!,Not(:date)]
            end

            if any(map(x->contains(x,"date"),names(df2)))
                df2 = df2[!,Not(:date)]
            end

            try
                df1 = select!(df1, col1)
                df2 = select!(df2, col2)
            catch
                @error "col not found!"
                return
            end

            

            # if ncol(df)>2
            #     df = df[!,1:2]
            # end

            r2 = round(cor(df1[!,1], df2[!,1])^2, digits=3)
            
            p = qqplot(df1[!,1], df2[!,1], 
                #title = "R² = "*string(ti),
                qqline = :fit)
                #color = :grays) # erstellt ein QQ-Diagramm <- black
            xlabel!(p,names(df1)[1])
            ylabel!(p,names(df2)[1])
            annotate!(p,:bottomright, text("R² = "*string(r2), :black))
                    
        end

        function tocb(s::Union{String,Regex})
            """
            greps first from current dir Regex and copies to clipboard
            """
            if isa(s,Regex)
                s = Regex(s,"i")
            end
            y = first(filter(file -> occursin(s,file), readdir()))
            println("abspath of $y in clipboard!")
            abspath(y)|>cb
        end


        function readroute(x::Union{String,AbstractString})
            """
            reads from wasim routing table 
            usually named "route.txt"
            from ... 
            infile = ctl()|>first|>split|>last|>x->replace(x,r"ctl.*"=>"ctl")
            ofl = "route.txt"
            routeg(infile, ofl)
            df = readroute(ofl)
            ...        
            """
            df = CSV.read(x,DataFrame,header=false,
                skipto=8,delim="\t",footerskip=1,lazystrings=false)
            rename!(df,1=>"sim",2=>"obs",3=>"name")
            df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
            df.name=map(x->replace(x,r"_>.*" => ""),df.name)
            sort!(df, :sim)
            return df
        end

        function fp(x::Regex;)
            """
            """
            x = glob(x)
            filter!(z->endswith(z,".nc"),x)
            x = first(x)
            r = read(Raster(x,missingval=0;lazy=true))
            println(descr(r))
            Plots.plot(
                r;
                #c=cgrad(:thermal),
                c=cgrad(:matter),
                xlabel="",
                ylabel="",
                title = name(r)
                )
        end

        function rename_columns(df::DataFrame, name_mapping::DataFrame)
            v = map(x->Grep.grep(x,names(df)),Regex.(string.(name_mapping.sim)))
            for (i, name) in enumerate(reduce(vcat, v))
                if occursin(string.(name_mapping.sim[i]),name)
                    println("renaming ",name, " <-> ",name_mapping.sim[i], " => " ,name_mapping.name[i])
                    rename!(df, name => name_mapping.name[i])
                end
            end
        end

        function tline(df::DataFrame, date_col::Symbol)
            """
            adds trendlines to plot
            """
            # Get the date column and column names for trendlines
            date_data = df[!, date_col]
            trendline_cols = setdiff(names(df), [string.(date_col)])
        
            p = Plots.plot()
            
            for y_col in trendline_cols
                y_data = df[!, y_col]
        
                # Perform linear regression to get the slope and intercept
                X = hcat(ones(length(date_data)), 1:length(date_data))
                y = y_data
                β = X \ y  # Linear regression
        
                # Extract the intercept and slope
                intercept, slope = β[1], β[2]
        
                # Generate the trendline using the linear equation
                trendline(x) = intercept + slope * x
        
                # Add the trendline to the plot
                plot!(p ,date_data, 
                    trendline.(1:length(date_data)), 
                    #label="$y_col Trendline", 
                    label=false,
                    linewidth=2, linestyle=:dash)
            end
        
            return p
        end

        function dtrange(start, stop, step)
            """
            dtrange(Date(2004),Date(2005),Month(4))
            """
            collect(start:step:stop)
        end

        
        function dtrange(start, stop, step)
            """
            dtrange(Date(2004),Date(2005),Month(4))
            """
            collect(start:step:stop)
        end

        function cdof(x::Union{String,DataFrame})
            if x isa DataFrame
                try
                    d = collect(DataFrames.metadata(x))[1][2]|>dirname
                    cd(dirname(d))
                catch
                    @error "no basename in $x !"
                end
            else
                cd(dirname(x))
            end
            d=pwd()
            println("current dir: $d")
        end

        function read_df(s::String)
            """
            reads a DataFrame from a file w dlm and tryparse subsetting
            """
            data, colnames = DelimitedFiles.readdlm(s, '\t', String, '\n', header=true)
            df = DataFrame(Matrix{Any}(data), :auto)
            rename!(df,Symbol.(colnames)|>vec)
            df = df[map(x->!isnothing(x),tryparse.(Int,df[!,1])),:]
            for i in 1:size(df, 2)
                df[!, i] = map(x -> parse(Float64, x), df[!, i])
            end
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
            df = df[:,Not(1:4)]
            metadata!(df, "filename", s, style=:note);
            return df
        end

        function kernelplot(df::Union{String,DataFrame})
            """
            using KernelDensity
            usage: kernelplot("route.txt")
            """
            if df isa String
                ofl = "route.txt"
                df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
                rename!(df,1=>"sim",2=>"obs",3=>"name")
                df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
                df.name=map(x->replace(x,r"_>.*" => ""),df.name)
            end
        
        
            #dsu = wa.qbb()|>last #recursive
            dsu = wa.qba()
            M = parse.(Float64,dsu[!,Cols(2,4)])|>y->subset(y,1=> ByRow(>(0.)))|>Matrix
            B = kde(M)
            plot(B);
            dsu.basin  .= df.name
            name_mapping = df    
            msk = @rsubset(parse.(Float64,dsu[!,Cols(2,4)]) .> 0)
            nmn = dsu[msk[!,1],:]
            nmn.basin=map(x->replace(x,r"_" => " "),nmn.basin)
            Plots.annotate!([(M[i,1], M[i,2], 
                Plots.text(
                nmn.basin[i], 
                8, :left, halign=:center, rotation=-15.0)) for i in 1:size(nmn, 1)])
            plot!(legend=false)
        end

        function corrbar(a::Vector{Float64}, b::Vector{Float64})
            # for specific subsets of dfs....
            df_A = a
            df_B = b
            #ti = "$a vs $b"
            ti = raw""

            #cor(df_A, df_B)^2
        
            #Calculate correlations and replace NaN with 0
            correlations = Vector{Float64}(undef, size(df_A,1))
            for i in 1:size(df_A, 1)
                correlations[i] = cor(df_A, df_B)^2
            end
            replace!(correlations, NaN => 0)
        
            p0 = Plots.bar(1:size(df_A, 1), correlations,
                legend = false,
                title = ti,
                fillcolor = ifelse.(correlations .> 0.35, 
                    "cornflowerblue", "coral2"),
                xticks = (1:size(df_A, 1), propertynames(df_A)),
                xrotation = 45,
                xlabel = "",
                ylabel = "Correlation R²",
                left_margin = 10mm,
                bottom_margin = 2mm);
        
            ann = map(x->string.(round(x;sigdigits=2)),correlations)
        
            for i in 1:size(df_A, 1)
                Plots.annotate!(i,correlations[i],
                (ann[i],9,:center,:top,:black))
                # println("R² "*ann[i]*" of Basin 
                # "*(df_A)[i]*" added")
            end
        
        return p0
        end
        
        #function cmplot(;temp::Regex,prec::Regex,col=1)
        function cmplot(;temp=r"^so_temper",prec=r"^so_prec",col=1)
            """
            selects first dfcol (for so)
            dt annotations
            wa.cmplot(;temp=r"^so_temper",prec=r"^so_prec")
            see also:
            climateplot(r"^tem",r"^pre")
            col = subbasin of interest
            wa.climateplot(r"^temp",r"pre";col="tot_average")
            ws. prc and temp tauschen un opacity einstellen....
            """
            temp,prec,col = (r"^so_temper",r"^so_pre",1)
            #eig. braucht man das col argument bei so_ input nicht.

            prec = waread2(prec)
            
            yrs = year.(prec.date)|>unique|>length
            prec = monsum(prec)
            if ncol(prec) > 2
                select!(prec, Not(:month))
                prec = prec[:,col]
                precvec = vec(Matrix(prec))
            else
                precvec = vec(Matrix(select(prec, Not(:month))))
            end
            
            
            precvec = precvec ./ yrs
            
            temp = waread2(temp)|>monmean
            if ncol(temp) > 2
                select!(temp, Not(:month))
                temp = temp[:,col]
                tempvec = vec(Matrix(temp))
            else
                tempvec = vec(Matrix(select(temp, Not(:month))))
            end
            
            month_abbr = ["Jan", "Feb", "Mär", "Apr", 
                "Mai", "Jun", "Jul", "Aug", "Sep", 
                    "Okt", "Nov", "Dez"];

            p1 = Plots.bar(prec.month, precvec, 
                color=:cornflowerblue, 
                #xlabel="", 
                xflip=false,
                ylabel="Niederschlag [mm]", 
                legend=false, yflip=true);
            xticks!(1:12, month_abbr)
            for i in prec.month
                val = round(precvec[i]; digits=1)
                annotate!(i, precvec[i], text("$(val)",8,:bottom))
            end

            # val = round(maximum(tempvec); sigdigits=2)
            # i = findmax(tempvec)[2]
            # Plots.annotate!(i, tempvec[i], 
            #             text("max: $(val) °C",10,:top))
            p2 = twinx();
            ann2 = map(x->
                text(
                string.(round(x; digits=1))*"°", 
                    7,:left, 
                    :red),
                    tempvec)
            
            plot!(p2, tempvec, xlabel="", 
                ylabel="Temperatur [°C]", color=:coral2,
                annotations = (temp.month .+ 0.125, 
                    tempvec .+ 0.1, ann2, :center),
                label=false, 
                linestyle = :dashdot,
                linewidth=3);      

            return p1
        end

        function selt(x::DataFrame, col::Any)
            """
            selects from dataframe date an column
            """
            df = select(x,Cols(col,:date))
            println(names(df))
            return df   #vec(Matrix(df))
        end

        function hydromon(x::Union{Regex,String,DataFrame}; leg = :outertopright, logy = false)
            """
            hydrgraph monthly mean plot
            """
            if isa(x,DataFrame)
                df = (x)
            else
                df = waread(x)
            end
            
            ti = try
                DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "No basename in metadata!"
                raw""
            end
            s = Symbol.(filter(x -> !occursin(r"date|year|month"i, x), names(df)))
            years = unique(Dates.year.(df[!, :date]))
            ylog = logy ? :log : :identity
        
            begin
                su = @rsubset df year.(:date)==years[1]
                su = monmean(su)
                @df su Plots.plot(:month,cols(s),
                    label = years[1], 
                    title=ti,
                    yaxis=ylog,
                    legend=leg)
                for yr in years[2:end]
                    su = @rsubset df year.(:date)===yr
                    su = monmean(su)
                    @df su Plots.plot!(:month,
                        cols(s),
                        yaxis=ylog,
                        label = yr)
                end
                month_abbr = ["Jan", "Feb", "Mär", "Apr", 
                            "Mai", "Jun", "Jul", "Aug", "Sep", 
                                "Okt", "Nov", "Dez"];
                #xticks!(0.5:11.5 , month_abbr)
                xticks!(1:12, month_abbr)
                plot!()
            end
        end

        function hydro(x::Union{Regex,String,DataFrame}; 
                col = 1, #Union{Symbol,Int64}, 
                leg = :outertopright, logy = false)
            """
            hydrgraph plot, selection of first column by default
            """
            if isa(x,DataFrame)
                df = (x)
            else
                df = waread(x)
            end
            
            ti = try
                DataFrames.metadata(df)|>values|>only
            catch
                @warn "No basename in metadata!"
                raw""
            end

            #df = selt(df,col)
            df = select(df,Cols(col,:date))
            println("Selection:",names(df))

            colname = names(df[!,Not(:date)])|>only
            #ti = colname*" - "*ti
        
            s = filter(x -> !occursin(r"(?i)date|year|month", 
                string(x)), names(df))
            years = unique(year.(df.date))
            years_str = string.(years)
        
            ylog = logy ? :log : :identity
            
            mn = [ monthabbr(x) for x in unique(month.(df.date)) ]
        
            hp1 = plot(#xticks = mn,
                    xrot = 45,
                    xminorticks = false,
                    xmajorticks = false,
                    yaxis = ylog,
                    #xlim = extrema(df.date),
                    title = ti,
                    legendtitle = colname,
                    #formatter = :plain,
                    legend = leg)
            #title!("Legende", legendtitle = "Titel")
        
            for yr in years #[2:end]
                su = filter(row -> year(row.date) == yr, df)
                hp1 = plot!(vec(Matrix(select(su, Not(:date)))),
                            label = yr) #,formatter = :plain)
            end
            lng = size(df,1) ./ length(years) ./ 12
            nln = size(df,1) ./ length(years)
            #lng = size(df,1) ./ 2
            #monthday.(df.date)|>unique #das müsste besser gehen...
            #st = 15:31:size(df,1)
            #st = 15:lng:size(df,1)
            st = 15:lng:nln
            xticks!(st, mn)
            #xlabel!(mn;xrot=45)
            #xlabel!(string(mn),rotation=45)
        
            return hp1
        end

        cdinto = cdof

        function strans(fn::String;src=EPSG(25832),dst=EPSG(4326))
            """
            reads from stationdata, reprojects to EPSG
            #ProjString("+proj=longlat +datum=WGS84 +no_defs")
            """
            fl = CSV.read(fn,DataFrame;limit=4)
            xc = fl[2,5:end]|>collect
            yc = fl[3,5:end]|>collect
            pts = ArchGDAL.IGeometry[]
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                pt = ArchGDAL.reproject(pt,src,dst)
                push!(pts,pt)
            end
            od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
            return od
            
        end

        function getnames(dfs::DataFrame,rg::String)
            nms = names(dfs)
            m = filter(x->occursin(Regex(rg,"i"),x),nms)
            m = length(m)===1 ? only(Symbol.(m)) : Symbol.(m)
            return(m)
        end

        function bargroup(x::Union{Regex,String,DataFrame};leg = :topright,fun=wa.monsum, lcols=1)
            """
            tst of grouped barplot
            """
            if isa(x,DataFrame)
                df = (x)
            else
                df = waread(x)
            end

            v = map(
                (x->occursin(r"date", x) & !occursin(r"year|month", x)),
                (names(df))
                )

            if fun==false
                df = hcat(df[:,Cols(r"date")],df[!,Not(Cols(r"date"))])
                @info "no function applied!\n"
            else
                @info "applying $fun ..."
                df = fun(df)
            end

            if ncol(df)<2
                @error "only one column available!"
                return
            end
            
            s = Symbol.(filter(x -> !(occursin(r"year|date|month"i, x)), names(df)))
            dt = first(Symbol.(filter(x -> (occursin(r"year|date|month"i, x)), names(df))))

            ti = try
                DataFrames.metadata(df)|>values|>only|>basename
            catch
                @warn "No basename in metadata!"
                raw""
            end
            xt = wa.tovec(df,dt)
            mx = round(findmax(df[!,Not(Cols(dt))]|>Matrix|>collect)[1];digits=0)
            minx = round(findmin(df[!,Not(Cols(dt))]|>Matrix|>collect)[1];digits=0)
            
            step = nrow(df) #.* 2
            lab = string.(minx:step:mx)
            lab = lab.*" [mm]"
            @df df groupedbar(xt,cols(s), 
                legend = :outertopright,
                xticks = xt,
                xrotation = 45,
                xlabel = "", 
                ylabel = "[mm]", 
                
                #yticks = (minx:step:mx,lab),
                # guidefonthalign = :right,
                # yguidefonthalign = :right,
                # yguidefontrotation = -90,
                #yguidefontvalign = :top,
                #guide = "[mm]", #both axis labels
                yguidefontsize = 8,
                minorgrid = true,
                tick_direction = :out,
                legendtitle = ti*"\n",
                margin = 0.5cm, #def 1mm
                legend_column = lcols,
                legend_font_pointsize = 6,
                #legend_font = Plots.Font("sans-serif", 6, 
                #    :hcenter, :vcenter, 0.0,colorant"black"),
                #legend_title_font_valign = :top,
                #legend_position = :outertopright, #`:best`
                #legend_title_font_halign = :left
                ## see https://docs.juliaplots.org/latest/generated/attributes_subplot/
                )
        end

        function stplot!(fn::Union{String,DataFrame})
            """
            reads from stationdata, reprojects to 4326 and plots to existing
            also check addplot
            """
            if isa(fn,DataFrame)
                fl = fn
                xc = fl.xc
                yc = fl.yc
                pts = ArchGDAL.IGeometry[]
                for i in 1:length(xc)
                    pt = ArchGDAL.createpoint([xc[i],yc[i]])
                    pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                    push!(pts,pt)
                end
                #trim chars to length 10
                len = 10
                mns = map(y->((length(y)> len)  ? y[1:len] : y),string.(fl.name))
                od = DataFrame(geometry=pts, name=mns, xc=xc, yc=yc)
            else
                fl = CSV.read(fn,DataFrame;limit=4)
                xc = fl[2,5:end]|>collect
                yc = fl[3,5:end]|>collect
                pts = ArchGDAL.IGeometry[]
                for i in 1:length(xc)
                    pt = ArchGDAL.createpoint([xc[i],yc[i]])
                    pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                    push!(pts,pt)
                end
                od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
                
            end
            
            Plots.plot!(od.geometry;c=:black,shape = :x)
            for (i, pt) in enumerate(od.geometry)
                x = ArchGDAL.getx(od.geometry[i], 0)
                y = ArchGDAL.gety(od.geometry[i], 0)
                name = od.name[i]
                annotate!(x, y, text(name, 8, :black, :bottom, :left))
            end
            Plots.plot!()
        end

        function pe()
            try
                inp = clipboard()
                wpath = replace(inp, "\\" => "/", "\"" => "")
                cmd = `wslpath -ua $wpath`
                ot = readchomp(pipeline(cmd))
                clipboard("$ot")
                println("$ot in clipboard!")
            return string(ot)
            catch e
                println("Error: $e")
                println("Failed to translate to wslpath.")
            end
        end

        function pew()
            try
                in = clipboard()
                wpath = replace(in, "\\" => "/")
                println("$wpath in clipboard!")
                clipboard("pt=$wpath")
                return wpath
            catch e
                @error "smth errored $e"
            end
        end
        
        function isroute()
            """
            checks state of wasim routing table
            """
            pt=src_path*"/win/routeg.jl"
            include(pt)
        end

        function lpro(x::Union{String,DataFrame})
            if x isa String
                fl = CSV.read(x,DataFrame;limit=4)
            else
                fl = x
            end
            xc = fl[2,5:end]|>collect
            yc = fl[3,5:end]|>collect
            pts = ArchGDAL.IGeometry[]
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                pt = ArchGDAL.reproject(pt,
                GeoFormatTypes.EPSG(25832),
                GeoFormatTypes.EPSG(4326))
                push!(pts,pt)
            end
            od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
        
            return od    
        end

        function fread(z::Union{Regex,String})
            """
            only for wasim output files
            """
            if isa(z,Regex)
                v = filter(file -> occursin(z,file), readdir());
                z = v[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt|log"i),v)]|>first
            end 
            println("loading $z ...")   
            m = map(x->(x=>Int64),[:YY,:MM,:HH,:DD])
            #df = CSV.read(z,DataFrame;types=Dict(m))|>f->dropmissing(f,:YY)
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
        
        function rsq(df::DataFrame;x1=1,x2=2)
            x=df[:,x1]
            y=df[:,x2]
            cor(x, y)^2
        end
        
        function dubplot(pt::Union{String,DataFrame})
            
            if pt isa String
                df = waread2(pt)
            else
                df = pt
            end
                
            
            #v = nonunique(df,:date)
            yl = try
                collect(DataFrames.metadata(df))[1][2]|>dirname|>splitpath|>last    
            catch
                @info "No basename in metadata!"
                raw""
            end
            
            
            dubs = df[findall(x -> count(==(x), df.date) > 1, df.date),:]
            grp = groupby(dubs,:date)
            v = [size(yx,1) for yx in grp]|>unique
            println("unique counts: $(v)")
            if length(v) == 0
                @info("Warning: no unique count!")
                return
            else
                println("unique counts: $(v)")
            end
            
            plot()
            
            # i = 0;
            for group in grp
                # cntz = size(group,1)
                # i =+ 1
                @df group violin!(
                    cols(propertynames(dubs)[1]), 
                    ylabel=yl, 
                    xlabel="Group", 
                    legend=false)
            end
            for i in 1:size(grp, 1)
                Plots.annotate!(i,
                maximum(grp[i][!,1]),    
                (size(grp[i],1),9,:center,:top,:black))
            end
            title!("Bandwith of duplicate values")
        end

        """
        plots a dataframe with info on dublicated date vals
        """
        function dubplot(pt::Union{String,DataFrame})
            
            if pt isa String
                df = waread2(pt)
            else
                df = pt
            end
                   
            
            #v = nonunique(df,:date)
            yl = try
                collect(DataFrames.metadata(df))[1][2]|>dirname|>splitpath|>last    
            catch
                @info "No basename in metadata!"
                raw""
            end
            
            
            dubs = df[findall(x -> count(==(x), df.date) > 1, df.date),:]
            grp = groupby(dubs,:date)
            v = [size(yx,1) for yx in grp]|>unique
            println("unique counts: $(v)")
            if length(v) == 0
                @info("Warning: no unique count!")
                return
            else
                println("unique counts: $(v)")
            end
            
            StatsPlots.plot()
            
            # i = 0;
            for group in grp
                # cntz = size(group,1)
                # i =+ 1
                @df group StatsPlots.violin!(
                    cols(propertynames(dubs)[1]), 
                    # annotations = Plots.text("no: $cntz"),
                    # annotations = (i,
                    # maximum(group[!,1]),
                    # "no: $cntz", :top),
                    ylabel=yl, 
                    xlabel="Group", 
                    legend=false)
            end
            for i in 1:size(grp, 1)
                Plots.annotate!(i,
                maximum(grp[i][!,1]),    
                (size(grp[i],1),9,:center,:top,:black))
            end
            title!("Bandwith of duplicate values")
        end
    
        """
        prints pretty_table of df
        and copies to clipboard
        kwars are passed to pretty_table
        tblcb(kd;backend = Val(:latex))
        """
        function tblcb(kd::DataFrame;kwargs...) 
            pretty_table(kd,header=uppercasefirst.(names(kd)); kwargs...)
            txt_str = sprint(io -> pretty_table(io, kd, header=uppercasefirst.(names(kd)); kwargs...))
            return clipboard(txt_str) 
        end
    
        """
        applies to quarterly data
        qrtr(pt::Union{String,DataFrame},
            fun=sum;agg=quarterofyear)
        example 
         wa.qrtr(df;agg=year)|>dfp
         wa.qrtr(df;fun=mean)|>dfp
         wa.qrtr(wa.skipyr(df);fun=mean)|>dfp
        """
        function qrtr(pt::Union{String,DataFrame};fun=sum,agg=quarterofyear)
    
            if pt isa String
                x = waread2(pt)
            else
                x = pt
            end
    
            yl = try
                collect(DataFrames.metadata(df))[1][2]|>dirname|>splitpath|>last    
            catch
                @info "No basename in metadata!"
                raw""
            end
            
            df = copy(x)
            y = filter(x->!occursin("date",x), names(df))
            s = map(y -> Symbol(y),y)
            df[!, :date] .= agg.(df[!,:date]);
            df_agg = DataFrames.combine(groupby(df, :date), 
                y .=> fun .=> y);
            return(df_agg)
        end
    
        """
        sim obs plot rglob regex
        x::Union{Regex,DataFrame}; 
            yscale::Symbol = :log,
            fun::Function = sum, 
            freq::String="monthly",simcol=1,obscol=2
            freq can also shortened, like D,M,Q,Y
        """
        function hyeval(
            x::Union{Regex,DataFrame}; 
            yscale::Symbol = :log,
            fun::Function = sum, 
            freq::String="monthly",
            ylab::String="[mm/$freq]",
            simcol=1,obscol=2)
            if x isa Regex
                x = first(dfonly(x))
                ndf = waread2(x;silencewarnings=true)
            else
                ndf = reorder_df(x)
            end
                    
            ndf = select(ndf,Cols(simcol,obscol,r"date|month|year"))
            dropmissing!(ndf)
    
            @info "aggregation (sum) freq: $freq"
            try 
            # Resample the DataFrame based on the specified freq
                if freq == "daily" || freq == "D" || freq == "day"
                    #ndf = ndf[hour.(ndf.date) .== 0, :]
                    ndf = ndf
                elseif freq == "monthly" || freq == "M" || freq == "mon" || freq == "month"
                    #df[day.(df.date) .== 1, :]
                    ndf = qrtr(ndf;fun=fun,agg=month) #monsum(ndf)
                elseif freq == "quarterly" || freq == "Q" 
                    #df[day.(df.date) .== 1 && month.(df.date) % 3 == 1, :]
                    ndf = qrtr(ndf;fun=fun)
                elseif freq == "seasonal"
                    #df[day.(df.date) .== 1 && month.(df.date) % 3 == 1, :]
                    ndf = qrtr(ndf;fun=fun)
                elseif freq == "yearly" || freq == "Y" || freq == "yr" || freq == "year"
                    #df[day.(df.date) .== 1 && month.(df.date) == 1, :]
                    ndf = qrtr(ndf;fun=fun,agg=year)  #yrsum(ndf)
                    #DataFrames.combine(groupby(df, year.(df.date)), y .=> mean .=> y);
                    #DataFrames.combine(groupby(df, quarterofyear.(df.date)), y .=> mean .=> y);
                end
            catch e
                @error("smth went wrong",e)
                return
            end
    
            ndf = hcat(ndf[!,Not(Cols(r"date|month|year"))],
                ndf[:,Cols(r"date|month|year")])
        
            #printstyled(names(ndf)...,bold=true,color=:red)
            #split(names(ndf)...,' '),bold=true,color=:red)
            printstyled(names(ndf),bold=true,color=:red)
            
            rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
            
            printstyled(" changed to :Simulated,:Observed,:Date ! \n",bold=true,color=:red)
            
            overall_pearson_r = cor(ndf[!, :Observed], 
                ndf[!, :Simulated])
            r2 = overall_pearson_r^2
            #nse(simulations, evaluation)
            nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
            kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
            ve = round(vef(ndf[!, :Simulated],ndf[!, :Observed]), digits=2)
            subs = "RSQ: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))\nVE: $ve"
            
    
            ti = try 
                first(split(basename(x),"_"))
            catch
                @warn "No basename in metadata!"
                raw""
            end
            
            fr = replace(freq,"ly"=>"")
            p = Plots.plot(title=ti, ylabel=ylab, 
                xlabel="", 
                yaxis = yscale, 
                legend=:topleft)
    
            Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], 
                    color=:red, label="Modeled")
            Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], 
            line=:dash, color=:blue,label="Observed")
                
            
            Plots.annotate!(
                :topright,
                text("$subs", 10, :black, :right)
            )
    
            if freq == "quarterly" || freq == "Q" 
                Plots.xticks!(
                    1:4, ["Q1", "Q2", "Q3", "Q4"])
            end
                    
            if freq == "monthly" || freq == "M" || freq == "mon" || freq == "month"
                Plots.xticks!(
                    1:12, monthabbr.(1:12))
            end
            
            if freq == "yearly" || freq == "Y" || freq == "yr" || freq == "year"
                Plots.xticks!(
                    #1:nrow(ndf), 
                    ndf.Date, 
                    string.(ndf.Date))
            end
            
    
            return p
        end
    
        """
        finds all ctl files in a directory -> NO subdirectories by default
        findctl(snippet::Union{String,Regex};recurse=false, dir="D:/Wasim/regio/control",suffix=".ctl")
        """
        function findctl(snippet::Union{String,Regex};recurse=false, dir="D:/Wasim/regio/control",suffix=".ctl")
            owd = abspath(pwd())
            platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
            
            if platform != "windows"
                wsl_cmd = `wslpath $(dir)`
                dir = readchomp(pipeline(wsl_cmd))
            end
            
            cd(dir)
            if recurse
                #files = rglob(Regex(suffix,"i"))
                res = []
                for (looproot, dirs, filenames) in walkdir(dir)
                    for filename in filenames
                        if (occursin(snippet,filename) && endswith(filename, suffix))
                            push!(res, joinpath(looproot, filename)) 
                        end
                    end
                end
            else
                files = filter(file -> endswith(file, suffix), readdir(dir,join=true))
                res = filter(x->occursin(snippet,x),files)
            end    
            
            cd(owd)
            return res
        end
    
        """
        sums up Conda.ROOTENV
        win.julia_conda:            3.06 GB
        """
        function condasize()
            cwd=Conda.ROOTENV
            osize = 0
            n = 0
            for (root, dirs, files) in walkdir(cwd)
             for file in files
                 osize += stat(joinpath(root, file)).size
                 n += 1
             end
             for dir in dirs
                printstyled("check dir: $dir\n",color=:light_red)
             end
            end 
            println("$(n) files in directory")
            @printf("%-40s %15.2f GB\n","$(cwd):",osize/1024^3)
        end
    
        """
        daily df input, yearly kge nse ve as DF
        outd = readall(r"qoutjl\$")
        ks = map(byear,outd)
    
        for z in dfonly(r"qoutjl\$")
            wa.byear(z)|>println
        end    
        """
        function byear(x::Union{String,Regex,DataFrame};)
            if x isa String
                printstyled("reading $x\n",color=:light_red)
                df = waread2(x;silencewarnings=false)
                dropmissing!(df)    
            elseif x isa Regex
                x = first(dfonly(x))
                printstyled("reading $x\n",color=:light_red)
                df = waread2(x;silencewarnings=false)
                dropmissing!(df)
            else
                #df = copy(x)
                df = reorder_df(x)
                dropmissing!(df)
            end
            
            df[!, :year] = year.(df[!,:date]);
            nm = collect(DataFrames.metadata(df))[1][2]|>basename;
            nm = replace.(nm,
                r"_" => " ",
                r"-qoutjl|qout" => "",
                r"#" => "")
    
            grouped_df = groupby(df, :year)
            DataFrames.combine(grouped_df) do group
                simulated, observed = vec(Matrix(group[!,Cols(1)])),vec(Matrix(group[!,Cols(2)]))
                kge = wa.kge2(simulated, observed)
                ve = wa.vef(simulated, observed)
                nse = wa.nse(simulated, observed)
                
                #grouping key :year is returned as first column
                dout = DataFrame(kge=kge,nse=nse,ve=ve,nm=nm)
                return dout
            end
        end
    
        """
        lat lon reverse
        """
        function reverse_coords(polygon)
            # Extract the coordinates from the polygon
            # crds = GeoInterface.coordinates.(polygon)
            # # Reverse each pair of coordinates
            # reversed_crds = [[[(lon, lat) for (lat, lon) in ring] for ring in polygon] for polygon in crds]
            if polygon isa DataFrame
                out = []
                for i in 1:length(polygon.geometry)
                    # Get the coordinates of the polygon
                    coords = GeoInterface.coordinates(polygon.geometry[i])
                    # Flatten the list of points
                    ptc = coords[1]
                    ptc_tuples = [tuple(pt...) for pt in ptc]
                    prev = [(Y,X) for (X,Y) in ptc_tuples]
                    reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
                    push!(out,reversed_polygon)       
                end
                pout = vcat(out...)   
                return pout
            end
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
    
        function heat(x::Union{String,Regex,DataFrame};ann=true,kw...)
            if x isa String
                printstyled("reading $x\n",color=:light_red)
                df = waread2(x;silencewarnings=false)
                dropmissing!(df)    
            elseif x isa Regex
                x = first(dfonly(x))
                printstyled("reading $x\n",color=:light_red)
                df = waread2(x;silencewarnings=false)
                dropmissing!(df)
            else
                df = copy(x)
                dropmissing!(df)
            end
    
            if ncol(df)<2
                @error "only one column available!"
                return
            end
    
            if ncol(df)>20
                @warn "more than 20 columns, heatmap will be unreadable!"
                return
            end
            #md = cor(Matrix(df[:, Not(:date)])).^2
            #md = Matrix(select(df, Not(Cols(r"date|year|month|day"))))
            md = cor(Matrix(select(df, Not(Cols(r"date|year|month|day"))))).^2
            md = convert(Matrix{Union{Float64, Missing}}, md)
            #md = cor(md)^2
            replace!(md, 1.0 => missing)
            replace!(md, 0.0 => missing)
    
            p = heatmap(md, 
                #c=:balance, 
                #c=:lightrainbow, grep(r"gr",cs)
                #c=:grays,
                #c = Plots.colormap("RdBu",nrow(df);logscale=true),
                c = Plots.colormap("RdBu",size(md, 1),mid=0.2),
                linewidth=0.8, 
                cbar=true,
                colorbar_title = "RSQ",
                margins=5mm,
                xticks =false,
                yticks =false,
                grid = :false,
                # legend = :outertopright,
                # legendtext = string(names(df[:, Not(:date)])),#xlabel="", #ylabel="", 
                title="Correlation Heatmap";kw...)
            if ann
                for i in 1:size(md, 1)
                    for j in 1:size(md, 2)
                        value = round(md[i, j], digits=2)
                        color = ismissing(value) ? :transparent : :black
                        #color = value==1.0 ? :transparent : :black
                        Plots.annotate!(j, i, 
                            Plots.text(string(value), 9, color, :center, 
                            halign=:center, rotation=-35.0))
                        # annotate!(j - 0.5, i,
                        #     text(value, 
                        #         8, color, :center))
                    end
                end
            end
    
            nm = names(df[!,Not(Cols(r"date|year|month|day"))])
            nm = replace.(nm,
                r"_" => " ",
                r"qoutjl|qout" => "",
                #"-qoutjl"=>"",
                r"#" => "")
            Plots.xticks!(1:length(nm), nm) #reverse(nm)
            Plots.yticks!(1:length(nm), nm;rotation=-35)
    
            return(p)
        end
    
        """
        returns the position of each element in a vector
        """
        function findalls(x::Vector)
            map(e -> e => findall(==(e), x), unique(x))
        end
    
        function tbx(df::DataFrame)
            ti = try
                DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end
        
    
            month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
            if any(x->occursin("year",x),names(df))
                #df = df[!,Not("year")]
                s = Symbol.(filter(x->!occursin(r"date|year|month|day",x),names(df)))
                #str = unique(df.year)
                @df df StatsPlots.violin(cols(s), linewidth=0.1)
                #xticks!(0.5:11.5, month_abbr)
                title!(ti)
            else
                str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                fnt = Plots.font(
                    family="sans-serif",
                    pointsize=8,
                    valign=:bottom,
                    rotation=0.0,
                    color=:black
                )
        
                dx = wa.monsum(df)
                y_offset = dx[!,2]
                anns = map(x->string.(round(x,digits=2)), y_offset)
                xans = map(x->Plots.text(x, fnt), anns)
                #col_annotations = (Vector(dx[!, 1]) .- 0.5, y_offset, xans) # x y val
                #col_annotations = (Vector(dx[!, 1]) .- 0.5, y_offset, xans) # x y val
                col_annotations = (Vector(dx[!, 1]) .+ 0.5, y_offset, xans) # x y val
        
                # Create a colormap from the `monmean` values
                #k=wa.colorfunction(dx[!,2])
                #k = cf2(dx[!,2])
                #mat = reshape(k, (1,nrow(dx)))
                #colors = colormap.(dx[!,2])
        
                @df df StatsPlots.boxplot(
                    str,
                    cols(s),
                    group = str, # group by month, needed for colormap
                    linewidth=0.1,
                    outliers = false,
                    whisker_width = :match,
                    #colors=map(x -> (x > 3 ? :green : :red), y_offset),
                    #fillcolor = mat|>reverse,
                    #fillcolor = Colors.gray.(1:12),
                    #fillcolor = :grey,
                    colors = Plots.colormap("RdBu",size(dx, 2),mid=0.2),
                    annotations = col_annotations,
                    legend=false,
                )
        
                xticks!(0.5:1.4:16, month_abbr)
                #xticks!(0.5:1.5:12.5, month_abbr)
                title!(ti)
            end
        end
    
        """
        moisture_plot_with_confidence(df, timestep=month; kwargs...)
        """
        function moisture_plot_with_confidence(df, timestep=month; kwargs...)
            # Extracting relevant columns
            dk = select(df, [:date, :tot_average])
            if timestep != :day
                dk[!, :agg_date] .= timestep.(dk[!,:date])
            else
                dk[!, :agg_date] .= dk[!,:date]
            end    
    
            # Grouping by the specified timestep and calculating the minimum and maximum values of tot_average
            su = DataFrames.combine(groupby(dk, :agg_date), 
                :tot_average .=> (minimum, maximum) .=> 
                [:min_moisture, :max_moisture])
    
                lab = uppercase(string(timestep))
            # Creating the plot with confidence bands
            plot(su.agg_date, [su.min_moisture su.max_moisture], ribbon=true, fillalpha=0.2,
                xlabel=string(lab), 
                ylabel="Total Average Moisture",
                title="Moisture Plot with Confidence Bands",
                label=["Min Moisture" "Max Moisture"]; kwargs...)
        end
    
        function moncol(x::DataFrame)
            dmean = copy(x)
            s = map(Symbol, filter(x -> !occursin("date", x), names(dmean)))
            dmean[!, :month] = month.(dmean[!,:date])
            select!(dmean, Not(:date))
        
            # Calculate monthly mean and confidence interval using combine and groupby
            dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
                s .=> (dmean -> (mean(dmean), 1.96 * std(dmean) / sqrt(length(dmean)))) 
                .=> s)
        
            # Create separate columns for mean, lower bound, and upper bound
            for col in s
                dmean[!, string(col)*"_mean"]  .= [first(only(x)) for x in eachrow(dmean[!, col])]
                dmean[!, string(col)*"_lower"] .= [first(only(x)) - last(only(x)) for x in eachrow(dmean[!, col])]
                dmean[!, string(col)*"_upper"] .= [first(only(x)) + last(only(x)) for x in eachrow(dmean[!, col])]
            end
        
            # # Drop the original columns representing the mean and standard deviation
            # select!(dmean, Not.(:month, s...))
        
            return dmean
        end
    
        """
        aggregate and add confidence intervals
        df::DataFrame; confidence_level=0.95
        """
        function monc(df::DataFrame; confidence_level=0.95)
            # Make a copy of the input DataFrame
            dmean = copy(df)
            
            # Extract the columns of interest
            columns = filter(x -> !occursin("date", x), names(dmean))
            
            # Add a 'month' column based on the 'date' column
            dmean[!, :month] = month.(dmean[!, :date])
            
            # Drop the 'date' column
            select!(dmean, Not(:date))
            
            # Calculate monthly mean and confidence interval using combine and groupby
            dmean = DataFrames.combine(DataFrames.groupby(dmean, :month)) do group
                result = DataFrame(month = group.month)
                for col in columns
                    values = group[!, col]
                    mean_value = mean(values)
                    confidence_interval = (mean_value, quantile(Normal(), 0.5 + confidence_level / 2) * std(values) / sqrt(length(values)))
                    mv = string(col) * "_mean"
                    lb = string(col) * "_lower"
                    ub = string(col) * "_upper"
                    result.mean .= mean_value
                    result.lower .= confidence_interval[1] - mean_value
                    result.upper .= confidence_interval[2] - mean_value
                end
                return result
            end
        
            return dmean
        end
    
        function eeread(filename)
            df = CSV.read(filename, DataFrame)
            df = unique(df, 1) # remove duplicate date rows
            rename!(df, 1 => :date)
            df.date .= [Date(d, dateformat"u d, y") for d in df.date]
            z = basename(filename)
            DataFrames.metadata!(df, "filename", z, style=:note)
            return df
        end
    
        """
        monthly mean plot with ribbon
        """
        function dfrib(x::Union{String,Regex,DataFrame};col::Any=1)
            if x isa String
                printstyled("reading $x\n",color=:light_red)
                df = waread2(x;silencewarnings=false)
                dropmissing!(df)    
            elseif x isa Regex
                x = first(dfonly(x))
                printstyled("reading $x\n",color=:light_red)
                df = waread2(x;silencewarnings=false)
                dropmissing!(df)
            else
                df = x #reorder_df(x)
                dropmissing!(df)
            end
    
            ti = try
                DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "No basename in metadata!"
                ti = raw""
            end       
    
            if !("date" in names(df))
                @error "no :date column found in $x !"
                return
            end
    
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
    
            #colname = string(names(df)[1]) #Not(:date)
            #colname = string(names(df[!, Not(:date)]))
            ##has to be in this order...
            colname = only(names(select(df,Cols(col))))
            if length(ti)>1
                lab = ti[1:4]*colname
            else
                lab = colname
            end
    
            df = select(df,Cols(col,:date))
            @info "aggregation monthly mean of col: $col"
            # if startswith(colname,"_")
            #     colname = "basin"*colname
            # end
    
            # #dmean = copy(df)
            # dmean = df
            # s = map(Symbol, filter(x -> !occursin("date", x), names(dmean)))
            # dmean[!, :month] = month.(dmean[!,:date])
            # select!(dmean, Not(:date))
            
            # # # Calculate monthly mean and confidence interval using combine and groupby
            # # dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
            # #     s .=> (dmean -> (mean(dmean), 1.96 * std(dmean) / sqrt(length(dmean)))) 
            # #     .=> s)
    
            # # Calculate monthly mean and narrower confidence interval using combine and groupby
            # dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
            #     s .=> (
            #         dmean -> 
            #         (mean(dmean), 1.1 * std(dmean) / sqrt(length(dmean)))) 
            #     .=> s)
        
        
            # # Create separate columns for mean, lower bound, and upper bound
            # for col in s
            #     dmean[!, string(col)*"_mean"]  .= [first(only(x)) for x in eachrow(dmean[!, col])]
            #     dmean[!, string(col)*"_lower"] .= [first(only(x)) - last(only(x)) for x in eachrow(dmean[!, col])]
            #     dmean[!, string(col)*"_upper"] .= [first(only(x)) + last(only(x)) for x in eachrow(dmean[!, col])]
            # end
    
            dmean = wa.monc(df) #s.o. confidence_interval 95%
        
            plt = @df dmean plot(:month, Matrix(dmean[!,Cols(r"mean")]), 
                #ribbon=(dmean.Lai_lower, dmean.Lai_upper), 
                ribbon=(Matrix(dmean[!,Cols(r"lower")]), Matrix(dmean[!,Cols(r"upper")])), 
                labels = lab, #see above at ti   # "Lai_mean", 
                xlabel = "", #"Month", 
                ylabel = "",  #"Lai", 
                title = ti, #"Lai Mean with Ribbon",
                legend = :topleft,
                xticks = (1:12, 
                [ monthabbr(x) for x in unique(dmean.month) ]),
                xrotation = 35
                )    
            return plt
        end
    
        
        """
        usage like: kernelplot("route.txt") but with linlog option
        reads with  wa.qba() and names of df or (route.txt)
        lin: cols=Cols(2,4)
        log: cols=Cols(3,5)
        """
        function klog(df::Union{String,DataFrame};lin=false)
            if lin 
                cols=Cols(2,4)
            else
                cols=Cols(3,5)
            end
    
            if df isa String
                ofl = df #should be "route.txt"
                df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
                rename!(df,1=>"sim",2=>"obs",3=>"name")
                df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
                df.name=map(x->replace(x,r"_>.*" => ""),df.name)
            end
            #dsu = wa.qbb()|>last #recursive
            dsu = wa.qba()
            M = parse.(Float64,dsu[!,cols])|>y->subset(y,1=> ByRow(>(0.)))|>Matrix
            mns = names(dsu[!,cols])
            B = kde(M)
            Plots.plot(B,legend=false,xlabel=mns[1],ylabel=mns[2])
            dsu.basin  .= df.name
            name_mapping = df    
            msk = @rsubset(parse.(Float64,dsu[!,cols]) .> 0)
            nmn = dsu[msk[!,1],:]
            nmn.basin=map(x->replace(x,r"_" => " "),nmn.basin)
            Plots.annotate!([(M[i,1], M[i,2], 
                Plots.text(
                nmn.basin[i], 
                8, :left, 
                halign=:center, 
                rotation=-15.0)) for i in 1:size(nmn, 1)])
            Plots.plot!(legend=false)
        end
    
    
        """
        pycall function to polygonize a raster
        polygonize_raster(input_raster_path::String, output_shapefile_path::String;epsg=25832)
        """
        function polygonize_raster(input_raster_path::String, output_shapefile_path::String;epsg=25832)
            gdal = pyimport("osgeo.gdal")
            ogr = pyimport("osgeo.ogr")
            osr = pyimport("osgeo.osr")
    
            # Open the raster dataset
            dataset = gdal.Open(input_raster_path)
    
            # Get the first band
            band = dataset.GetRasterBand(1)
    
            # Get the "ESRI Shapefile" driver
            driver = gdal.GetDriverByName("ESRI Shapefile")
    
            # Create a new shapefile dataset
            out_ds = driver.Create(output_shapefile_path, 0, 0, 0, gdal.GDT_Unknown)
    
            # Create a spatial reference object
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(epsg)
    
            # Create a new layer
            layer = out_ds.CreateLayer("polygonized", srs, ogr.wkbPolygon)
    
            # Polygonize the raster
            try
                gdal.Polygonize(band, py"None", layer, -1, [], callback=py"None")
                @info "new shapefile created at: $output_shapefile_path ..."
            catch
                @error "gdal.Polygonize failed ..."
            end
            
            # Close the dataset to write it to the disk
            out_ds = py"None"
        end
    
    
    
        """
        reads controlfile
        uses Grep.grep to select lines
        returns a DataFrame
        """
        function ctlook(infile::String)
            data = open(filename) do io
                readbetween(io, "soil_table", "substance_transport")
            end
            data = data[3:end]
            numbers = extract_numbers(data)
            # Remove comments and unnecessary characters
            data = broadcast(x -> replace(x,    
                    r"^#.*" => "",
                    r"^[[].*" => "",
                    r"method" => "",
                    r"MultipleHorizons" => "",
                    r"}" => "",
                    r" = " => "=",
                    #r"[?*.{=;]" => "",      #problems
                    r";" => "" ), data)
            data = Grep.grep(r"^(?! Evap.*$|^[0-9].*$)",data)
            data = strip.(data)
            filter!(x -> !isempty(x), data)
            sel = Grep.grep(r"^Name", data)
            # Split each string into "Name" and "Value"
            split_name_value = split.(sel, '=', limit=2)
            
            ks = split.(Grep.grep(r"^ksat", data), '=', limit=2)
            ths = split.(Grep.grep(r"^theta_sat", data), '=', limit=2)
            thr = split.(Grep.grep(r"^theta_res", data), '=', limit=2)
            thk = split.(Grep.grep(r"^thickness", data), '=', limit=2)
            parn = split.(Grep.grep(r"^Par_n", data), '=', limit=2)
            #ks = Dict(first(ks)[1] => getindex.(ks, 2))
            # Create a DataFrame
            df = DataFrame( id = numbers,
                            Name = getindex.(split_name_value, 2),
                            ksat = getindex.(ks, 2),
                            theta_sat = getindex.(ths, 2),
                            theta_res = getindex.(thr, 2),
                            Par_n = getindex.(parn, 2),
                            Thick = getindex.(thk, 2))
            return df
        end
    
        """
        reads controlfile
        select landuse table
        returns a Vector of DataFrames
        """
        function read_landuse_data2(filename::AbstractString)
            data = open(filename) do io
                readbetween(io, "landuse_table", "special_output")
            end
            data = broadcast(x -> replace(x,
            r"\\t" => "",        #strip comments
            r"^#.*" => "",        #strip comments
                                r";.*$" => "",         # match everything after semicolon until the end the line
                                r" = " => "=",      #strip spaces around "="
                                r";" => "" ), data)
                                #r"^[[].*" => "",      #strip module names
                                #r"method" => "",     #startflag
                                #r"}" => "",          #this is needed for end flag
                                
                                #r"[?*.{=;]" => "",      #problems
    
    
            class_dataframes = []
            io = open(filename, "r")
    
            for line in data
                if occursin(r"(?i)method", line)
                    classmatch = strip(line[3:15])  # first 10 chars, stripping leading/trailing whitespaces
                    println(classmatch)
                    filtered_lines = readbetween(io, Regex(classmatch), r"}$")
                    filtered_lines = replace.(filtered_lines, r";.*" => "")
                    filtered_lines = filter(x -> !isempty(x), filtered_lines)
                    # Add a header
                    header = ["Parameter=Value"]
                    filtered_lines = vcat(header, filtered_lines)
                    # Create a CSV File from the filtered lines
                    csv_file = CSV.File(IOBuffer(join(filtered_lines, '\n')), delim='=', quotechar=' ', ignorerepeated=true)
                    # Convert the CSV File to a DataFrame
                    class_df = DataFrame(csv_file)
                    push!(class_dataframes, class_df)
                end
            end
            close(io)
            return class_dataframes
        end
    
        """
        scatterplot of landusedata
        dfs::Vector=dfs,x::String="rs_evaporation"
        """
        function luscatter(dfs::Vector=dfs,x::String="rs_evaporation")
            month_abbr = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]
            p1 = plot(title=x);
            for i in 1:size(dfs,1)
                value_str = findindf(dfs[i], x).Value
                #value_str = replace.(value_str, r" " => "")
                value_str = string.(value_str)
                #value_str = strip.(value_str)
                vs = split.(value_str)
                z0 = [parse.(Float64, z) for z in vs]
                cls = split(dfs[i][1, 1], r"\W")[2] 
                scatter!(month_abbr, first(z0), label=cls,
                    yaxis=:identity) #log10
                #push!(plots, plot(month_abbr, only(z0), label=cls))
            end
            #title!("LAI in [m2/m2]")
            Plots.show(p1)
            return p1    
        end
    
        """
        heatmap of landusedata 
        dfs::Vector=dfs,x::String="LAI"
        """
        function luheat(dfs::Vector=dfs,x::String="LAI";colors=:grays)
            month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
            ti = findindf(dfs[1],Regex(x,"i"))[1,1]|>first|>strip
            #ti = findindf(dfs[1],Regex(x,"i")).Parameter|>first|>strip
            a = map(k->findindf(k,Regex(x,"i")),dfs)
            a = vcat(a...)|>DataFrame
            z = vcat(map(x->split(replace(
                string(x[1, 1]),
                r"\s+"=>" ",
                r"{"=>"",
                r"method"=>"")),
                dfs)...)
            a.nm = filter(x->length(x)>2,z)
            a = a[:, sort(names(a),rev=true)] #sort columns
            a.nm .= replace.(a.nm, r"_" => " ")
            begin 
                value_str = string.(a.Value)
                vs = split.(value_str)
                z0 = [parse.(Float64, z) for z in vs]
                m = hcat(z0...)'
                #myc = colormap("Grays",10;logscale=true)
                #myc = colormap("Grays",Int(round(findmax(m)[1]));logscale=false)
                #heatmap(m,c=:grays)
                heatmap(m,c=colors)
                #heatmap(m,c=myc[5:end-1])
                xticks!(1:12,month_abbr)
                yticks!(1:size(m,1),a.nm)
                
                title!("$ti")
    
                #add annotation for each field
                qv = quantile(m)[4] #|>round
                for i in 1:size(m,1)
                    for j in 1:size(m,2)
                        #tcolor = m[i,j] < qv ? :black : :white
                        tcolor = m[i,j] < qv ? :white : :black
                        annotate!(j,i,text(string(
                            round(m[i, j]; digits=2)
                            ),7,tcolor,rotation=-45))
                    end
                end
                return plot!()
            end
        end
        
        
    
        """
        runs perl wsl \$bsp x x 4 5 or 
        \$f1 \$f2 \$c1 \$c2 ∇
        pers(x::Union{String,Regex}=r"qoutjl";c1=4,c2=5,f2=nothing,tofile=false)
        """
        function pers(x::Union{String,Regex}=r"qoutjl$";c1=4,c2=5,f2=nothing,tofile=false)
            bsp = "/mnt/c/Users/Public/Documents/Python_Scripts/bashscripts/pearson_wsl.pl"
            dirpt = last(split(dirname(pwd()),"\\"))
            res = "pers_"*dirpt*".txt"
            if f2 != nothing
                f1 = isa(x,String) ? x : first(x)
                if tofile
                    #res = readchomp(pipeline(cmd))
                    #run(`wsl sh -c "(echo bla) | tee z.txt" `)
                    run(`wsl sh -c "$bsp $f1 $f2 $c1 $c2 | tee $res" `)
                else
                    run(`wsl $bsp $f1 $f2 $c1 $c2`)
                end
            else   
            files = glob(x)
                for x in files
                    if tofile
                        res = "pers_"*dirpt*".txt"
                        res = replace(res,"pers"=>x)
                        #run(`wsl sh -c "($cmd) | tee $res" `)
                        run(`wsl sh -c "$bsp $x $x $c1 $c2 | tee $res" `)
                    else
                        run(`wsl $bsp $x $x $c1 $c2`)
                    end
                    
                end
            end
        end
    
        """
        takes Vector{DataFrame} from:
        dataframes = map(byear,outd)
        plot_grouped_metrics(dataframes::Vector{DataFrame};col=:ve,all=false,kw...)
        """
        function plot_grouped_metrics(dataframes::Vector{DataFrame};col=:ve,usethresold=true,threshold = -0.41,all=false,kw...)
            #nam = wa.getnames(dataframes)
            # Create a new plot
            plot()
    
            # Iterate through each DataFrame in the input
            for df in dataframes
                # Extract year and metrics from the DataFrame
                year = df.year
                if all
                    nam=collect(DataFrames.metadata(df))[1][2]|>basename
                    nam=replace(nam,"-qoutjl"=>"",r"_" => " ",r"#" => "",r"qout" => "")
                    # if usethresold
                    #     kge = ifelse.(df.kge .<= threshold, -1, df.kge)
                    #     nse = ifelse.(df.nse .<= threshold, -1, df.nse)
                    #     ve = ifelse.(df.ve .<= threshold, -1, df.ve)
                    #     yaxis!((threshold,1))
                    # else
                    #     kge = df.kge
                    #     nse = df.nse
                    #     ve = df.ve
                    # end
                    if usethresold
                        yaxis!((threshold,1))
                    end
                    kge = df.kge
                    nse = df.nse
                    ve = df.ve
                    # Plot each metric with a different color
                    plot!(year, kge, label="KGE $nam", marker=:circle)
                    plot!(year, nse, label="NSE $nam", marker=:square)
                    plot!(year, ve, label="VE $nam", 
                        marker=:diamond,
                        legend=:outertopright,
                        kw...)
                    title!("Grouped Metrics by Year")
                else
                    #dfilter = DataFrames.subset(df, col => ByRow(<=(threshold)))|>Matrix|>vec
                    #DataFrames.transform(df, names(df)[2:end-1] .=> cor,renamecols=false)
                    # if usethresold
                    #     df[!, col] .= ifelse.(df[!, col] .<= threshold, -1, df[!, col])
                    #     yaxis!((threshold,1))
                    # end
                    if usethresold
                        yaxis!((threshold,1))
                    end
                    ve = select(df,col)|>Matrix|>vec
                    ti = select(df,col)|>names|>first|>uppercase
                    nam=collect(DataFrames.metadata(df))[1][2]|>basename
                    nam=replace(nam,r"-qoutjl"=>"",r"_" => " ",r"#" => "",r"qout" => "")
                    plot!(year, ve, label=nam, 
                        title =ti, 
                        marker=:diamond,
                        legend=:outertopright,
                        kw...)
                    Plots.xticks!(df[:, :year])
                end
            end
            # Customize the plot
            #xlabel!("Year")
            
            xx = uppercase(string(col))
            ylabel!("Score $xx")
            # Add annotation from DataFrames.metadata filename
            #annotate!(Plots.text(size=(10, 10), metadata_filename, 2012, 1.5, :center, :center))
        end
    
        """
        takes Vector{DataFrame} from:
        dataframes = map(byear,outd)
        annoated scatterplot of metrics
        pxm(dataframes::Vector{DataFrame}; col=:ve, usethresold=true, threshold=-0.410, all=false, kw...)
        """
        function pxm(dataframes::Vector{DataFrame}; col=:ve, usethresold=true, threshold=-0.410, all=false, kw...)
            p1=Plots.plot()
            # Iterate through each DataFrame in the input
            for df in dataframes
                # Extract year and metrics from the DataFrame
                year = df.year
                if all
                    nam = collect(DataFrames.metadata(df))[1][2] |> basename
                    nam = replace(nam, "-qoutjl" => "", r"_" => " ", r"#" => "", r"qout" => "")
                    if usethresold
                        yaxis!((threshold, 1))
                    end
                    kge = df.kge
                    nse = df.nse
                    ve = df.ve
                    # Plotting with annotations
                    Plots.scatter!(year, kge, label="KGE $nam", marker=:circle)
                    for (x, y) in zip(year, kge)
                        if y > threshold
                            Plots.annotate!([(x, y, 
                            Plots.text("$(round(y, digits=2))", 8, 
                            :black, :top))])
                        end
                    end
                    Plots.scatter!(year, nse, label="NSE $nam", marker=:square)
                    for (x, y) in zip(year, nse)
                        if y > threshold
                            Plots.annotate!([(x, y, 
                            Plots.text("$(round(y, digits=2))", 8, 
                            :black, :top))])
                        end
                    end
                    Plots.scatter!(year, ve, label="VE $nam", marker=:diamond, legend=:outertopright, kw...)
                    for (x, y) in zip(year, nse)
                        if y > threshold
                            Plots.annotate!([(x, y, 
                            Plots.text("$(round(y, digits=2))", 8, 
                            :black, :top))])
                        end
                    end
                    title!("Grouped Metrics by Year")
                    Plots.xticks!(df[:, :year])
                    ylabel!("Score")
                else
                    if usethresold
                        yaxis!((threshold, 1))
                    end
                    
                    ve = select(df, col) |> Matrix |> vec
                    ti = select(df, col) |> names |> first |> uppercase
                    nam = collect(DataFrames.metadata(df))[1][2] |> basename
                    nam = replace(nam, r"-qoutjl" => "", r"_" => " ", r"#" => "", r"qout" => "")
                    
                    # Plotting with annotations
                    Plots.scatter!(year, ve, 
                        #label=nam, 
                        label=false,
                        title=ti, 
                        marker=:diamond, legend=:outertopright, kw...)
                    
                    # Adding annotations for values above threshold
                    for (x, y) in zip(year, ve)
                        if y > threshold
                            Plots.annotate!([(x, y, 
                            Plots.text("$(round(y, digits=2)) "*nam[1:3], 
                            8, :black, :top))])
                        end
                    end
                    
                    Plots.xticks!(df[:, :year])
                end
                xx = uppercase(string(col))
                ylabel!("Score $xx")
            end
            return p1
        end
        
    
    
        """
        merges and saves to qoutjl
        outd = pout(infile;sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/",ofl = "route.txt")
        """
        function pout(infile;sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/",ofl = "route.txt")
            routeg(infile, ofl)
            sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
            specfile=joinpath(sfpt,sfn)
            obs = readdf(specfile)
            df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
            rename!(df,1=>"sim",2=>"obs",3=>"name")
            df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
            df.name=map(x->replace(x,r"_>.*" => ""),df.name)
            sort!(df, :sim)
            sim = r"qges"|>glob|>first|>waread
            @info "innerjoin has to be sim -> obs!"
            outd = []
            for i in eachrow(df)
                println(i[1],"->",i[3])
                try
                    dm =innerjoin(
                        sim[!,Cols("C"*string(i[1]),end)],
                        obs[!,Cols(Regex(i[3],"i"),end)],    
                        on=:date)
                    onam = i[3]*"-qoutjl"
                    wawrite(dm,onam)
                    println("$onam saved!")
                    DataFrames.metadata!(dm, "filename", onam, style=:note);
                    push!(outd,dm)
                    println(names(dm)," on $onam pushed!")
                catch
                    onam = i[3]*"-qoutjl"
                    @warn "merge is empty for $onam ! ..."
                    continue
                end
            end
            vz = filter(xx->xx[2]==5,wa.cntcolv("outjl"))
            if length(vz)>0
                @warn "some files are not properly merged..."
            end
            return outd
        end
    
        """
        dfroute(;ofl="route.txt")
        reads from wa.routeg(infile, ofl) and returns a DataFrame with the following columns:
            - sim: simulated flow
            - obs: observed flow
            - name: name of the station
        """
        function dfroute(;ofl="route.txt")
            df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
            rename!(df,1=>"sim",2=>"obs",3=>"name")
            df.name=map(x->replace(x,r"#" => "",r" " => "",r"_>.*" => "",r"-" => "_"),df[:,3])
            sort!(df, :sim)
            return df    
        end
    
        """
        The F1 Score ranges from 0 to 1, where a higher value indicates a better balance between precision and recall
        F1= 2⋅(Precision⋅Recall) / Precision+Recall 
        """
        function f1_score(data::DataFrame, 
            sim_col::Union{Symbol,Int64} = 1, 
            obs_col::Union{Symbol,Int64} = 2; threshold = 0.5, minobs=true)
            # Extracting simulation and observation columns
            #sim_values = data[:, sim_col]
            #obs_values = data[:, obs_col]
            sim_nm=only(names(select(data, sim_col)))
            obs_nm=only(names(select(data, obs_col)))
            @info "simcol is $sim_nm, obscol is $obs_nm"
            sim_values = select(data, sim_col)|>Matrix|>vec
            obs_values = select(data, obs_col)|>Matrix|>vec
    
            if minobs
                threshold = minimum(obs_values)
                @info "threshold is set to minimum of observations: $threshold"
            else
                # Threshold for classification (you may need to adjust this based on your problem)
                threshold = threshold
            end
    
            # Convert to binary classification (1 if greater than threshold, 0 otherwise)
            sim_binary = sim_values .> threshold
            obs_binary = obs_values .> threshold
    
            # True Positives, False Positives, True Negatives, False Negatives
            tp = sum(sim_binary .& obs_binary)
            fp = sum(sim_binary .& .!obs_binary)
            tn = sum(.!sim_binary .& .!obs_binary)
            fn = sum(.!sim_binary .& obs_binary)
    
            # Precision, Recall, and F1 Score
            precision = tp / (tp + fp)
            recall = tp / (tp + fn)
            
            # F1 Score is the harmonic mean of precision and recall
            f1_score = 2 * (precision * recall) / (precision + recall)
    
            return f1_score
        end
    
        """
        recursive rmeq; keeps otherdata
        """
        function rmopt()
            # Get the list of files
            files = wa.rglob("")
    
            #!occursin(r"pl|sh|csv|html|xml|fzt|ftz|log|ini|^wq|yrly|nc|tif|jpeg|png|svg|txt", file)
            # py|R|ftz_0|tex
            # Define the file extensions to keep
            keep_exts = [".ipynb", ".py", ".R", ".Rmd", ".log",
                ".tif", ".jpeg", ".png", ".svg",
                ".cpg", ".shx", ".dbf", ".prj", ".shp", ".tex", 
                ".csv", 
                ".html", ".ftz", ".ftz_0", ".txt", 
                ".list", ".nc", ".xml", ".sh", ".grd", ".yrly"]
                
            files = filter(file -> 
                !occursin(r"(wq_|pl$|sh$|fzt|ftz|log$|ini|otherdata|intern)", file)
                , files)
            
                # Process each file
            for file in files
                # Check if the file extension is in the list of extensions to keep
                if !(splitext(file)[2] in keep_exts)
                    try
                    # Read the file as a DataFrame
                    println("Processing $file ...")
                    ms = ["-9999", "lin", "log", "--"]
                    #df = CSV.read(file, DataFrame)
                    df = CSV.File(file; 
                        delim="\t", 
                        header=1,
                        silencewarnings=true, 
                        normalizenames=false, 
                        missingstring=ms, 
                        types=Float64) |> DataFrame
                    dropmissing!(df,ncol(df))
                    # dropmissing!(df,4)
                    # Calculate the sums of the 5th column and the last column
                    y = sum(df[:, 5])
                    x = sum(df[:, end])
    
                    # Check the conditions
                    if x <= (nrow(df) - 5) * -9999 && y <= (nrow(df) - 5) * -9999 || x == 0 && y == 0
                        # Remove the file
                        rm(file)
                        printstyled("$file removed ...\n", color=:red)
                    end
                    catch e
                        println("Error processing $file: $e")
                    end
                end
            end
        end
    
    
        
#    end ##end of precompile
    
end ##end of module