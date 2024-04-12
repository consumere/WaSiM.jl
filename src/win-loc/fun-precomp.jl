#functions

#cd("/home/ubu/.julia/packages/WaSiM")
#generate WaSiM
#cp -v /mnt/c/Users/Public/Documents/Python_Scripts/julia/fun-precomp.jl /home/ubu/.julia/packages/WaSiM/src/WaSiM.jl
#using PkgTemplates; 
#t = Template(; user="consumere", authors=["christian-schäfer"], plugins=[License(name="MIT"), Git(), GitHubActions()], )

module WaSiM
    # include("src/WaSiM.jl")
    # end
    
    using SnoopPrecompile    # this is a small dependency
    import Printf, DataFrames, CSV, Statistics, Dates, Rasters, Distributions, StatsPlots, PlotlyJS
    
    @precompile_setup begin
        # Putting some things in `setup` can reduce the size of the
        # precompile file and potentially make loading faster.
        #list = [OtherType("hello"), OtherType("world!")]
        @precompile_all_calls begin
            # all calls in this block will be precompiled, regardless of whether
            # they belong to your package or not (on Julia 1.8 and higher)
            # d = Dict(MyType(1) => list)
            # x = get(d, MyType(2), nothing)
            # last(d[MyType(1)])
    
            #using Query
    
            function setup()
                #include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/functions.jl")
                include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/fun-precomp.jl")
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
    
            function qgk(rootdir=".", prefix="")
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
    
            #Note that the function assumes that P is defined 
            #and represents the atmospheric pressure (kPa).
    
    
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
            #    dirs = readdir(".")
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
                 #if isfile(file) && contains(file, ext)
                 if isfile(file) && occursin(Regex(ext),file)
                 #if isfile(file) && occursin(ext,file)
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
    
            #list(x) = Any[i for i ∈ x]
    
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
    
    
            function denseplot(regex::AbstractString,dfs::Vector{DataFrame})
                "selects first match and plots..."
                df = dfs[map(n->occursin(Regex(regex,"i"),n),
                     map(x->basename(only(DataFrames.metadata(x))[2]),
                     dfs))] |> first
                     s = propertynames(df)[Not(end)];
                     o = DataFrames.metadata(df)|>collect
                     ti = basename(o[1][2])
                     @df df density(cols(s),title=ti,legend = :topright) 
            end
    
    
            function denseplot(df::DataFrame)
                #println(propertynames(df))
                s = propertynames(df)[Not(end)] #masks last column == date     #[1:end-1]
                #,propertynames(df)[end]
                o = DataFrames.metadata(df)|>collect
                ti = basename(o[1][2])
                @df df density(cols(s), legend = :topright, title=ti)
            end
    
            function denseplot(df::String)
                df=readdf(df)
                s = propertynames(df)[Not(end)]
                o = DataFrames.metadata(df)|>collect
                ti = basename(o[1][2])
                @df df density(cols(s), legend = :topright, title=ti)
            end
    
            function dfp(df::DataFrame)
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
    
            # plotf("tem")
    
            function readdf(x::Regex)
                """
                readdf(x::Regex)
                reads first match of regex wasim timeseries
                """
                rootdir="."
                results = []
                for (looproot, dirs, filenames) in walkdir(rootdir)
                    for filename in filenames
                        if (occursin(x,filename)) && (!occursin(r"txt|yrly|nc|png|svg|grd",filename))
                            push!(results, joinpath(looproot, filename)) 
                        end
                    end
                end
                x = first(results)
                ms=["-9999","lin","log"]
                df = CSV.read(x,DataFrame,missingstring=ms,
                ntasks=4,
                limit=typemax(Int),
                types = Float64,
                delim="\t",
                silencewarnings=true,
                normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
                df.YY=map(x ->Int(x),df.YY);
                df.MM=map(x ->Int(x),df.MM);
                df.DD=map(x ->Int(x),df.DD);
                df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
                df=df[:,Not(1:3)]
                metadata!(df, "filename", x, style=:note);
            end
    
    
            function readdf(x::AbstractString)
        #ms="-9999"
        ms=["-9999","lin","log"]
        df = CSV.read(x,
        DataFrame,
        missingstring=ms,
        #skipto=4, 
        ntasks=4,
        limit=typemax(Int),
        types = Float64,
        delim="\t",
        #comment=r"[A-z]",
        silencewarnings=true,
        normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,Not(1:3)]
        metadata!(df, "filename", x, style=:note);
            end
    
            readmeteo = readdf
            loaddf = readdf
            # function readmeteo(x::AbstractString)
            #     df = DataFrame(CSV.File(x; missingstring="-9999",
            #                         skipto=6,
            #                         limit=typemax(Int),
            #                         comment="#",
            #                         stringtype=String))
            #     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
            #     df=df[:,Not(1:4)]
            #     metadata!(df, "filename", x, style=:note);
            # end
    
            # function loaddf(path::AbstractString)
            #     ms=["-999","-9999","lin","log","LIN","LOG"]
            #     df = CSV.read(path,DataFrame,
            #     missingstring=ms,
            #     delim="\t",comment="-",
            #     silencewarnings=false,
            #     ntasks=4,downcast=true,	    
            #     normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
            #     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
            #     df=df[:,Not(1:3)]
            #     metadata!(df, "filename", path, style=:note);
            # end
    
            # mx=ct("so_a")
            # x=mx[2]
            # df = CSV.read(x, DataFrame,missingstring=["-9999"], delim="\t",
            # skipto=4, #bei nich so: 3
            #         #comment="[A-z]",
            #         comment="-",
            #         silencewarnings=true) |>dropmissing
            # df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
            # df = select!(df, Not(names(df)[1:4]))
            # plotf(df)
    
            function loadso(path::AbstractString, prefix::AbstractString)
                files = readdir(path)
                dfs = DataFrame[]
                for file in files
                    if isfile(file) && occursin(Regex(prefix),file) && (!occursin(r"txt|yrly|nc|png|svg|ftz_0|ftz",file))
                       file_path = joinpath(path, file)
            	   println("reading ",file_path,"...")
            	   p1 = readdf(file_path)
            	   push!(dfs, p1)
                    end
                end
                return(dfs)
            end
    
            #md = loadso(pwd(),"so")
            #map(names,md)
            #size(md) 
            #for i in md; print(size(i)) ; end
            #broadcast(x -> innerjoin(x,on = :date),md)
            #innerjoin(md,on = :date)
    
            # using StatsPlots;
            # dfs = loadso(pwd(),"te")
            # vibx(dfs[2])
            # vibx(dfs[3])
    
            # size(dfs)
            # i=dfs[16]
            # s = propertynames(i)[1:end-1]
            # @df i plot(:date,cols(s))
    
            # Plots.gr()
    
    
            #by(df, :month, broadcast(x->ln[x]=>sum,ln))
    
            # df[!, :month] = month.(df[!, :date]);
            # by(df, :month, ln[1]=>sum)
            # for i in ln;println(by(df, :month, i=>sum));end     
    
            # u=by(df, :month, ln[1]=>sum);
            # nm = propertynames(u)[Not(1)]
            # #str = propertynames(u)[1]  
            # str = [ @sprintf("%02i", x) for x in u[1] ];
            # @df u StatsPlots.boxplot(str,cols(nm),fillalpha=0.75, linewidth=0.25);
    
            #broadcast(x->findall("date",x),names(df))  
            function toyrsum(df::DataFrame)
            od = []
            df[!, :year] = year.(df[!,:date]);
            ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            for i in ln;
            #x=(by(df,:year,i=>sum));
            x = combine(df, :year, AsTable(i) => sum, renamecols=false) 
            push!(od,x)
            end ;
            return(od)
            end
    
    
            # #od = DataFrame[]
            # df[!, :year] = year.(df[!,:date]);
    
            ##komplett mean:
            function fullmean(df::DataFrame)
            df[!, :year] = year.(df[!,:date]);
            combine(df, :, AsTable(Not([:date,:year])) => mean, renamecols=false)
            end
    
    
            function cattoyrsum(df::DataFrame)
            ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            df[!, :year] = year.(df[!,:date]);
            it=[]
            od=(by(df,:year,ln[1]=>sum));
            #od = combine(df, :year, AsTable(ln[1]) => sum, renamecols=false) 
            for i in ln[2:end];
            x=(by(df[Not(:date)],:year,i=>sum));
            #x = combine(df, :year, AsTable(i) => sum, renamecols=false)
            push!(it,x[end])      
            #push!(it,x[!,end])      
            end ;
            # x = combine(df, :year, AsTable(:) => sum, renamecols=false)
            # combine(groupby(df,:year),:=>sum)
            ot = hcat(od,DataFrame(it))  
            #ot=join(od,x[end],:year,makeuniuqe=true)
            #ot = hcat(od,it)
            return(ot)
            end
    
            # md=cattoyrsum(df)
            # @df md plot(:year,cols(propertynames(md)[2:end]),yaxis=:log)      
    
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
    
            # md=cattoyrmean(df)
            # @df md plot(:year,cols(propertynames(md)[2:end])) 
            # function tsyr(df::DataFrame)
            	# str = [ @sprintf("%02i", x) for x in (year.(df.date)) ];
            	# #tsn = "date";
            	# #ln = propertynames(df)[Not(tsn)];
            	# ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            	# @df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
            	# @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
            	# @df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
            # end
    
    
            function qpl(df::DataFrame)
    	ti = DataFrames.metadata(df)|>only|>last|>basename
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
        ti = DataFrames.metadata(df)|>only|>last|>basename
        s = names(df)[1:2]
    	t2 = string.(ti,"\n",s[1],"|",s[2],ti)
        StatsPlots.plot( 
    	qqplot(df[!,1],df[!,2], qqline = :fit), 
    	qqplot(Cauchy,df[!,2]), 
    	qqnorm(df[!,2], qqline = :R),
        title = t2)
            end
    
    
            function vibx(df::DataFrame)
            	str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            	#ln = propertynames(df[end-1])
            	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            	@df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
            	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
            	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
            end
    
            function vibx(df::String)
            	df = readdf(df)
                str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            	@df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
            	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
            	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.5,marker=(:black,stroke(1)),legend=false)
            end
    
            function vio(df::DataFrame)
            	str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            	#ln = propertynames(df[end-1])
            	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            	@df df StatsPlots.violin(str,cols(ln),linewidth=0.1)
            	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
            #	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
            end
    
            # gr()
            # default(show = true)
            # df = readdf(r"qgk")
            # nms = unique!(map(monthabbr,(month.(df.date))))
            # str = (map(monthabbr,(month.(df.date))))|>sort
            # #str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            # ln = Symbol.(filter(x->!occursin("date",x),names(df)))
            # ti=DataFrames.metadata(df)|>only|>last|>basename
            # #df[!, :year] = year.(df[!,:date]);
            # @df df StatsPlots.violin(str,cols(ln),
            #     xlabel="Months",
            #     #xaxis=nms,
            #     linewidth=0.1,
            #     title=ti)
            # #xlabel!(nms|>only)
            # @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false)
    
            # gr()
            # default(show = true)
            # df = readdf(r"qgk")
            # #nms = unique!(map(monthabbr,(month.(df.date))))
            # y = Symbol.(filter(x->!occursin("date",x),names(df)))
            # df[!, :month] = month.(df[!,:date]);
            # #regand(v=names(df),"date","month")
            # dfm = combine(groupby(df, :month), y .=> sum .=> y);
            # str = (map(monthabbr,(dfm.month)))
            # dfm.mab = str
            # #str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            # ln = Symbol.(filter(x->!occursin("month",x),names(dfm)))
            # ti=DataFrames.metadata(dfm)|>only|>last|>basename
            # @df dfm StatsPlots.violin(cols(ln[Not(end)]),
            #     xlabel="Months",
            #     linewidth=0.1,
            #     title=ti)
            # #xaxis!(dfm.mab)
    
            # #df[!, :year] = year.(df[!,:date]);
            # @df dfm StatsPlots.violin(str,cols(ln),
            #     xlabel="Months",
            #     #xaxis=nms,
            #     linewidth=0.1,
            #     title=ti)
            # #xlabel!(nms|>only)
            # @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false)
    
    
            	# t = month.(df.date);
            	# @df df StatsPlots.violin(  string.(broadcast(x -> Dates.monthname(x),t))), :tot_average, linewidth=0.01,legend=false);
            	# @df df StatsPlots.boxplot!(string.(broadcast(x -> Dates.monthname(x),t))), :tot_average, fillalpha=0.75, linewidth=0.25,legend=false);
            	# @df df StatsPlots.dotplot!(string.(broadcast(x -> Dates.monthname(x),t))), :tot_average, fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
            # end
    
            function vio(df::String)
                df = readdf(df)
            	str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
                #ln = Symbol.(filter(x->!occursin("date",x),names(df)))
                ti=DataFrames.metadata(df)|>only|>last|>basename
            	@df df StatsPlots.violin(str,cols(ln),linewidth=0.1,title=ti)
            	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.55, linewidth=0.25,legend=false);
            	#@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.15,legend=false)
                #marker=(:black,stroke(1)),legend=false)
            end
    
            function vio(regex::AbstractString,dfs::Vector{DataFrame})
                "selects first match and plots..."
                df = dfs[map(n->occursin(Regex(regex,"i"),n),
                     map(x->basename(only(DataFrames.metadata(x))[2]),
                     dfs))] |> first
                str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
                ln = Symbol.(filter(x->!occursin("date",x),names(df)))
                @df df StatsPlots.violin(str,cols(ln),linewidth=0.1)
                @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.55, linewidth=0.25,legend=true);
            end
    
    
            function ldf(path::AbstractString, prefix::AbstractString)
                files = readdir(path)
                dfs = DataFrame[]
            #    outname = []
                for file in files
                    #if isfile(file_path) && endswith(file, ext)
            	#if isfile(file_path) && (startswith(file, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",file))
                    if isfile(file) && occursin(Regex(prefix),file)&& (!occursin(r"txt|yrly|nc|png|svg",file))
                       file_path = joinpath(path, file)
            	   println("reading",file_path)
            	   p1 = loaddf(file_path)
            	   push!(dfs, p1)
            #	   m = match(r".*[.]",basename(file_path))
            #	   nm = contains(basename(file_path),".") ? string(m.match,"png") : basename(file_path)*".png"
            #	   push!(outname,nm)
                    end
                end
                return(dfs)
            #    return(outname)
            end
    
    
            # a = true
            # b = 1
            # c = 2
            # a ? b : c # 1
            #outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"
            #lyr=2
            #xlyr = length(lyr)!=1 ? 1 : lyr
    
    
            function pline(path::AbstractString)
                ms=["-999","-9999","lin","log","LIN","LOG"]
                df = CSV.read(path,DataFrame,
                #missingstring="-9999", #also windows
                missingstring=ms,
                delim="\t",comment="-",
                silencewarnings=false,
                ntasks=4,downcast=true, # got unsupported keyword arguments "ntasks", "downcast" @windows                                          
                normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
                df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
                df=df[:,Not(1:3)]
                # ms=["-9999","lin","log","LIN","LOG","--"] #comment="-",
                # #df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",comment="-",ignorerepeated=true,silencewarnings=true,typemap=Dict(Int64=>String))  |> @dropna() |> DataFrame
                # df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",ignorerepeated=true,silencewarnings=true,typemap=Dict(String=>Int64))
                # df = df[completecases(df), :]
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
                relayout!(p,height=600*1.5,width=900*1.5,title_text="Series of "*basename(path))
                p
            end
    
            # nrows=size(df)[2]-1
            # o = DataFrames.metadata(df)|>collect
            # ti = basename(o[1][2])
            # px = make_subplots(rows=1,cols=nrows, 
            # shared_xaxes=true,
            # shared_yaxes=true);
            # for i in 1:nrows;
            #     add_trace!(px, 
            #     PlotlyJS.scatter(x=df.date, y=df[:,i],
            #     name=names(df)[i]),   row=1,     col=i);
            # end
            # px
    
            # function dfplotjs(df::DataFrame)
            #     nrows=size(df)[2]-1 
            #     #length(names(df))-1
            #     o = DataFrames.metadata(df)|>collect
            #     ti = basename(o[1][2])
            #     fig = make_subplots(
            #         shared_xaxes=true, 
            #         shared_yaxes=true    
            #     #rows=2, cols=2
            #         );
            #     for i in 1:nrows;
            #         add_trace!(fig, 
            #         PlotlyJS.scatter(x=df.date, y=df[:,i],
            #         name=names(df)[i]));
            #     end
            #     fact=1;
            #     PlotlyJS.relayout!(fig,
            #     height=600*fact,width=900*fact,
            #     title_text="Series of "*ti)
            #     display(fig)
            # end
    
            # pt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/met0/temp_1970.txt"
            # pt="/mnt/d/Wasim/Goldbach/revision/fab150/gbfab150/qgesfab150.p1.2018"
            # df = readmeteo(pt)
    
            # function dfplotjs(df::DataFrame;log::Bool)
            #     nrows=size(df)[2]-1 
            #     #length(names(df))-1
            #     o = DataFrames.metadata(df)|>collect
            #     ti = basename(o[1][2])
            #     fig = make_subplots(
            #         shared_xaxes=true, 
            #         shared_yaxes=true    
            #     #rows=2, cols=2
            #         );
            #     for i in 1:nrows;
            #         add_trace!(fig, 
            #         PlotlyJS.scatter(x=df.date, y=df[:,i],
            #         name=names(df)[i]));
            #     end
            #     fact=1;
            #     PlotlyJS.relayout!(fig,yaxis_type="log",
            #     height=600*fact,width=900*fact,
            #     title_text="Series of "*ti)
            #     display(fig)
            # end
    
            #Wenn Sie ein Schlüsselwort-Argument nicht zugewiesen haben wollen, können Sie es weglassen oder nothing als Wert verwenden. Zum Beispiel:
            #function dfplotjs(df::DataFrame;logy::Any,fact::Any)
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
            # dfplotjs(df;logy=true,fact=0.66)
            # dfplotjs(df;logy=false,fact=0.7)
            # dfplotjs(df;logy=true)
            # dfplotjs(df)
    
            #ps=xx("qg")
    
            # dx=readdf(pt)
            # dfplotjs(df;logy=false,fact=0.7)
            function dfplotjs(df::AbstractString;logy::Bool,fact::Float64)
                df=readmeteo(df)
                nrows=size(df)[2]-1
                #length(names(df))-1
                o = DataFrames.metadata(df)|>collect
                ti = basename(o[1][2])
                fig = make_subplots(
                    shared_xaxes=true, 
                    shared_yaxes=true    
                    );
                for i in 1:nrows;
                    add_trace!(fig, 
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
    
            # Erstellen Sie einige zufällige Zeitreihendaten
            # dates = Date(2020):Day(1):Date(2021)
            # y1 = cumsum(randn(length(dates)))
            # y2 = cumsum(randn(length(dates)))
            # y3 = cumsum(randn(length(dates)))
            # y4 = cumsum(randn(length(dates)))
            # fig = make_subplots(rows=2, cols=2)
            # # Erstellen Sie einige Zeitreihendiagramme
            # trace1 = PlotlyJS.scatter(x=dates, y=y1)
            # trace2 = PlotlyJS.scatter(x=dates, y=y2)
            # trace3 = PlotlyJS.scatter(x=dates, y=y3)
            # trace4 = PlotlyJS.scatter(x=dates, y=y4)
            # add_trace!(fig, trace1, row=1, col=1)
            # add_trace!(fig, trace2, row=1, col=2)
            # add_trace!(fig, trace3, row=2, col=1)
            # add_trace!(fig, trace4, row=2, col=2)
            # display(fig)
    
    
    
            #px=PlotlyJS.scatter(x=df.date, y=df[:,Not(:date)])
            #PlotlyJS.relayout!(px;"Series of "*ti)
            #PlotlyJS.relayout!(px,height=600*1.5,width=900*1.5,title_text="Series of "*ti)
            #px.show()
    
            #PlotlyJS.plot(PlotlyJS.scatter(x=df.date, y=df[!,6]))
    
            #df = dataset(DataFrame,"iris")
            #plot(df,x=:sepal_length,y=:sepal_width,z=:petal_width,color=:species,type="scatter3d",mode="markers")
    
            function dfplot(df::DataFrame)
                nrows=size(df)[2]-1
                st=[]
                for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
                p = PlotlyJS.make_subplots(rows=nrows, cols=1, 
                shared_xaxes=true, 
                shared_yaxes=false,
                vertical_spacing=0.05,
                )
                for i in 1:nrows;
                        add_trace!(p, 
                        PlotlyJS.scatter(x=df.date, y=df[:,i],
                        name=st[i]),   row=i,     col=1);
                end
                PlotlyJS.relayout!(p,height=600*1.5,width=900*1.5)
                return(p)
            end
    
    
            function dfplot(df::AbstractString)
                df=readmeteo(df)
                nrows=size(df)[2]-1
                st=[]
                for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
                p = PlotlyJS.make_subplots(rows=nrows, cols=1, 
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
    
    
            function kge2(observed::Vector{Float64}, simulated::Vector{Float64})
                r = cor(observed, simulated)
                α = std(simulated) / std(observed)
                β = mean(simulated) / mean(observed)
                return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
            end
    
            function kge2(df::DataFrame)
                observed, simulated = df[:,5],df[:,6]
                r = cor(observed, simulated)
                α = std(simulated) / std(observed)
                β = mean(simulated) / mean(observed)
                return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
            end
    
            function nse(predictions::Vector{Float64}, targets::Vector{Float64})
                return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
            end
    
            function nse(df::DataFrame)
                observed, simulated = df[:,5],df[:,6]
                return (1 - (sum((simulated .- observed).^2) / sum((simulated .- mean(observed)).^2)))
            end
    
            # function nse(df::DataFrame;kw...)
            #     observed, simulated = df[:,5],df[:,6]
            #     return (1 - (sum((simulated .- observed).^2) / sum((simulated .- mean(observed)).^2)))
            #     if more
            #         return(getfield(df[:,5:6],:colindex))
            #     end
            #     #more=(names(df[:,5:6]))
            #     #more=propertynames(df[:,5:6])
            # end
    
            # function nse(df::DataFrame;more::AbstractString)
            #     observed, simulated = df[:,5],df[:,6]
            #     nse=return (1 - (sum((simulated .- observed).^2) / sum((simulated .- mean(observed)).^2)))
            #     print("NSE $getfield(df[:,5:6],:colindex) of is $nse")
            #     #more=(names(df[:,5:6]))
            #     #more=propertynames(df[:,5:6])
            # end
    
    
            function kge_read(path::AbstractString, ext::AbstractString)
                # function kge2(observed, simulated)
                #     r = cor(observed, simulated)
                #     α = std(simulated) / std(observed)
                #     β = mean(simulated) / mean(observed)
                #     return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
                # end
                files = readdir(path)
                for file in files
                    file_path = joinpath(path, file)
                    if isfile(file_path) && endswith(file, ext)
                        dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
                        observed  = dd[:,5]
                        simulated = dd[:,6]
                        kge_value = kge2(observed, simulated)
                        println(replace("KGE value is $kge_value on $file_path", "\\"  => "/"))
                    elseif isdir(file_path)
                        dfs_in_subdir = kge_read(file_path, ext)
                    end
                end
            end
    
            function qgg()
                kge_read(pwd(),"out");
            end
    
            function kge_read()
                kge_read(pwd(),"out");
            end
    
            # function kge_read(ext::AbstractString)
            #     path = pwd()
            #     files = readdir(path)
            #     v=[]
            #     for file in files
            #         file_path = joinpath(path, file)
            #         if isfile(file_path) && endswith(file, ext)
            #             dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
            #             observed  = dd[:,5]
            #             simulated = dd[:,6]
            #             kge_value = kge2(observed, simulated)
            #             nm = basename(file_path)
            #             println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
            #             push!(v,Dict(nm=>kge_value))
            #             #push!(v,(kge_value,nm))
            #         elseif isdir(file_path)
            #             dfs_in_subdir = kge_read(file_path, ext)
            #         end
            #     end
            #     return(v)
            # end
    
            function kge_read(ext::AbstractString)
                path = pwd()
                files = readdir(path)
                v=[]
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
                        println(replace("NSE value is $nse_value on $nm", "\\"  => "/"))
                        push!(v,Dict(nm=>kge_value))
                        #push!(v,(kge_value,nm))
                    elseif isdir(file_path)
                        dfs_in_subdir = kge_read(file_path, ext)
                    end
                end
                return(v)
            end
    
    
            #https://stackoverflow.com/questions/42499528/julia-convention-for-optional-arguments
            #. Mostly you'll define two methods for your function:
            function rp(file::AbstractString, lyr::Int)
                #xlyr = length(lyr)!=1 ? 1 : lyr
                ts=read(Raster(file,missingval=0))
                x = ts[t=lyr]
                contourf(x; c=cgrad(:thermal),size=(1200, 800))
            end
    
            function rp(file::AbstractString)
                ts=read(Raster(file,missingval=0))
                contourf(ts; c=cgrad(:thermal),size=(1200, 800))
            end
    
    
            function rp(file::AbstractString, lyr::Int)
                #xlyr = length(lyr)!=1 ? 1 : lyr
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
    
            function rpall(file::AbstractString)
            	xr = read(Raster(file;crs=EPSG(25832),missingval=0))
            	#xr = read(Raster(file;crs=EPSG(25832),missingval=-9999))
            	Plots.plot(xr;c=cgrad(:thermal),
                xlabel="",
                ylabel="",
                size=(1200*.8, 800*.8))
            end
    
            function cpl(file::AbstractString,lyr::Int,msk::Int)
                x=read(Raster(file,missingval=0)) #-9999
                x=x[t=lyr]
                zm = x .< msk
                Plots.contourf(Rasters.mask(x; with=zm); 
                c=cgrad(:thermal),
                xlabel="",
                ylabel="",
                size=(1200, 800))
            end
    
    
            # function crpl(file::AbstractString,lyr::Int,vmin::Int,vmax::Int)
            #     x=read(Raster(file,missingval=-9999)) #
            #     x=x[t=lyr]
            #     tmp=(x .> vmin);
            #     mm=Rasters.mask(x; with=tmp);
            #     zm =(x .< vmax);
            #     Plots.contourf(Rasters.mask(mm; with=zm); c=cgrad(:thermal),size=(1200, 800),
            #     xlabel="",ylabel="")
            # end
    
            # function crpl(file::AbstractString,lyr::Int,vmin::Int,vmax::Int)
            #     x=read(Raster(file,missingval=-9999)) #
            #     x=x[t=lyr]
            #     x=Rasters.rebuild(x;missingval=vmin)
            #     #tmp=(x .> vmin);
            #     #mm=Rasters.mask(x; with=tmp);
            #     zm =(x .< vmax);
            #     Plots.contourf(Rasters.mask(x; with=zm); c=cgrad(:thermal),size=(1200, 800),
            #     xlabel="",ylabel="")
            # end
    
    
            function cntplt(file::AbstractString)
                x=read(Raster(file,missingval=0))
                #describe(x)
                Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
            end
    
            # getf("sum")[2] |> cntplt
            # v = getf("sum")[3]
            # cntplt(v)
            # rv = Raster(v)
            # cntplt(rv)
            # #describe(ar[2])  
    
            function readras(file::AbstractString)
                x=read(Raster(file,missingval=0)) #read all in RAM
                describe(x)
                return(x)
            end
    
            function readras(path::Regex)
                "reads first match"
                v::Vector{String} = readdir();
                v = v[broadcast(x->endswith(x,"nc"),v)];
                file = v[(broadcast(x->occursin(path,x),v))] |>first;
                x::Raster = read(Raster(file,missingval=0))
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
    
    
            function cntplt(file::Raster{Union{Missing, Float64}, 2, Tuple{Y{Mapped{Float64, Vector{Float64}, DimensionalData.Dimensions.LookupArrays.ReverseOrdered, DimensionalData.Dimensions.LookupArrays.Regular{Float64}, DimensionalData.Dimensions.LookupArrays.Points, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}, EPSG, EPSG, Y{Colon}}}, X{Mapped{Float64, Vector{Float64}, DimensionalData.Dimensions.LookupArrays.ForwardOrdered, DimensionalData.Dimensions.LookupArrays.Regular{Float64}, DimensionalData.Dimensions.LookupArrays.Points, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}, EPSG, EPSG, X{Colon}}}}, Tuple{Dim{:t, DimensionalData.Dimensions.LookupArrays.Sampled{Float64, Vector{Float64}, DimensionalData.Dimensions.LookupArrays.ForwardOrdered, DimensionalData.Dimensions.LookupArrays.Regular{Float64}, DimensionalData.Dimensions.LookupArrays.Points, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}}}}, Matrix{Union{Missing, Float64}}, Symbol, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}, Missing})
                x=file
                describe(x)
                Plots.contourf(x; c=cgrad(:thermal))
                #,size=(1200, 800))
            end
    
            function cplt(file::AbstractString)
                x=read(Raster(file))
                Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
            end
    
            # v=ct("ev")
            # v=ct("gws")
            # lst=[]
            # nconly=v[broadcast(x->endswith(x,".nc"),v)]
            # for i in nconly;push!(lst,cplt(i));end   
            # #savefig([plot,] filename)     
            # #for i in nconly; #println(nconly[i]); 
            # for i in 1:length(nconly);fn=nconly[i];outname=replace(fn,"nc"=>"jl.png");
            # savefig(lst[i],outname);println(outname,"saved!");end
    
    
            # function xx(ext::AbstractString)
                # path = pwd()
                # v = []
                # files = readdir(path)
            	# for file in files
                       # i = joinpath(path, file)
                       # if isfile(i) &&  occursin(Regex(ext),file) && endswith(file, ".nc")
            	   # #println(i);
            	   # i=basename(i)
            	   # osize = stat(i).size
            	   # @printf("%-40s %15.2f MB\n","$(i):",osize/1024^2);
            	   # #outname=replace(i,"nc"=>"jl.png");
                       # #println(outname," saved!");
            	   # push!(v,i)
                    # end
            	# return(v)
               # end
            # end	
    
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
    
            # flags = Dict(
                # :s_srs => "epsg:25832",
                # :t_srs => "epsg:4326",
                # :tr => [100,100],
                # :r => :near,
            # )
            # flags = Dict(
                # :tr => [100,100],
                # :r => :near,
            # )
            # warp(r[t=3],flags)  |> Plots.plot
    
            function cc(ext::AbstractString)
                path = pwd()
                files = readdir(path)
            	for file in files
                       i = joinpath(path, file)
                       if isfile(i) && occursin(Regex(ext),file) && endswith(file, ".nc")
            	   outname=replace(i,"nc"=>"jl.png");
                       #println(outname)
                       r=read(Raster(i,missingval=-9999,mappedcrs=EPSG(25832)));
                       #p=Plots.contourf(r;
                       p=Plots.plot(r;
            		title=replace(basename(i),".nc"=>""), #split(outname,"/")[end], #basename(i)
            		c=cgrad(:thermal),
            		size=(1200, 800));
                       savefig(p,outname)
                       println(basename(outname)," saved!");
                    end
               end
            end
    
            	# nconly=v[broadcast(x->endswith(x,".nc"),v)]
            	# for i in nconly;push!(lst,
            		# Plots.contourf(read(Raster(i,missingval=0)); c=cgrad(:thermal),size=(1200, 800)));
            	# end 
            	# for i in 1:length(nconly);
            		# fn=nconly[i];
            		# outname=replace(fn,"nc"=>"jl.png");
            		# savefig(lst[i],outname);println(outname,"saved!");
            	# end   
            # end
            function fd()
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
    
            function fd(cwd::AbstractString)  
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
    
            # dfs=loadalldfs(p)
    
            # T = typeof(dfs[2])
            # for (name, typ) in zip(fieldnames(T), T.types)
            #     println("type of the fieldname $name is $typ")
            # end
    
            # z=map(x->DataFrames.metadata(x)|>collect,dfs)
            # z=map(x->basename(x[1][2]),z)
            # findall("al",z[2])
            # map(x->occursin(r"wol",x),z)
            # map(x->findall(r"wol",x),z)
            # #function grepl(df::DataFrame)
            # using Grep
            # grep("wlf",z)
            # grep("wlf",z)|>getindex
    
            #dfs[grep("wlf",z)]
            #plotf()
    
    
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
                o = DataFrames.metadata(df)|>collect
                ti = basename(o[1][2])
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
    
            #collect(DataFrames.metadata(df))[1][2]
            #for i in dfs;collect(DataFrames.metadata(i))[1][2]|>basename|>println;end
    
    
            function aplot(df::DataFrame)
                df[!,:year]=year.(df[!,:date]) ;
                s = Symbol.(filter(x->!occursin("date",x),names(df)))
                o = DataFrames.metadata(df)|>collect
                ti = "AndrewsPlot of "*basename(o[1][2])
                @df df andrewsplot(:year, cols(s), legend = :topleft,title=ti)
            end
    
    
            function getnames(dfs::Vector)
            	nms = [];
            	for i in dfs;	
            		x=collect(DataFrames.metadata(i))[1][2]|>basename
            		push!(nms, x)
            	end
            	return(nms)
            end
            function getnames(dfs::DataFrame)
            	x=collect(DataFrames.metadata(dfs))[1][2]|>basename
            	return(x)
            end
    
    
            function loadalldfs(path::Vector{Any})
                files = path
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
    
            function loadalldfs(path::Vector{String})
                files = path
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
    
            function loadalldfs(path::Regex)
                v = readdir();
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
    
            readall = loadalldfs
    
    
            function loadalldfs(path::AbstractString)
                files = readdir(path)
                dfs = DataFrame[]
                #nms = []
                for file in files #&& occursin(Regex(prefix),file)
                    if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",file))
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
    
    
            #kge_read(pwd(),"out")
            #println("try ",kge_read(pwd(),"out"))
    
            # md = loadso(p,"so")
            # md[end-1]|>vio
    
            # z = loadalldfs(p)
            # vv = listdfs(p)
            # px = z[end-1]|>vio;
            # plot!(title=vv[end-1])
    
            # function xd()
            # rootdir=pwd()
            # for (looproot, dirs, filenames) in walkdir(rootdir)
                    # for dir in dirs
                        # if isdir(dir)
                        # @printf("%-8s\t|","$dir");
                    # end
                 # end
              # end
            # end
    
            #REGEX AND!!!  
            #z = v[(broadcast(x->occursin(r".*(th)+.*(nc)+.*",x),v))] 
            ##geht gut.
    
            #if (occursin(Regex(prefix,"i"),filename))
            function regand(v::Vector{String},x1::AbstractString,y1::AbstractString)
                needle=join([x1,y1],"+.*");
                z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))] 
            return(z)
            end
    
            # regand(v,"sum","nc\\E")
            # x1,y1="sum","nc\\E"
            # needle=join([x1,y1],"+.*");
            # Regex(needle,"i")
            #xv=("utm_rcm", "rcm-c4")    
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
    
    
            function regand(v::Vector{String},xv::Regex)
                "here you can put any regex to filter the Vector"
                z = v[(broadcast(x->occursin(xv,x),v))] 
            return(z)
            end
    
    
    
            function nconly(x1::AbstractString)
            v::Vector{String} = readdir();
            v = v[broadcast(x->endswith(x,"nc"),v)];
            z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
            return(z)
            end
    
            #tats. map is same here
            #v[map(x->endswith(x,"nc"),v)]  
            #https://stackoverflow.com/questions/52892726/julia-whats-the-difference-between-map-and-broadcast
    
            function nconly(Any)
            v = readdir();
            z = v[broadcast(x->endswith(x,"nc"),v)];
            return(z)
            end
    
    
    
            function readallras(path::AbstractString)
                v = readdir(path);
                v = v[broadcast(x->endswith(x,"nc"),v)];
                z=[];
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
                z=[];
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
                z=[];
                for s in v; 
                #if contains(x1,s) & occursin(r"nc$",s)
                ts=read(Raster(s,missingval=0))
                push!(z,ts);
                end
                return(z)
            end
    
    
    
            function hometeo()
                cd("/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/")
                println("you are here: ",pwd())
            end
    
            function homreg()
                cd("/mnt/d/Wasim/regio/out/");
                println("you are here: ",pwd())
                fd()
            end
    
    
            function writewa(file::AbstractString, df::DataFrame)
                dout = df
                dout.YY = map(x ->year(x),dout.date)
                dout.MM = map(x ->month(x),dout.date)
                dout.DD = map(x ->day(x),dout.date)
                dout[!, "HH"] .= 24
                #df = select!(df,Symbol.(filter(x->!occursin("date",x), names(df))))
                #dout = select(df, Not(:date))
                #dout = dout[!,Cols([:YY,:MM,:HH,:DD],1:end-4)]
                dout = dout[!,Cols([:YY,:MM,:HH,:DD],Not(Cols(r"date")))]
                #cls = propertynames(df)|>sort|>reverse
                #df = df[!,cls[2:end]] 
                CSV.write(file, dout, transform = (col, val) -> something(val, missing),delim="\t")  
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
    
    
            function vgg(regex::AbstractString, ending::AbstractString)
                cmd = `grep --color=always -C2 -rIHn -E "$regex" --include="*.$ending"`
                run(cmd)
            end
    
            ##julia with no regex
    
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
    
            #faster and pure julia:
            function vgrep(regex, file_ending)
        # list files that start with "qgko" and end with file_ending
        #files = filter(file -> startswith(file, "qgko") && endswith(file, file_ending), readdir())
        files = filter(file -> endswith(file, file_ending), readdir())
        # loop over each file
        for file in files
            # open the file and read each line
            #xlines = readlines(file)
            #filter(z -> true, xlines) |> (x -> for i in 1:length(x) println("$i\t$(x[i])") end)
    
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    # check if the line matches the regex
                    if occursin(Regex(regex,"i"), line)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                        # print the file name, line number and line content
                        #println("$file:$(f.lineno):$line") <-nope
                        #m=match(regex, line)
                        #m=count(_ -> true, line) #das zählt die linechars
                        #println("$file: $counter:\t $line")
                    end
                end
            end
        end
            end
    
            #cd("/mnt/d/Wasim/regio/out/");
            #cd("/mnt/d/Wasim/streu/out");
    
            #using GeoArrays
            function gplot(r::AbstractString)
                geoarray = GeoArrays.GeoArray(ArchGDAL.readraster(lk))
                geoarray|>plot
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
    
            function fdf(df::DataFrame)
                nrows=size(df)[2]-1 
                o = DataFrames.metadata(df)|>collect
                ti = basename(o[1][2])
                fig = make_subplots(
                    shared_xaxes=true, 
                    shared_yaxes=true    
                    );
                for i in 1:nrows;
                    add_trace!(fig, 
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
                println(describe(df))
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
                    end
                end
                return(s)
            end
    
            #fdd()
    
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
    
            #dfilter("cl",dfs)
            #typeof(dfs)
            #dfs
            #regex="qout"
            #z::Vector{Raster} = 
            function filterplot(regex::AbstractString,z::Vector{Raster})
                "selects first match and plots..."
                rr::Raster = z[map(n->occursin(Regex(regex,"i"),n),
                               #map(x->map(String,name(x)),z)
                               map(String,map(name,z))
                               )                   
                               ] |> first
                #Plots.contourf(rr; c=cgrad(:thermal),size=(1200*.7, 800*.7))
                ti::Symbol = name(rr)
                fct::AbstractFloat = 0.5
                Plots.plot(rr; 
                    c=cgrad(:thermal),
                    title=ti,
                    size=(1200*fct, 800*fct)
                    )
            end
    
            #filterplot("tem",rds)
            #filterplot("vap",rds)
    
            function filterplot(regex::AbstractString,dfs::Vector{DataFrame})
                "selects first match and plots..."
                df::DataFrame = dfs[map(n->occursin(Regex(regex,"i"),n),
                map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
                )] |> first
                dfp(df)
                #indexin(1:length(dfs),
                #map(x->basename(only(DataFrames.metadata(x))[2]),dfs))
            end
    
            function filterplot!(regex::AbstractString,dfs::Vector{DataFrame})
                "selects first match and add to plot..."
                df::DataFrame = dfs[map(n->occursin(Regex(regex,"i"),n),
                map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
                )] |> first
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                #Symbol(names(df))
                #s = propertynames(df)[Not(end)] #geht auch, aber positionsabhängig
                @df df Plots.plot!(:date,cols(s),legend = :topright)
            end
    
            # filterplot("win",dfs)
            # filterplot!("qg",dfs)
            # filterplot!("qout",dfs)
            #dfs[5] |>dfp
            #convert(DataFrame,dfs[5]|>DataFrames.metadata|>only)
    
            function fread(ext::AbstractString)
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
    
            # function dfpall(dfs::Vector{DataFrame})
            #     for xx in 2:length(dfs)
            #         dfp(first(dfs))
            #         y = filter(x->!occursin("date",x),
            #         names(dfs[xx]))
            #         s = map(y -> Symbol(y),y)
            #         px=@df dfs[xx] Plots.plot!(:date,cols(s),legend = :topright)
            #         return(px)
            #     end
            # end
    
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
    
    
            # function dfpall(dfs::Vector{DataFrame})
            #     odf=DataFrame[]
            #     odf = map(x->innerjoin(dfs[1],x,on="date",makeunique=true),dfs)
    
            #     # for i in dfs; 
            #     #     push!(odf,innerjoin(dfs[1],i,on="date",makeunique=true))
            #     # end 
            #     out = reduce(hcat, odf) 
            # end
    
    
            # z=fread("eva")
            # dfpall(z)
            #filterplot("ev",z)
    
            # dfs = z
            # for df in 2:length(dfs)
            #     y = filter(x->!occursin("date",x),
            #     names(dfs[df]))
            #     s = map(y -> Symbol(y),y)
            #     println(s,df,dfs[df])
            # end
    
    
            #xdf(df)
            #fdf(df)
    
    
            function p()
                return(pwd())
            end
    
            function yrsum(x::String)
                df = readdf(x)
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                #ti=DataFrames.metadata(df)|>only|>last|>basename
                df[!, :year] = year.(df[!,:date]);
                df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
                return(df_yearsum)
            end
    
            function yrsum(x::DataFrame)
                    df = x
                    y = filter(x->!occursin("date",x),names(df))
                    s = map(y -> Symbol(y),y)
                    #ti=DataFrames.metadata(df)|>only|>last|>basename
                    df[!, :year] = year.(df[!,:date]);
                    df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
                    return(df_yearsum)
            end
    
            function yrmean(x::String)
                    df = readdf(x)
                    y = filter(x->!occursin("date",x),names(df))
                    s = map(y -> Symbol(y),y)
                    #ti=DataFrames.metadata(df)|>only|>last|>basename
                    df[!, :year] = year.(df[!,:date]);
                    df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
                    return(df_yearsum)
            end
    
            function yrmean(x::DataFrame)
                df = x
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                #ti=DataFrames.metadata(df)|>only|>last|>basename
                df[!, :year] = year.(df[!,:date]);
                df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
                return(df_yearsum)
            end
    
            function bardf(x::String)
                "with String"
                df = readdf(x)
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
                ti=DataFrames.metadata(df)|>only|>last|>basename
                df[!, :year] = year.(df[!,:date]);
                df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
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
                ti=DataFrames.metadata(df)|>only|>last|>basename
                df[!, :year] = year.(df[!,:date]);
                df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
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
            ti=DataFrames.metadata(df)|>only|>last|>basename
            df[!, :year] = year.(df[!,:date]);
            df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
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
        ti=DataFrames.metadata(df)|>only|>last|>basename
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
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
                    ti=DataFrames.metadata(df)|>only|>last|>basename
                    df[!, :year] = year.(df[!,:date]);
                    df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
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
                    ti=DataFrames.metadata(df)|>only|>last|>basename
                    df[!, :year] = year.(df[!,:date]);
                    df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
                    @df df_yearsum Plots.plot(:year,
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
    
            ##https://chifi.dev/weve-been-writing-julia-wrong-speed-up-julia-with-annotations-9144845c3e24
            # route_from_dir(".")
            # ll()
    
            mutable struct LinearRegression{A<:AbstractFloat, B<:AbstractFloat, P<:Function}
                a::A
                b::B
                predict::P
                function LinearRegression(x::Array,y::Array)
                    # a = ((∑y)(∑x^2)-(∑x)(∑xy)) / (n(∑x^2) - (∑x)^2)
                    # b = (x(∑xy) - (∑x)(∑y)) / n(∑x^2) - (∑x)^2
                    if length(x) != length(y)
                        throw(ArgumentError("The array shape does not match!"))
                    end
                    Σx::Float64 = sum(x)
                    Σy::Float64 = sum(y)
                    xy::Array = x .* y
                    Σxy::Float64 = sum(xy)
                    x2::Array{Float64} = x .^ 2
                    Σx2::Float64 = sum(x2)
                    n::Int64 = length(x)
                    # Calculate a
                    a::Float64 = (((Σy) * (Σx2)) - ((Σx * (Σxy)))) / ((n * (Σx2))-(Σx^2))
                    b::Float64 = ((n*(Σxy)) - (Σx * Σy)) / ((n * (Σx2)) - (Σx ^ 2))
                    predict(xt::Array) = (xt = [i = a + (b * i) for i in xt]::Array)
                    return new{Float64, Float64, Function}(a::Float64, b::Float64, predict::Function)
                end
            end
    
            # x = randn(50000000)
            # y = randn(50000000)
            # @time LinearRegression(x, y).predict(y)
    
    
            function dfyrs(df::DataFrame;)
                ti = DataFrames.metadata(df)|>only|>last|>basename
                fact,logy = 1,0
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                df[!, :year] = year.(df[!,:date]);
                df = combine(groupby(df, :year), y .=> sum .=> y);
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
                ti = DataFrames.metadata(df)|>only|>last|>basename
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
    
    
    
            function monsum(x::String)
                df = readdf(x)
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                df[!, :month] = month.(df[!,:date]);
                df_monthsum = combine(groupby(df, :month), y .=> sum .=> y);
                return(df_monthsum)
            end
    
            function monsum(x::DataFrame)
                df = x
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                df[!, :month] = month.(df[!,:date]);
                df_monthsum = combine(groupby(df, :month), y .=> sum .=> y);
                return(df_monthsum)
            end
    
            function monmean(x::String)
                df = readdf(x)
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                df[!, :month] = month.(df[!,:date]);
                df_monthsum = combine(groupby(df, :month), y .=> mean .=> y);
                return(df_monthsum)
            end
    
            function monmean(x::DataFrame)
                df = x
                y = filter(x->!occursin("date",x),names(df))
                s = map(y -> Symbol(y),y)
                df[!, :month] = month.(df[!,:date]);
                df_monthsum = combine(groupby(df, :month), y .=> mean .=> y);
                return(df_monthsum)
            end
    
            function barp(x::DataFrame)
                "with DataFrame input"
                    df = x
                    ti=DataFrames.metadata(df)|>only|>last|>basename
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
                df=df[:,Cols(4,end)] #only 1stlevel and date col.
                DataFrames.metadata!(df, "filename", x, style=:note);
            end
        end
    end
    
    println("you are here: ",p())
    fd(pwd())
    #set backend from gr() to 
    #plotlyjs()

end

# #using Printf
# using Printf
# module wasim
#     # include("src/WaSiM.jl")
#     # end
    
#     using SnoopPrecompile    # this is a small dependency
#     import Printf, DataFrames, CSV, Statistics, Dates, Distributions
    
#     @precompile_setup begin
#         # Putting some things in `setup` can reduce the size of the
#         # precompile file and potentially make loading faster.
#         #list = [OtherType("hello"), OtherType("world!")]
#         @precompile_all_calls begin
#             # all calls in this block will be precompiled, regardless of whether
#             # they belong to your package or not (on Julia 1.8 and higher)
#             # d = Dict(MyType(1) => list)
#             # x = get(d, MyType(2), nothing)
#             # last(d[MyType(1)])
    
#             #using Query
    
#             function setup()
#                 #include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/functions.jl")
#                 include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/fun-precomp.jl")
#             end
    
    
#             function recursive_glob_prfx(rootdir=".", prefix="")
#                 results = []
#                 for (looproot, dirs, filenames) in walkdir(rootdir)
#                     for filename in filenames
#                         if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                             push!(results, joinpath(looproot, filename)) 
#                         end
#                     end
#                 end
#                 return results
#             end
    
#             function qgk(rootdir=".", prefix="")
#                 files = []
#                 for (looproot, dirs, filenames) in walkdir(rootdir)
#                     for filename in filenames
#                         if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                             push!(files, joinpath(looproot, filename)) 
#                         end
#                     end
#                 end
    
#                 for z in files
#                     println(raw"file:	",basename(z),"...")
#                     m = filter(line -> occursin(r"^[LIN. R]|^[LOG. R]|^CO",line), readlines(open(z)))
#                     for l in m
#                         x = replace(l, r"\s+" => "\t")
#                         x = replace(x, ".\t" => " ")
#                         println(x)
#                     end
#                 end
#                 return nothing
#             end
    
    
#             function penman_monteith(ETo, G, T, Td, u2, es, ea, Ra)
#                 """
#                 Calculates the potential evapotranspiration (PET) using the Penman-Monteith equation.
            
#                 Parameters
#                 ----------
#                 ETo : Float64
#                     Reference evapotranspiration (mm/day).
#                 G : Float64
#                     Soil heat flux density (mm/day).
#                 T : Float64
#                     Air temperature (°C).
#                 Td : Float64
#                     Dew point temperature (°C).
#                 u2 : Float64
#                     Wind speed at 2 m height (m/s).
#                 es : Float64
#                     Saturation vapor pressure (kPa).
#                 ea : Float64
#                     Actual vapor pressure (kPa).
#                 Ra : Float64
#                     Aerodynamic resistance (s/m).
            
#                 Returns
#                 -------
#                 PET : Float64
#                     Potential evapotranspiration (mm/day).
#                 """
#                 # Constants
#                 R = 8.314 # J/mol/K
#                 cp = 1.013e-3 # kJ/g/K
            
#                 # Latent heat of vaporization (MJ/kg)
#                 Lambda = 2.501 - 0.002361 * T
            
#                 # Psychrometric constant (kPa/°C)
#                 gamma = cp * P / (0.622 * Lambda)
            
#                 # Slope of the saturation vapor pressure curve (kPa/°C)
#                 delta = 4098 * es / (T + 237.3) ^ 2
            
#                 # Net radiation (MJ/m2/day)
#                 Rn = (1 - 0.23) * ETo
            
#                 # Air density (kg/m3)
#                 rho = P * 1000 / (R * (T + 273.15))
            
#                 # Specific heat of air (kJ/kg/K)
#                 cpa = 1.013 * rho ^ -0.0065 * 1000
            
#                 # Delta term (MJ/m2/day/°C)
#                 delta_term = (delta / (delta + gamma)) * (Rn - G)
            
#                 # Psi term (MJ/m2/day)
#                 psi_term = (gamma / (delta + gamma)) * rho * cp * (es - ea) / Ra * u2
            
#                 # Potential evapotranspiration (mm/day)
#                 PET = (delta_term + psi_term) / Lambda
            
#                 return PET
#             end
    

    
#             ##w endswith
#             function lg(path::AbstractString, ext::AbstractString)
#                 files = readdir(path)
#                 v=[]
#                 for file in files
#                     file_path = joinpath(path, file)
#                     if isfile(file_path) && endswith(file, ext)
#                        println(file_path)
#             	   push!(v,file_path)
#                     end
#                 end
#                 return(v)
#             end
    
#             #list(x) = Any[i for i ∈ x]
    
#             function ddense(path::String,skip::Int,start::Int,stop::Int)
#                 ms=["-999","-9999","lin","log","LIN","LOG"]
#                 df = CSV.read(path,DataFrame,skipto=skip,
#                 missingstring=ms,delim="\t",comment="-",silencewarnings=false,
#                 ntasks=4,downcast=true,normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
#                 df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#                 df=df[:,Not(1:3)]
#                 #nrows=size(df)[2]-1
#                 println(propertynames(df))
#                 #@df df density(:_11, group = (:tot_average, :date), legend = :topleft)
#                 #@df df density(:tot_average, legend = :topleft)
#                 @df df density(cols(start:stop), legend = :topleft)
#             end
    
#             function denselog(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match and plots..."
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                      map(x->basename(only(DataFrames.metadata(x))[2]),
#                      dfs))] |> first
#                      s = propertynames(df)[Not(end)];
#                      o = DataFrames.metadata(df)|>collect
#                      ti = basename(o[1][2])
#                      @df df density(
#                         cols(s),
#                         title=ti,
#                         yaxis=:log,
#                         legend = :topright) 
#             end
    
    
#             function denseplot(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match and plots..."
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                      map(x->basename(only(DataFrames.metadata(x))[2]),
#                      dfs))] |> first
#                      s = propertynames(df)[Not(end)];
#                      o = DataFrames.metadata(df)|>collect
#                      ti = basename(o[1][2])
#                      @df df density(cols(s),title=ti,legend = :topright) 
#             end
    
    
#             function denseplot(df::DataFrame)
#                 #println(propertynames(df))
#                 s = propertynames(df)[Not(end)] #masks last column == date     #[1:end-1]
#                 #,propertynames(df)[end]
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 @df df density(cols(s), legend = :topright, title=ti)
#             end
    
#             function denseplot(df::String)
#                 df=readdf(df)
#                 s = propertynames(df)[Not(end)]
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 @df df density(cols(s), legend = :topright, title=ti)
#             end
    
#             function dfp(df::DataFrame)
#         o = DataFrames.metadata(df)|>collect
#         ti = basename(o[1][2])
#         if (any(x->occursin("year",x),names(df)))
#             s = Symbol.(filter(x->!occursin("year",x),names(df)))
#             @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
#         else    
#         s = Symbol.(filter(x->!occursin("date",x),names(df)))
#         @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
#         end
#             end
    
#             function dfp(df::String)
#                 df=readdf(df)
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 if (any(x->occursin("year",x),names(df)))
#                     s = Symbol.(filter(x->!occursin("year",x),names(df)))
#                     @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
#                 else    
#                 s = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
#                 end
#             end
    
#             function dfp(mm::Regex)
#                 """
#                 plots wasim timeseries
#                 """
#                 df=readdf(mm)
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 if (any(x->occursin("year",x),names(df)))
#                     s = Symbol.(filter(x->!occursin("year",x),names(df)))
#                     @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
#                 else    
#                 s = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
#                 end
#             end
    
    
#             function dfp(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match and plots..."
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                      map(x->basename(only(DataFrames.metadata(x))[2]),
#                      dfs))] |> first
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 if (any(x->occursin("year",x),names(df)))
#                     s = Symbol.(filter(x->!occursin("year",x),names(df)))
#                     @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
#                 else    
#                 s = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
#                 end
#             end
    
#             function getf(ext::AbstractString)
#                 cwd = pwd() 
#                 m = []
#                 for (root, dirs, files) in walkdir(cwd)
#                  for file in files
#                  if isfile(file) && occursin(Regex(ext),file)
#             	 nm=joinpath(root, file)
#             	 push!(m,(nm))
#                  end
#                 end 
#                 end 
#                  return(m)
#             end 
    
#             function getdf(ext::AbstractString)
#                 cwd = pwd() 
#                 m = []
#                 for (root, dirs, files) in walkdir(cwd)
#                  for file in files
#                  if isfile(file) && occursin(Regex(ext),file)&&(!occursin(r"txt|yrly|nc|png|svg",file))
#             	 nm=joinpath(root, file)
#             	 push!(m,(nm))
#                  end
#                 end 
#                 end 
#                  return(m)
#             end 
    
    
#             # getf(".*(^th)+.*(nc)+.*")  
#             # #SAME
#             # getf("^th+.*nc")
#             # ###lookbehind	
#             # #getf("stack?+.*nc") 
#             # #getf("!stack?+.*nc") 
    
#             function plotf(ext::AbstractString)
#                 cwd = pwd() 
#                 m = []
#                 for (root, dirs, files) in walkdir(cwd)
#                  for file in files
#                  if isfile(file) && occursin(Regex(ext),file)&&(!occursin(r"txt|yrly|nc|png|svg",file))
#             	 nm=joinpath(root, file)
#             	 push!(m,(nm))
#                  end
#                 end 
#                 end 
#                  return(
#                  dfp(readdf(m[1])))
#             end 
    
#             function plotf(ext::String)
#                 dfp(readdf(ext))
#                 plot!(title=basename(ext))
#             end 
    
#             function plotf(ext::DataFrame)
#             dfp(ext)
#             end 
    
#             # plotf("tem")
    
#             function readdf(x::Regex)
#                 """
#                 readdf(x::Regex)
#                 reads first match of regex wasim timeseries
#                 """
#                 rootdir="."
#                 results = []
#                 for (looproot, dirs, filenames) in walkdir(rootdir)
#                     for filename in filenames
#                         if (occursin(x,filename)) && (!occursin(r"txt|yrly|nc|png|svg|grd",filename))
#                             push!(results, joinpath(looproot, filename)) 
#                         end
#                     end
#                 end
#                 x = first(results)
#                 ms=["-9999","lin","log"]
#                 df = CSV.read(x,DataFrame,missingstring=ms,
#                 ntasks=4,
#                 limit=typemax(Int),
#                 types = Float64,
#                 delim="\t",
#                 silencewarnings=true,
#                 normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
#                 df.YY=map(x ->Int(x),df.YY);
#                 df.MM=map(x ->Int(x),df.MM);
#                 df.DD=map(x ->Int(x),df.DD);
#                 df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#                 df=df[:,Not(1:3)]
#                 metadata!(df, "filename", x, style=:note);
#             end
    
    
#             function readdf(x::AbstractString)
#         #ms="-9999"
#         ms=["-9999","lin","log"]
#         df = CSV.read(x,
#         DataFrame,
#         missingstring=ms,
#         #skipto=4, 
#         ntasks=4,
#         limit=typemax(Int),
#         types = Float64,
#         delim="\t",
#         #comment=r"[A-z]",
#         silencewarnings=true,
#         normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
#         df.YY=map(x ->Int(x),df.YY);
#         df.MM=map(x ->Int(x),df.MM);
#         df.DD=map(x ->Int(x),df.DD);
#         df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#         df=df[:,Not(1:3)]
#         metadata!(df, "filename", x, style=:note);
#             end
    
#             readmeteo = readdf
#             loaddf = readdf
#             # function readmeteo(x::AbstractString)
#             #     df = DataFrame(CSV.File(x; missingstring="-9999",
#             #                         skipto=6,
#             #                         limit=typemax(Int),
#             #                         comment="#",
#             #                         stringtype=String))
#             #     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
#             #     df=df[:,Not(1:4)]
#             #     metadata!(df, "filename", x, style=:note);
#             # end
    
#             # function loaddf(path::AbstractString)
#             #     ms=["-999","-9999","lin","log","LIN","LOG"]
#             #     df = CSV.read(path,DataFrame,
#             #     missingstring=ms,
#             #     delim="\t",comment="-",
#             #     silencewarnings=false,
#             #     ntasks=4,downcast=true,	    
#             #     normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
#             #     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#             #     df=df[:,Not(1:3)]
#             #     metadata!(df, "filename", path, style=:note);
#             # end
    
#             # mx=ct("so_a")
#             # x=mx[2]
#             # df = CSV.read(x, DataFrame,missingstring=["-9999"], delim="\t",
#             # skipto=4, #bei nich so: 3
#             #         #comment="[A-z]",
#             #         comment="-",
#             #         silencewarnings=true) |>dropmissing
#             # df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#             # df = select!(df, Not(names(df)[1:4]))
#             # plotf(df)
    
#             function loadso(path::AbstractString, prefix::AbstractString)
#                 files = readdir(path)
#                 dfs = DataFrame[]
#                 for file in files
#                     if isfile(file) && occursin(Regex(prefix),file) && (!occursin(r"txt|yrly|nc|png|svg|ftz_0|ftz",file))
#                        file_path = joinpath(path, file)
#             	   println("reading ",file_path,"...")
#             	   p1 = readdf(file_path)
#             	   push!(dfs, p1)
#                     end
#                 end
#                 return(dfs)
#             end
    
#             #md = loadso(pwd(),"so")
#             #map(names,md)
#             #size(md) 
#             #for i in md; print(size(i)) ; end
#             #broadcast(x -> innerjoin(x,on = :date),md)
#             #innerjoin(md,on = :date)
    
#             # using StatsPlots;
#             # dfs = loadso(pwd(),"te")
#             # vibx(dfs[2])
#             # vibx(dfs[3])
    
#             # size(dfs)
#             # i=dfs[16]
#             # s = propertynames(i)[1:end-1]
#             # @df i plot(:date,cols(s))
    
#             # Plots.gr()
    
    
#             #by(df, :month, broadcast(x->ln[x]=>sum,ln))
    
#             # df[!, :month] = month.(df[!, :date]);
#             # by(df, :month, ln[1]=>sum)
#             # for i in ln;println(by(df, :month, i=>sum));end     
    
#             # u=by(df, :month, ln[1]=>sum);
#             # nm = propertynames(u)[Not(1)]
#             # #str = propertynames(u)[1]  
#             # str = [ @sprintf("%02i", x) for x in u[1] ];
#             # @df u StatsPlots.boxplot(str,cols(nm),fillalpha=0.75, linewidth=0.25);
    
#             #broadcast(x->findall("date",x),names(df))  
#             function toyrsum(df::DataFrame)
#             od = []
#             df[!, :year] = year.(df[!,:date]);
#             ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             for i in ln;
#             #x=(by(df,:year,i=>sum));
#             x = combine(df, :year, AsTable(i) => sum, renamecols=false) 
#             push!(od,x)
#             end ;
#             return(od)
#             end
    
    
#             # #od = DataFrame[]
#             # df[!, :year] = year.(df[!,:date]);
    
#             ##komplett mean:
#             function fullmean(df::DataFrame)
#             df[!, :year] = year.(df[!,:date]);
#             combine(df, :, AsTable(Not([:date,:year])) => mean, renamecols=false)
#             end
    
    
#             function cattoyrsum(df::DataFrame)
#             ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             df[!, :year] = year.(df[!,:date]);
#             it=[]
#             od=(by(df,:year,ln[1]=>sum));
#             #od = combine(df, :year, AsTable(ln[1]) => sum, renamecols=false) 
#             for i in ln[2:end];
#             x=(by(df[Not(:date)],:year,i=>sum));
#             #x = combine(df, :year, AsTable(i) => sum, renamecols=false)
#             push!(it,x[end])      
#             #push!(it,x[!,end])      
#             end ;
#             # x = combine(df, :year, AsTable(:) => sum, renamecols=false)
#             # combine(groupby(df,:year),:=>sum)
#             ot = hcat(od,DataFrame(it))  
#             #ot=join(od,x[end],:year,makeuniuqe=true)
#             #ot = hcat(od,it)
#             return(ot)
#             end
    
#             # md=cattoyrsum(df)
#             # @df md plot(:year,cols(propertynames(md)[2:end]),yaxis=:log)      
    
#             function cattoyrmean(df::DataFrame)
#             ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             df[!, :year] = year.(df[!,:date]);
#             it=[]
#             od=(by(df[Not(:date)],:year,ln[1]=>mean));
#             for i in ln[2:end];
#             x=(by(df,:year,i=>mean));
#             push!(it,x[end])      
#             end ;
#             ot = hcat(od,DataFrame(it))  
#             return(ot)
#             end
    
#             # md=cattoyrmean(df)
#             # @df md plot(:year,cols(propertynames(md)[2:end])) 
#             # function tsyr(df::DataFrame)
#             	# str = [ @sprintf("%02i", x) for x in (year.(df.date)) ];
#             	# #tsn = "date";
#             	# #ln = propertynames(df)[Not(tsn)];
#             	# ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             	# @df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
#             	# @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
#             	# @df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
#             # end
    
    
#             function qpl(df::DataFrame)
#     	ti = DataFrames.metadata(df)|>only|>last|>basename
#         s = names(df)[1:2]
#     	t2 = string.(ti,"\n",s[1],"|",s[2],ti)
#         StatsPlots.plot( 
#     	qqplot(df[!,1],df[!,2], qqline = :fit), 
#     	qqplot(Cauchy,df[!,2]), 
#     	qqnorm(df[!,2], qqline = :R),
#         title = t2)
#             end
    
#             function qpl(x::AbstractString)
#     	df = readdf(x)
#         ti = DataFrames.metadata(df)|>only|>last|>basename
#         s = names(df)[1:2]
#     	t2 = string.(ti,"\n",s[1],"|",s[2],ti)
#         StatsPlots.plot( 
#     	qqplot(df[!,1],df[!,2], qqline = :fit), 
#     	qqplot(Cauchy,df[!,2]), 
#     	qqnorm(df[!,2], qqline = :R),
#         title = t2)
#             end
    
    
#             function vibx(df::DataFrame)
#             	str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
#             	#ln = propertynames(df[end-1])
#             	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             	@df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
#             	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
#             	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
#             end
    
#             function vibx(df::String)
#             	df = readdf(df)
#                 str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
#             	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             	@df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
#             	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
#             	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.5,marker=(:black,stroke(1)),legend=false)
#             end
    
#             function vio(df::DataFrame)
#             	str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
#             	#ln = propertynames(df[end-1])
#             	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             	@df df StatsPlots.violin(str,cols(ln),linewidth=0.1)
#             	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
#             #	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
#             end
    
#             # gr()
#             # default(show = true)
#             # df = readdf(r"qgk")
#             # nms = unique!(map(monthabbr,(month.(df.date))))
#             # str = (map(monthabbr,(month.(df.date))))|>sort
#             # #str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
#             # ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#             # ti=DataFrames.metadata(df)|>only|>last|>basename
#             # #df[!, :year] = year.(df[!,:date]);
#             # @df df StatsPlots.violin(str,cols(ln),
#             #     xlabel="Months",
#             #     #xaxis=nms,
#             #     linewidth=0.1,
#             #     title=ti)
#             # #xlabel!(nms|>only)
#             # @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false)
    
#             # gr()
#             # default(show = true)
#             # df = readdf(r"qgk")
#             # #nms = unique!(map(monthabbr,(month.(df.date))))
#             # y = Symbol.(filter(x->!occursin("date",x),names(df)))
#             # df[!, :month] = month.(df[!,:date]);
#             # #regand(v=names(df),"date","month")
#             # dfm = combine(groupby(df, :month), y .=> sum .=> y);
#             # str = (map(monthabbr,(dfm.month)))
#             # dfm.mab = str
#             # #str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
#             # ln = Symbol.(filter(x->!occursin("month",x),names(dfm)))
#             # ti=DataFrames.metadata(dfm)|>only|>last|>basename
#             # @df dfm StatsPlots.violin(cols(ln[Not(end)]),
#             #     xlabel="Months",
#             #     linewidth=0.1,
#             #     title=ti)
#             # #xaxis!(dfm.mab)
    
#             # #df[!, :year] = year.(df[!,:date]);
#             # @df dfm StatsPlots.violin(str,cols(ln),
#             #     xlabel="Months",
#             #     #xaxis=nms,
#             #     linewidth=0.1,
#             #     title=ti)
#             # #xlabel!(nms|>only)
#             # @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false)
    
    
#             	# t = month.(df.date);
#             	# @df df StatsPlots.violin(  string.(broadcast(x -> Dates.monthname(x),t))), :tot_average, linewidth=0.01,legend=false);
#             	# @df df StatsPlots.boxplot!(string.(broadcast(x -> Dates.monthname(x),t))), :tot_average, fillalpha=0.75, linewidth=0.25,legend=false);
#             	# @df df StatsPlots.dotplot!(string.(broadcast(x -> Dates.monthname(x),t))), :tot_average, fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
#             # end
    
#             function vio(df::String)
#                 df = readdf(df)
#             	str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
#             	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 #ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 ti=DataFrames.metadata(df)|>only|>last|>basename
#             	@df df StatsPlots.violin(str,cols(ln),linewidth=0.1,title=ti)
#             	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.55, linewidth=0.25,legend=false);
#             	#@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.15,legend=false)
#                 #marker=(:black,stroke(1)),legend=false)
#             end
    
#             function vio(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match and plots..."
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                      map(x->basename(only(DataFrames.metadata(x))[2]),
#                      dfs))] |> first
#                 str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
#                 ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 @df df StatsPlots.violin(str,cols(ln),linewidth=0.1)
#                 @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.55, linewidth=0.25,legend=true);
#             end
    
    
#             function ldf(path::AbstractString, prefix::AbstractString)
#                 files = readdir(path)
#                 dfs = DataFrame[]
#             #    outname = []
#                 for file in files
#                     #if isfile(file_path) && endswith(file, ext)
#             	#if isfile(file_path) && (startswith(file, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",file))
#                     if isfile(file) && occursin(Regex(prefix),file)&& (!occursin(r"txt|yrly|nc|png|svg",file))
#                        file_path = joinpath(path, file)
#             	   println("reading",file_path)
#             	   p1 = loaddf(file_path)
#             	   push!(dfs, p1)
#             #	   m = match(r".*[.]",basename(file_path))
#             #	   nm = contains(basename(file_path),".") ? string(m.match,"png") : basename(file_path)*".png"
#             #	   push!(outname,nm)
#                     end
#                 end
#                 return(dfs)
#             #    return(outname)
#             end
    
    
#             # a = true
#             # b = 1
#             # c = 2
#             # a ? b : c # 1
#             #outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"
#             #lyr=2
#             #xlyr = length(lyr)!=1 ? 1 : lyr
    
    
#             function pline(path::AbstractString)
#                 ms=["-999","-9999","lin","log","LIN","LOG"]
#                 df = CSV.read(path,DataFrame,
#                 #missingstring="-9999", #also windows
#                 missingstring=ms,
#                 delim="\t",comment="-",
#                 silencewarnings=false,
#                 ntasks=4,downcast=true, # got unsupported keyword arguments "ntasks", "downcast" @windows                                          
#                 normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
#                 df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#                 df=df[:,Not(1:3)]
#                 # ms=["-9999","lin","log","LIN","LOG","--"] #comment="-",
#                 # #df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",comment="-",ignorerepeated=true,silencewarnings=true,typemap=Dict(Int64=>String))  |> @dropna() |> DataFrame
#                 # df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",ignorerepeated=true,silencewarnings=true,typemap=Dict(String=>Int64))
#                 # df = df[completecases(df), :]
#                 # #df = filter( [2]=> x -> !any(f -> f(x), (ismissing)), df)
#                 # #df = filter( [5]=> x -> isnumeric, df)
#                 # #parse.(Date, df[:,1:4])
#                 # #parse.(Date, string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH")
#                 # df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
#                 # df=df[:,Not(1:4)]
#                 nrows=size(df)[2]-1
#                 st=[]
#                 for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
#                 p = make_subplots(rows=nrows, cols=1, 
#                 shared_xaxes=true, 
#                 shared_yaxes=false,
#                 vertical_spacing=0.05,
#                 #subplot_titles= st;
#                 )
#                 for i in 1:nrows;
#                         add_trace!(p, 
#                         PlotlyJS.scatter(x=df.date, y=df[:,i],
#                         name=st[i]),   row=i,     col=1);
#                 end
#                 #relayout!(p,height=600*2,width=900*2,title_text="Series of "*basename(path))
#                 relayout!(p,height=600*1.5,width=900*1.5,title_text="Series of "*basename(path))
#                 p
#             end
    
#             # nrows=size(df)[2]-1
#             # o = DataFrames.metadata(df)|>collect
#             # ti = basename(o[1][2])
#             # px = make_subplots(rows=1,cols=nrows, 
#             # shared_xaxes=true,
#             # shared_yaxes=true);
#             # for i in 1:nrows;
#             #     add_trace!(px, 
#             #     PlotlyJS.scatter(x=df.date, y=df[:,i],
#             #     name=names(df)[i]),   row=1,     col=i);
#             # end
#             # px
    
#             # function dfplotjs(df::DataFrame)
#             #     nrows=size(df)[2]-1 
#             #     #length(names(df))-1
#             #     o = DataFrames.metadata(df)|>collect
#             #     ti = basename(o[1][2])
#             #     fig = make_subplots(
#             #         shared_xaxes=true, 
#             #         shared_yaxes=true    
#             #     #rows=2, cols=2
#             #         );
#             #     for i in 1:nrows;
#             #         add_trace!(fig, 
#             #         PlotlyJS.scatter(x=df.date, y=df[:,i],
#             #         name=names(df)[i]));
#             #     end
#             #     fact=1;
#             #     PlotlyJS.relayout!(fig,
#             #     height=600*fact,width=900*fact,
#             #     title_text="Series of "*ti)
#             #     display(fig)
#             # end
    
#             # pt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/met0/temp_1970.txt"
#             # pt="/mnt/d/Wasim/Goldbach/revision/fab150/gbfab150/qgesfab150.p1.2018"
#             # df = readmeteo(pt)
    
#             # function dfplotjs(df::DataFrame;log::Bool)
#             #     nrows=size(df)[2]-1 
#             #     #length(names(df))-1
#             #     o = DataFrames.metadata(df)|>collect
#             #     ti = basename(o[1][2])
#             #     fig = make_subplots(
#             #         shared_xaxes=true, 
#             #         shared_yaxes=true    
#             #     #rows=2, cols=2
#             #         );
#             #     for i in 1:nrows;
#             #         add_trace!(fig, 
#             #         PlotlyJS.scatter(x=df.date, y=df[:,i],
#             #         name=names(df)[i]));
#             #     end
#             #     fact=1;
#             #     PlotlyJS.relayout!(fig,yaxis_type="log",
#             #     height=600*fact,width=900*fact,
#             #     title_text="Series of "*ti)
#             #     display(fig)
#             # end
    
#             #Wenn Sie ein Schlüsselwort-Argument nicht zugewiesen haben wollen, können Sie es weglassen oder nothing als Wert verwenden. Zum Beispiel:
#             #function dfplotjs(df::DataFrame;logy::Any,fact::Any)
#             function dfplotjs(df::DataFrame;logy::Bool,fact::Float64)
#                 nrows=size(df)[2]-1 
#                 #length(names(df))-1
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 fig = PlotlyJS.make_subplots(
#                     shared_xaxes=true, 
#                     shared_yaxes=true    
#                 #rows=2, cols=2
#                     );
#                 for i in 1:nrows;
#                     PlotlyJS.add_trace!(fig, 
#                     PlotlyJS.scatter(x=df.date, y=df[:,i],
#                     name=names(df)[i]));
#                 end
#                 fact = isnothing(fact) ? 1 : fact; #nice
#                 logy = isnothing(logy)==true ? logy==false : logy==true;
#                 if logy == true
#                     PlotlyJS.relayout!(fig,yaxis_type="log",
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 #elseif isnothing(log) 
#                 else
#                     PlotlyJS.relayout!(fig,
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 end
#                 display(fig)
#             end
    
#             function dfplotjs(df::AbstractString;logy::Bool,fact::Float64)
#                 df=readmeteo(df)
#                 nrows=size(df)[2]-1
#                 #length(names(df))-1
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 fig = PlotlyJS.make_subplots(
#                     shared_xaxes=true, 
#                     shared_yaxes=true    
#                     );
#                 for i in 1:nrows;
#                     PlotlyJS.add_trace!(fig, 
#                     PlotlyJS.scatter(x=df.date, y=df[:,i],
#                     name=names(df)[i]));
#                 end
#                 fact = isnothing(fact) ? 1 : fact; #nice
#                 logy = isnothing(logy)==true ? logy==false : logy==true;
#                 if logy == true
#                     PlotlyJS.relayout!(fig,yaxis_type="log",
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 else
#                     PlotlyJS.relayout!(fig,
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 end
#                 display(fig)
#             end
    
#             function dfplotjs(filepath::AbstractString)
#                 dfplotjs(filepath;logy=false,fact=1.0)
#             end
    
#             function dflogjs(filepath::AbstractString)
#                 dfplotjs(filepath;logy=true,fact=1.0)
#             end
#             # dfplotjs(df;logy=true,fact=0.66)
#             # dfplotjs(df;logy=false,fact=0.7)
#             # dfplotjs(df;logy=true)
#             # dfplotjs(df)
    
#             #ps=xx("qg")
    
#             # dx=readdf(pt)
#             # dfplotjs(df;logy=false,fact=0.7)
#             function dfplotjs(df::AbstractString;logy::Bool,fact::Float64)
#                 df=readmeteo(df)
#                 nrows=size(df)[2]-1
#                 #length(names(df))-1
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 fig = make_subplots(
#                     shared_xaxes=true, 
#                     shared_yaxes=true    
#                     );
#                 for i in 1:nrows;
#                     add_trace!(fig, 
#                     PlotlyJS.scatter(x=df.date, y=df[:,i],
#                     name=names(df)[i]));
#                 end
#                 fact = isnothing(fact) ? 1 : fact; #nice
#                 logy = isnothing(logy)==true ? logy==false : logy==true;
#                 if logy == true
#                     PlotlyJS.relayout!(fig,yaxis_type="log",
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 else
#                     PlotlyJS.relayout!(fig,
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 end
#                 display(fig)
#             end
    
#             function dfplotjs(filepath::AbstractString)
#                 dfplotjs(filepath;logy=false,fact=1.0)
#             end
    
#             function dflogjs(filepath::AbstractString)
#                 dfplotjs(filepath;logy=true,fact=1.0)
#             end
    
#             # Erstellen Sie einige zufällige Zeitreihendaten
#             # dates = Date(2020):Day(1):Date(2021)
#             # y1 = cumsum(randn(length(dates)))
#             # y2 = cumsum(randn(length(dates)))
#             # y3 = cumsum(randn(length(dates)))
#             # y4 = cumsum(randn(length(dates)))
#             # fig = make_subplots(rows=2, cols=2)
#             # # Erstellen Sie einige Zeitreihendiagramme
#             # trace1 = PlotlyJS.scatter(x=dates, y=y1)
#             # trace2 = PlotlyJS.scatter(x=dates, y=y2)
#             # trace3 = PlotlyJS.scatter(x=dates, y=y3)
#             # trace4 = PlotlyJS.scatter(x=dates, y=y4)
#             # add_trace!(fig, trace1, row=1, col=1)
#             # add_trace!(fig, trace2, row=1, col=2)
#             # add_trace!(fig, trace3, row=2, col=1)
#             # add_trace!(fig, trace4, row=2, col=2)
#             # display(fig)
    
    
    
#             #px=PlotlyJS.scatter(x=df.date, y=df[:,Not(:date)])
#             #PlotlyJS.relayout!(px;"Series of "*ti)
#             #PlotlyJS.relayout!(px,height=600*1.5,width=900*1.5,title_text="Series of "*ti)
#             #px.show()
    
#             #PlotlyJS.plot(PlotlyJS.scatter(x=df.date, y=df[!,6]))
    
#             #df = dataset(DataFrame,"iris")
#             #plot(df,x=:sepal_length,y=:sepal_width,z=:petal_width,color=:species,type="scatter3d",mode="markers")
    
#             function dfplot(df::DataFrame)
#                 nrows=size(df)[2]-1
#                 st=[]
#                 for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
#                 p = PlotlyJS.make_subplots(rows=nrows, cols=1, 
#                 shared_xaxes=true, 
#                 shared_yaxes=false,
#                 vertical_spacing=0.05,
#                 )
#                 for i in 1:nrows;
#                         add_trace!(p, 
#                         PlotlyJS.scatter(x=df.date, y=df[:,i],
#                         name=st[i]),   row=i,     col=1);
#                 end
#                 PlotlyJS.relayout!(p,height=600*1.5,width=900*1.5)
#                 return(p)
#             end
    
    
#             function dfplot(df::AbstractString)
#                 df=readmeteo(df)
#                 nrows=size(df)[2]-1
#                 st=[]
#                 for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
#                 p = PlotlyJS.make_subplots(rows=nrows, cols=1, 
#                 shared_xaxes=true, 
#                 shared_yaxes=false,
#                 vertical_spacing=0.05,
#                 )
#                 for i in 1:nrows;
#                         add_trace!(p, 
#                         PlotlyJS.scatter(x=df.date, y=df[:,i],
#                         name=st[i]),   row=i,     col=1);
#                 end
#                 PlotlyJS.relayout!(p,height=600,width=900)
#                 display(p)
#             end
    
    
#             function kge2(observed::Vector{Float64}, simulated::Vector{Float64})
#                 r = cor(observed, simulated)
#                 α = std(simulated) / std(observed)
#                 β = mean(simulated) / mean(observed)
#                 return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
#             end
    
#             function kge2(df::DataFrame)
#                 observed, simulated = df[:,5],df[:,6]
#                 r = cor(observed, simulated)
#                 α = std(simulated) / std(observed)
#                 β = mean(simulated) / mean(observed)
#                 return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
#             end
    
#             function nse(predictions::Vector{Float64}, targets::Vector{Float64})
#                 return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
#             end
    
#             function nse(df::DataFrame)
#                 observed, simulated = df[:,5],df[:,6]
#                 return (1 - (sum((simulated .- observed).^2) / sum((simulated .- mean(observed)).^2)))
#             end
    
#             # function nse(df::DataFrame;kw...)
#             #     observed, simulated = df[:,5],df[:,6]
#             #     return (1 - (sum((simulated .- observed).^2) / sum((simulated .- mean(observed)).^2)))
#             #     if more
#             #         return(getfield(df[:,5:6],:colindex))
#             #     end
#             #     #more=(names(df[:,5:6]))
#             #     #more=propertynames(df[:,5:6])
#             # end
    
#             # function nse(df::DataFrame;more::AbstractString)
#             #     observed, simulated = df[:,5],df[:,6]
#             #     nse=return (1 - (sum((simulated .- observed).^2) / sum((simulated .- mean(observed)).^2)))
#             #     print("NSE $getfield(df[:,5:6],:colindex) of is $nse")
#             #     #more=(names(df[:,5:6]))
#             #     #more=propertynames(df[:,5:6])
#             # end
    
    
#             function kge_read(path::AbstractString, ext::AbstractString)
#                 # function kge2(observed, simulated)
#                 #     r = cor(observed, simulated)
#                 #     α = std(simulated) / std(observed)
#                 #     β = mean(simulated) / mean(observed)
#                 #     return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
#                 # end
#                 files = readdir(path)
#                 for file in files
#                     file_path = joinpath(path, file)
#                     if isfile(file_path) && endswith(file, ext)
#                         dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
#                         observed  = dd[:,5]
#                         simulated = dd[:,6]
#                         kge_value = kge2(observed, simulated)
#                         println(replace("KGE value is $kge_value on $file_path", "\\"  => "/"))
#                     elseif isdir(file_path)
#                         dfs_in_subdir = kge_read(file_path, ext)
#                     end
#                 end
#             end
    
#             function qgg()
#                 kge_read(pwd(),"out");
#             end
    
#             function kge_read()
#                 kge_read(pwd(),"out");
#             end
    
#             # function kge_read(ext::AbstractString)
#             #     path = pwd()
#             #     files = readdir(path)
#             #     v=[]
#             #     for file in files
#             #         file_path = joinpath(path, file)
#             #         if isfile(file_path) && endswith(file, ext)
#             #             dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
#             #             observed  = dd[:,5]
#             #             simulated = dd[:,6]
#             #             kge_value = kge2(observed, simulated)
#             #             nm = basename(file_path)
#             #             println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
#             #             push!(v,Dict(nm=>kge_value))
#             #             #push!(v,(kge_value,nm))
#             #         elseif isdir(file_path)
#             #             dfs_in_subdir = kge_read(file_path, ext)
#             #         end
#             #     end
#             #     return(v)
#             # end
    
#             function kge_read(ext::AbstractString)
#                 path = pwd()
#                 files = readdir(path)
#                 v=[]
#                 for file in files
#                     file_path = joinpath(path, file)
#                     if isfile(file_path) && endswith(file, ext)
#                         dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
#                         observed  = dd[:,5]
#                         simulated = dd[:,6]
#                         kge_value = kge2(observed, simulated)
#                         nse_value = nse(observed, simulated)
#                         nm = basename(file_path)
#                         println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
#                         println(replace("NSE value is $nse_value on $nm", "\\"  => "/"))
#                         push!(v,Dict(nm=>kge_value))
#                         #push!(v,(kge_value,nm))
#                     elseif isdir(file_path)
#                         dfs_in_subdir = kge_read(file_path, ext)
#                     end
#                 end
#                 return(v)
#             end
    
    
#             #https://stackoverflow.com/questions/42499528/julia-convention-for-optional-arguments
#             #. Mostly you'll define two methods for your function:
#             function rp(file::AbstractString, lyr::Int)
#                 #xlyr = length(lyr)!=1 ? 1 : lyr
#                 ts=read(Raster(file,missingval=0))
#                 x = ts[t=lyr]
#                 contourf(x; c=cgrad(:thermal),size=(1200, 800))
#             end
    
#             function rp(file::AbstractString)
#                 ts=read(Raster(file,missingval=0))
#                 contourf(ts; c=cgrad(:thermal),size=(1200, 800))
#             end
    
    
#             function rp(file::AbstractString, lyr::Int)
#                 #xlyr = length(lyr)!=1 ? 1 : lyr
#                 ts=read(Raster(file,missingval=0))
#                 x = ts[t=lyr]
#                 contourf(x; c=cgrad(:thermal),size=(1200, 800))
#             end
    
#             function rplot(file::AbstractString)
#             	xr = read(Raster(file;crs=EPSG(25832),missingval=0))
#             	Plots.plot(xr;c=cgrad(:thermal),
#                 xlabel="",
#                 ylabel="",
#                 size=(1200*.66, 800*.66))   
#             end
    
    
#             function rplot(file::AbstractString, lyr::Int)
#             	xr = read(Raster(file;crs=EPSG(25832),missingval=0))
#             	Plots.plot(xr[t=lyr];c=cgrad(:thermal),size=(1200*.8, 800*.8))   
#             end
#             	#xr[t=20,cname="RdBl"]|>plot
#             	# c=:thermal]|>Plots.plot
    
#             function rplot(file::Raster, lyr::Int)
#                 xr = file
#                 Plots.plot(xr[t=lyr];c=cgrad(:thermal),size=(1200*.8, 800*.8))
#             end
    
#             function rpall(file::AbstractString)
#             	xr = read(Raster(file;crs=EPSG(25832),missingval=0))
#             	#xr = read(Raster(file;crs=EPSG(25832),missingval=-9999))
#             	Plots.plot(xr;c=cgrad(:thermal),
#                 xlabel="",
#                 ylabel="",
#                 size=(1200*.8, 800*.8))
#             end
    
#             function cpl(file::AbstractString,lyr::Int,msk::Int)
#                 x=read(Raster(file,missingval=0)) #-9999
#                 x=x[t=lyr]
#                 zm = x .< msk
#                 Plots.contourf(Rasters.mask(x; with=zm); 
#                 c=cgrad(:thermal),
#                 xlabel="",
#                 ylabel="",
#                 size=(1200, 800))
#             end
    
    
#             # function crpl(file::AbstractString,lyr::Int,vmin::Int,vmax::Int)
#             #     x=read(Raster(file,missingval=-9999)) #
#             #     x=x[t=lyr]
#             #     tmp=(x .> vmin);
#             #     mm=Rasters.mask(x; with=tmp);
#             #     zm =(x .< vmax);
#             #     Plots.contourf(Rasters.mask(mm; with=zm); c=cgrad(:thermal),size=(1200, 800),
#             #     xlabel="",ylabel="")
#             # end
    
#             # function crpl(file::AbstractString,lyr::Int,vmin::Int,vmax::Int)
#             #     x=read(Raster(file,missingval=-9999)) #
#             #     x=x[t=lyr]
#             #     x=Rasters.rebuild(x;missingval=vmin)
#             #     #tmp=(x .> vmin);
#             #     #mm=Rasters.mask(x; with=tmp);
#             #     zm =(x .< vmax);
#             #     Plots.contourf(Rasters.mask(x; with=zm); c=cgrad(:thermal),size=(1200, 800),
#             #     xlabel="",ylabel="")
#             # end
    
    
#             function cntplt(file::AbstractString)
#                 x=read(Raster(file,missingval=0))
#                 #describe(x)
#                 Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
#             end
    
#             # getf("sum")[2] |> cntplt
#             # v = getf("sum")[3]
#             # cntplt(v)
#             # rv = Raster(v)
#             # cntplt(rv)
#             # #describe(ar[2])  
    
#             function readras(file::AbstractString)
#                 x=read(Raster(file,missingval=0)) #read all in RAM
#                 describe(x)
#                 return(x)
#             end
    
#             function readras(path::Regex)
#                 "reads first match"
#                 v::Vector{String} = readdir();
#                 v = v[broadcast(x->endswith(x,"nc"),v)];
#                 file = v[(broadcast(x->occursin(path,x),v))] |>first;
#                 x::Raster = read(Raster(file,missingval=0))
#                 return(x)
#             end
    
#             function readrasrec(prefix::Regex)
#         """
#         readras(prefix::Regex)
#         reads first match of regex raster
#         """
#         rootdir="."
#         results = []
#         for (looproot, dirs, filenames) in walkdir(rootdir)
#             for filename in filenames
#                 #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                 if (occursin(prefix,filename)) && (endswith(filename,"nc"))
#                     push!(results, joinpath(looproot, filename)) 
#                 end
#             end
#         end
#         println(results)
#         file = first(results)
#         x = read(Raster(file,missingval=0)) #read all in RAM
#         describe(x)
#         return(x)
#             end
    
    
#             function cntplt(file::Raster{Union{Missing, Float64}, 2, Tuple{Y{Mapped{Float64, Vector{Float64}, DimensionalData.Dimensions.LookupArrays.ReverseOrdered, DimensionalData.Dimensions.LookupArrays.Regular{Float64}, DimensionalData.Dimensions.LookupArrays.Points, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}, EPSG, EPSG, Y{Colon}}}, X{Mapped{Float64, Vector{Float64}, DimensionalData.Dimensions.LookupArrays.ForwardOrdered, DimensionalData.Dimensions.LookupArrays.Regular{Float64}, DimensionalData.Dimensions.LookupArrays.Points, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}, EPSG, EPSG, X{Colon}}}}, Tuple{Dim{:t, DimensionalData.Dimensions.LookupArrays.Sampled{Float64, Vector{Float64}, DimensionalData.Dimensions.LookupArrays.ForwardOrdered, DimensionalData.Dimensions.LookupArrays.Regular{Float64}, DimensionalData.Dimensions.LookupArrays.Points, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}}}}, Matrix{Union{Missing, Float64}}, Symbol, DimensionalData.Dimensions.LookupArrays.Metadata{Rasters.NCDfile, Dict{String, Any}}, Missing})
#                 x=file
#                 describe(x)
#                 Plots.contourf(x; c=cgrad(:thermal))
#                 #,size=(1200, 800))
#             end
    
#             function cplt(file::AbstractString)
#                 x=read(Raster(file))
#                 Plots.contourf(x; c=cgrad(:thermal),size=(1200, 800))
#             end
    
#             function cpal(ext::AbstractString)
#                 path = pwd()
#                 files = readdir(path)
#             	for file in files
#                        i = joinpath(path, file)
#                        if isfile(i) && occursin(Regex(ext),file) && (!occursin("stack",file)) && endswith(file, ".nc")
#             	   outname=replace(i,"nc"=>"jl.png");
#                        #println(outname)
#                        r=read(Raster(i,missingval=0));
#                        p=Plots.contourf(r;
#             		title=replace(basename(i),".nc"=>""), #split(outname,"/")[end], #basename(i)
#             		c=cgrad(:thermal),
#             		size=(1200, 800));
#                        savefig(p,outname)
#                        println(basename(outname)," saved!");
#                     end
#                end
#             end
    
#             function stackplot(ext::AbstractString)
#                 path = pwd()
#                 files = readdir(path)
#             	for file in files
#                        i = joinpath(path, file)
#                        if isfile(i) && occursin(Regex(ext),file) && (occursin("stack",file)) && endswith(file, ".nc")
#             	   outname=replace(i,"nc"=>"jl.png");
#                        #println(outname)
#                        r=read(Raster(i,missingval=0,mappedcrs=EPSG(25832)));
#             	   #(i,missingval=-9999,mappedcrs=EPSG(25832))
#             	   ee = Int(r.dims[3][end])
#             	   rn = r[t=2:ee];    #subset till end
#                        p=Plots.plot(rn;
#             #		title=replace(basename(i),".nc"=>""), #no title cause problems
#             		c=cgrad(:thermal),
#             		size=(1200, 800));
#                        savefig(p,outname)
#                        println(basename(outname)," saved!");
#                     end
#                end
#             end
    
#             # flags = Dict(
#                 # :s_srs => "epsg:25832",
#                 # :t_srs => "epsg:4326",
#                 # :tr => [100,100],
#                 # :r => :near,
#             # )
#             # flags = Dict(
#                 # :tr => [100,100],
#                 # :r => :near,
#             # )
#             # warp(r[t=3],flags)  |> Plots.plot
    
#             function cc(ext::AbstractString)
#                 path = pwd()
#                 files = readdir(path)
#             	for file in files
#                        i = joinpath(path, file)
#                        if isfile(i) && occursin(Regex(ext),file) && endswith(file, ".nc")
#             	   outname=replace(i,"nc"=>"jl.png");
#                        #println(outname)
#                        r=read(Raster(i,missingval=-9999,mappedcrs=EPSG(25832)));
#                        #p=Plots.contourf(r;
#                        p=Plots.plot(r;
#             		title=replace(basename(i),".nc"=>""), #split(outname,"/")[end], #basename(i)
#             		c=cgrad(:thermal),
#             		size=(1200, 800));
#                        savefig(p,outname)
#                        println(basename(outname)," saved!");
#                     end
#                end
#             end


    
#             function ll()
#                 readdir()
#             end
    
#             # dfs=loadalldfs(p)
    
#             # T = typeof(dfs[2])
#             # for (name, typ) in zip(fieldnames(T), T.types)
#             #     println("type of the fieldname $name is $typ")
#             # end
    
#             # z=map(x->DataFrames.metadata(x)|>collect,dfs)
#             # z=map(x->basename(x[1][2]),z)
#             # findall("al",z[2])
#             # map(x->occursin(r"wol",x),z)
#             # map(x->findall(r"wol",x),z)
#             # #function grepl(df::DataFrame)
#             # using Grep
#             # grep("wlf",z)
#             # grep("wlf",z)|>getindex
    
#             #dfs[grep("wlf",z)]
#             #plotf()
    
    
#             function lplot(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match and plots..."
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                      map(x->basename(only(DataFrames.metadata(x))[2]),
#                      dfs))] |> first
#                      ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#                      nm = propertynames(df)[1:end-1];
#                      o = DataFrames.metadata(df)|>collect
#                      ti = basename(o[1][2])
#                      @df df Plots.plot(:date,cols(ln),yaxis=:log,title=ti)  
#             end
    
    
    
#             function lplot(df::DataFrame)
#                 nm = propertynames(df)[1:end-1];
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 @df df Plots.plot(:date,cols(ln),yaxis=:log,title=ti)     
#             end
    
#             function lplot(df::String)
#                 df=readdf(df)
#                 nm = propertynames(df)[1:end-1];
#                 o = collect(DataFrames.metadata(df))[1][2] |>basename
#                 ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 @df df Plots.plot(:date,cols(ln),yaxis=:log,title=o)     
#             end
    
#             function lplot(x::Regex)
#                 df=readdf(x)
#                 nm = propertynames(df)[1:end-1];
#                 o = collect(DataFrames.metadata(df))[1][2] |>basename
#                 ln = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 @df df Plots.plot(:date,cols(ln),yaxis=:log,title=o)     
#             end
    
#             #collect(DataFrames.metadata(df))[1][2]
#             #for i in dfs;collect(DataFrames.metadata(i))[1][2]|>basename|>println;end
    
    
#             function aplot(df::DataFrame)
#                 df[!,:year]=year.(df[!,:date]) ;
#                 s = Symbol.(filter(x->!occursin("date",x),names(df)))
#                 o = DataFrames.metadata(df)|>collect
#                 ti = "AndrewsPlot of "*basename(o[1][2])
#                 @df df andrewsplot(:year, cols(s), legend = :topleft,title=ti)
#             end
    
    
#             function getnames(dfs::Vector)
#             	nms = [];
#             	for i in dfs;	
#             		x=collect(DataFrames.metadata(i))[1][2]|>basename
#             		push!(nms, x)
#             	end
#             	return(nms)
#             end
#             function getnames(dfs::DataFrame)
#             	x=collect(DataFrames.metadata(dfs))[1][2]|>basename
#             	return(x)
#             end
    
    
#             function loadalldfs(path::Vector{Any})
#                 files = path
#                 dfs::Vector{DataFrame} = []
#                 for file in files
#                     if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
#                        file_path = file
#             	   println("reading ",file_path,"...")
#             	   p1 = readdf(file_path)
#             	   push!(dfs, p1)
#                     end
#                 end
#                 return(dfs)
#             end
    
#             function loadalldfs(path::Vector{String})
#                 files = path
#                 dfs::Vector{DataFrame} = []
#                 for file in files
#                     if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
#                        file_path = file
#             	   println("reading ",file_path,"...")
#             	   p1 = readdf(file_path)
#             	   push!(dfs, p1)
#                     end
#                 end
#                 return(dfs)
#             end
    
#             function loadalldfs(path::Regex)
#                 v = readdir();
#                 v = v[broadcast(x->!endswith(x,"nc"),v)];
#                 files = v[(broadcast(x->occursin(path,x),v))];
#                 dfs::Vector{DataFrame} = []
#                 for file in files
#                     if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
#                        file_path = file
#             	   println("reading ",file_path,"...")
#             	   p1 = readdf(file_path)
#             	   push!(dfs, p1)
#                     end
#                 end
#                 return(dfs)
#             end
    
#             readall = loadalldfs
    
    
#             function loadalldfs(path::AbstractString)
#                 files = readdir(path)
#                 dfs = DataFrame[]
#                 #nms = []
#                 for file in files #&& occursin(Regex(prefix),file)
#                     if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",file))
#                        file_path = joinpath(path, file)
#             	   println("reading ",file_path,"...")
#             	   p1 = readdf(file_path)
#             	   push!(dfs, p1)
#             	   #push!(nms, file)
#                     end
#                 end
#                 return(dfs)
#                 #return(nms)
#             end
    
#             function listdfs(path::AbstractString)
#         files = readdir(path)
#         #dfs = DataFrame[]
#         nms = []
#         for file in files #&& occursin(Regex(prefix),file)
#             if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
#                file_path = joinpath(path, file)
#     	   println("reading ",file,"...")
#     	   #p1 = readdf(file_path)
#     	   #push!(dfs, p1)
#     	   push!(nms, file)
#             end
#         end
#         #return(dfs)
#         return(nms)
#             end
    
#             function vars()
#                 varinfo()
#             end
    
#             function vars(pt::AbstractString)
#                 #varinfo(Core,r".*field.*")
#                 #varinfo(Main,r".*load*")
#                 varinfo(Main,Regex(".*pt*"))
#             end
    
    
    
#             #if (occursin(Regex(prefix,"i"),filename))
#             function regand(v::Vector{String},x1::AbstractString,y1::AbstractString)
#                 needle=join([x1,y1],"+.*");
#                 z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))] 
#             return(z)
#             end
    
#             # regand(v,"sum","nc\\E")
#             # x1,y1="sum","nc\\E"
#             # needle=join([x1,y1],"+.*");
#             # Regex(needle,"i")
#             #xv=("utm_rcm", "rcm-c4")    
#             function regand(v::Vector{String},xv::Tuple{String, String})
#                 needle=join([xv[1],xv[2]],"+.*");
#                 z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
#             return(z)
#             end
    
#             #function regand(v::Vector{String},xv::Tuple{Symbol,Symbol})
#             function regand(v::Vector{String},xv::Vector{Symbol})
#                 needle=join(xv,"+.*");
#                 z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
#             return(z)
#             end
    
#             function regand(v::Vector{String},xv::Vector{String})
#                 needle=join(xv,"+.*");
#                 z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
#             return(z)
#             end
    
    
#             function regand(v::Vector{String},xv::Regex)
#                 "here you can put any regex to filter the Vector"
#                 z = v[(broadcast(x->occursin(xv,x),v))] 
#             return(z)
#             end
    
    
    
#             function nconly(x1::AbstractString)
#             v::Vector{String} = readdir();
#             v = v[broadcast(x->endswith(x,"nc"),v)];
#             z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
#             return(z)
#             end
    
#             #tats. map is same here
#             #v[map(x->endswith(x,"nc"),v)]  
#             #https://stackoverflow.com/questions/52892726/julia-whats-the-difference-between-map-and-broadcast
    
#             function nconly(Any)
#             v = readdir();
#             z = v[broadcast(x->endswith(x,"nc"),v)];
#             return(z)
#             end
    
    
    
#             function readallras(path::AbstractString)
#                 v = readdir(path);
#                 v = v[broadcast(x->endswith(x,"nc"),v)];
#                 z=[];
#                 for s in v; 
#                 #if contains(x1,s) & occursin(r"nc$",s)
#                 ts=read(Raster(s,missingval=0))
#                 push!(z,ts);
#                 end
#                 return(z)
#             end
    
#             function readallras(path::AbstractString, ex::AbstractString)
#                 v = readdir(path);
#                 v = v[broadcast(x->endswith(x,"nc") & occursin(ex,x),v)];
#                 z=[];
#                 for s in v; 
#                 #if contains(x1,s) & occursin(r"nc$",s)
#                 ts=read(Raster(s,missingval=0))
#                 push!(z,ts);
#                 end
#                 return(z)
#             end
    
#             function readallras(ex::Regex)
#                 v = readdir(".");
#                 v = v[broadcast(x->endswith(x,"nc") & occursin(ex,x),v)];
#                 z=[];
#                 for s in v; 
#                 #if contains(x1,s) & occursin(r"nc$",s)
#                 ts=read(Raster(s,missingval=0))
#                 push!(z,ts);
#                 end
#                 return(z)
#             end
    
    
    
#             function hometeo()
#                 cd("/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/")
#                 println("you are here: ",pwd())
#             end
    
#             function homreg()
#                 cd("/mnt/d/Wasim/regio/out/");
#                 println("you are here: ",pwd())
#                 fd()
#             end
    
    
#             function writewa(file::AbstractString, df::DataFrame)
#                 dout = df
#                 dout.YY = map(x ->year(x),dout.date)
#                 dout.MM = map(x ->month(x),dout.date)
#                 dout.DD = map(x ->day(x),dout.date)
#                 dout[!, "HH"] .= 24
#                 #df = select!(df,Symbol.(filter(x->!occursin("date",x), names(df))))
#                 #dout = select(df, Not(:date))
#                 #dout = dout[!,Cols([:YY,:MM,:HH,:DD],1:end-4)]
#                 dout = dout[!,Cols([:YY,:MM,:HH,:DD],Not(Cols(r"date")))]
#                 #cls = propertynames(df)|>sort|>reverse
#                 #df = df[!,cls[2:end]] 
#                 CSV.write(file, dout, transform = (col, val) -> something(val, missing),delim="\t")  
#                 nothing
#             end
    
#             function writedf(file, table)
#                 CSV.write(file, table, transform = (col, val) -> something(val, missing),delim="\t")  
#                 nothing
#             end
    
#             function writedesc(file, table)
#                 CSV.write(file, describe(table), transform = (col, val) -> something(val, missing),delim="\t")  
#                 nothing
#             end
    
#             #wc -l in julia:
#             function wcl(file::AbstractString)
#                 open(file) do f
#                     println(count(_ -> true, eachline(f)))
#                 end
#             end
    
#             function wcl(file::AbstractString,Bool)
#                 open(file) do f
#                     ct=(count(_ -> true, eachline(f)))
#                     #println(file,ct)
#                     println("$file:\t $ct")
#                 end
#             end
    
#             #wcl(file,true)
    
    
#             function vgg(regex::AbstractString, ending::AbstractString)
#                 cmd = `grep --color=always -C2 -rIHn -E "$regex" --include="*.$ending"`
#                 run(cmd)
#             end
    
#             ##julia with no regex
    
#             function vg(snippet::AbstractString, file_ending::AbstractString)
#                 files = filter(file -> endswith(file, file_ending), readdir())
#                 # loop over each file
#                 for file in files
#                     open(file) do f
#                         counter = 0 # Zähler initialisieren
#                         for line in eachline(f)
#                             counter += 1 # Zähler erhöhen
#                             # check if the line matches the regex
#                             #if occursin(Regex(regex), line)
#                             if contains(line,snippet)
#             #                    println("$file: $counter:\t $line")
#                                 printstyled("$counter:\t",color=:light_red) 
#                                 printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
#                                 printstyled("$line\n",color=:green,bold=true) 
#                             end
#                         end
#                     end
#                 end
#             end
    
#             #faster and pure julia:
#             function vgrep(regex, file_ending)
#         # list files that start with "qgko" and end with file_ending
#         #files = filter(file -> startswith(file, "qgko") && endswith(file, file_ending), readdir())
#         files = filter(file -> endswith(file, file_ending), readdir())
#         # loop over each file
#         for file in files
#             # open the file and read each line
#             #xlines = readlines(file)
#             #filter(z -> true, xlines) |> (x -> for i in 1:length(x) println("$i\t$(x[i])") end)
    
#             open(file) do f
#                 counter = 0 # Zähler initialisieren
#                 for line in eachline(f)
#                     counter += 1 # Zähler erhöhen
#                     # check if the line matches the regex
#                     if occursin(Regex(regex,"i"), line)
#                         printstyled("$counter:\t",color=:light_red) 
#                         printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
#                         printstyled("$line\n",color=:green,bold=true) 
#                         # print the file name, line number and line content
#                         #println("$file:$(f.lineno):$line") <-nope
#                         #m=match(regex, line)
#                         #m=count(_ -> true, line) #das zählt die linechars
#                         #println("$file: $counter:\t $line")
#                     end
#                 end
#             end
#         end
#             end
    
#             #cd("/mnt/d/Wasim/regio/out/");
#             #cd("/mnt/d/Wasim/streu/out");
    
#             #using GeoArrays
#             function gplot(r::AbstractString)
#                 geoarray = GeoArrays.GeoArray(ArchGDAL.readraster(lk))
#                 geoarray|>plot
#             end
    
#             function rglob(prefix::AbstractString)
#                 rootdir="."
#                 results = []
#                 for (looproot, dirs, filenames) in walkdir(rootdir)
#                     for filename in filenames
#                         #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                         if (occursin(Regex(prefix,"i"),filename))
#                             push!(results, joinpath(looproot, filename)) 
#                         end
#                     end
#                 end
#                 return results
#             end
    
#             function fdf(df::DataFrame)
#                 nrows=size(df)[2]-1 
#                 o = DataFrames.metadata(df)|>collect
#                 ti = basename(o[1][2])
#                 fig = make_subplots(
#                     shared_xaxes=true, 
#                     shared_yaxes=true    
#                     );
#                 for i in 1:nrows;
#                     add_trace!(fig, 
#                     PlotlyJS.scatter(x=df.date, y=df[:,i],
#                     name=names(df)[i]));
#                 end
#                 fact = 0.7
#                 logy = true;
#                 if logy == true
#                     PlotlyJS.relayout!(fig,
#                     template="seaborn",
#                     yaxis_type="log",
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 else
#                     PlotlyJS.relayout!(fig,
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 end
#                 println(describe(df))
#                 println("showing plot...")
#                 display(fig)
#             end
    
    
#             function xdf(df::DataFrame)
#                 try
#                     nrows=size(df)[2]-1 
#                     fig = make_subplots(
#                         shared_xaxes=true, 
#                         shared_yaxes=true    
#                         );
#                     for i in 1:nrows;
#                         add_trace!(fig, 
#                         PlotlyJS.scatter(x=df.date, y=df[:,i],
#                         name=names(df)[i]));
#                     end
#                     PlotlyJS.relayout!(fig,
#                     template="plotly_dark",
#                     yaxis_type="log")
#                     display(fig)
#                 catch e
#                     println("An error occurred: ", e)
#                 finally
#                     println("showing plot...")
#                     println(describe(df))
#             end
#             end
    
#             #go dir up
#             function cdb()
#                 dirname(pwd())|>cd
#                 pwd()|>println
#             end
    

#             function vgr(regex, file_ending)
#                 rootdir=pwd()
#                 println("starting on: $rootdir...\n searching for >> $regex << with file ending >> $file_ending <<\n")
#                 files = []
#                 for (looproot, dirs, filenames) in walkdir(rootdir)
#                     for filename in filenames
#                         #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                         #if (occursin(Regex(prefix,"i"),filename))
#                         if (endswith(filename, file_ending))
#                             push!(files, joinpath(looproot, filename)) 
#                         end
#                     end
#                 end
#                 #files = filter(file -> endswith(file, file_ending), readdir())
#                 for file in files
#                     open(file) do f
#                         counter = 0 # Zähler initialisieren
#                         for line in eachline(f)
#                             counter += 1 # Zähler erhöhen
#                             if occursin(Regex(regex,"i"), line)
#                                 println("$file: $counter:\t $line")
#                             end
#                         end
#                     end
#                 end
#             end
    
#             function vjl(regex)
#                 # greps jl from julia folder
#                 pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia";
#                 file_ending=".jl"
#                 files = filter(file -> endswith(file, file_ending), readdir(pt,join=true))
#                 for file in files
#                     open(file) do f
#                         counter = 0
#                         for line in eachline(f)
#                             counter += 1
#                             if occursin(Regex(regex,"i"), line)
#                                 println("$file: $counter:\t $line")
#                             end
#                         end
#                     end
#                 end
#             end
    
#             function median_filter(ras::Raster)
#                 # Get the array and dimensions of the raster
#                 Z=Band(1)
#                 #arr = ras[:Z]
#                 arr = ras[Z]
#                 nx, ny = size(arr)
#                 # Create an output array with the same size and type
#                 out = similar(arr)
#                 # Loop over the pixels, excluding the borders
#                 for i in 2:nx-1, j in 2:ny-1
#                   # Get the values in the 3x3 window
#                   window = arr[i-1:i+1, j-1:j+1]
#                   # Calculate the median of the window
#                   out[i,j] = median(window)
#                 end
#                 # Return a new raster with the filtered array
#                 return rebuild(ras,out)
#             end
    
#             function dfilter(regex::AbstractString,dfs::Vector{DataFrame})
#                 filter(n->occursin(Regex(regex,"i"),n),
#                 map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
#                 )
#             end
    
#             #dfilter("cl",dfs)
#             #typeof(dfs)
#             #dfs
#             #regex="qout"
#             #z::Vector{Raster} = 
#             function filterplot(regex::AbstractString,z::Vector{Raster})
#                 "selects first match and plots..."
#                 rr::Raster = z[map(n->occursin(Regex(regex,"i"),n),
#                                #map(x->map(String,name(x)),z)
#                                map(String,map(name,z))
#                                )                   
#                                ] |> first
#                 #Plots.contourf(rr; c=cgrad(:thermal),size=(1200*.7, 800*.7))
#                 ti::Symbol = name(rr)
#                 fct::AbstractFloat = 0.5
#                 Plots.plot(rr; 
#                     c=cgrad(:thermal),
#                     title=ti,
#                     size=(1200*fct, 800*fct)
#                     )
#             end
    
#             #filterplot("tem",rds)
#             #filterplot("vap",rds)
    
#             function filterplot(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match and plots..."
#                 df::DataFrame = dfs[map(n->occursin(Regex(regex,"i"),n),
#                 map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
#                 )] |> first
#                 dfp(df)
#                 #indexin(1:length(dfs),
#                 #map(x->basename(only(DataFrames.metadata(x))[2]),dfs))
#             end
    
#             function filterplot!(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match and add to plot..."
#                 df::DataFrame = dfs[map(n->occursin(Regex(regex,"i"),n),
#                 map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
#                 )] |> first
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 #Symbol(names(df))
#                 #s = propertynames(df)[Not(end)] #geht auch, aber positionsabhängig
#                 @df df Plots.plot!(:date,cols(s),legend = :topright)
#             end
    
#             # filterplot("win",dfs)
#             # filterplot!("qg",dfs)
#             # filterplot!("qout",dfs)
#             #dfs[5] |>dfp
#             #convert(DataFrame,dfs[5]|>DataFrames.metadata|>only)
    
#             function fread(ext::AbstractString)
#                 cwd = pwd() 
#                 m = DataFrame[]
#                 for (root, dirs, files) in walkdir(cwd)
#                  for file in files
#                  if isfile(file) && occursin(Regex(ext),file)&&
#                     (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",file))
#             	 nm=joinpath(root, file)
#             	 push!(m,readdf(nm))
#                  end
#                 end 
#                 end 
#                  return(m)
#             end 
    
#             # function dfpall(dfs::Vector{DataFrame})
#             #     for xx in 2:length(dfs)
#             #         dfp(first(dfs))
#             #         y = filter(x->!occursin("date",x),
#             #         names(dfs[xx]))
#             #         s = map(y -> Symbol(y),y)
#             #         px=@df dfs[xx] Plots.plot!(:date,cols(s),legend = :topright)
#             #         return(px)
#             #     end
#             # end
    
#             function dfpall(files::Vector{Any})
#                 "reads, reduces + merges by date and plots"
#                 #files
#                 dfs = DataFrame[]
#                 for file in files
#                     if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
#                        file_path = file
#             	   println("reading ",file_path,"...")
#             	   p1 = readdf(file_path)
#             	   push!(dfs, p1)
#                     end
#                 end
#                 df = reduce((left, right) -> 
#                   innerjoin(left, right, on = :date,makeunique=true), 
#                   dfs)
#                 y = filter(x->!occursin("date",x), names(df))
#                 s = map(y -> Symbol(y),y)
#                 @df df Plots.plot(:date,
#                         cols(s),
#                         #yaxis = :log,
#                         #legend = :bottom)
#                         legend = false)
#             end
    
    
#             function dfpall(dfs::Vector{DataFrame})
#         "reduces + merges by date and plots"
#         df = reduce((left, right) -> 
#           innerjoin(left, right, on = :date,makeunique=true), 
#           dfs)
#         y = filter(x->!occursin("date",x), names(df))
#         s = map(y -> Symbol(y),y)
#         @df df Plots.plot(:date,
#                 cols(s),
#                 #yaxis = :log,
#                 legend = :bottom)
#             end
    
    
#             function mall(files::Vector{Any})
#                 "reads, reduces + merges by date"
#                 #files
#                 dfs = DataFrame[]
#                 for file in files
#                     if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
#                        file_path = file
#             	   println("reading ",file_path,"...")
#             	   p1 = readdf(file_path)
#             	   push!(dfs, p1)
#                     end
#                 end
#                 df = reduce((left, right) -> 
#                   innerjoin(left, right, on = :date,makeunique=true), 
#                   dfs)
#                 return(df)
#             end
    
#             function mall(files::Vector{DataFrame})
#                 "reduces + merges by date"
#                 df = reduce((left, right) -> 
#                   innerjoin(left, right, on = :date,makeunique=true), 
#                   files)
#                 return(df)
#             end
    
#             function getdf(regex::AbstractString,dfs::Vector{DataFrame})
#                 "selects first match..."
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                      map(x->basename(only(DataFrames.metadata(x))[2]),
#                      dfs))] |> first
#                 return(df)
#             end
    
#             function getdf(dfs::Vector{DataFrame},index::Integer)
#         "selects by index.."
#         df = getindex(dfs,index)
#         return(df)
#             end
    
#             function getdf(regex::AbstractString,dfs::Vector{Any})
#                 "selects by regex.."
#                 #df = getindex(dfs,index)
    
#                 #y = filter(x->!occursin("date",x),names(df))
#                 al = filter(x->occursin(Regex(regex,"i"),x),dfs)
    
#                 df = filter(x->!occursin(r"txt|yrly|nc|png|svg|ftz_0|ftz",x),al)|>
#                 first |>readdf
    
#                 return(df)
#             end
    
    
#             function filterdf(dfs::Vector{Any})
#                 "selects presumably dfs from vector..."
#                 df = filter(x->!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",x),dfs)
#                 df = loadalldfs(df) #geht
#                 return(df)
#             end
    
#             function filterdf(regex::AbstractString,dfs::Vector{DataFrame})
#                 """
#                 for namestring, see as dfilter
#                 selects df from dfvector...
#                 same as getdf
#                 example:
#                 filterdf("clou",dfs)|>bardfm 
#                 """
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                      map(x->basename(only(DataFrames.metadata(x))[2]),
#                      dfs))] |> first
                
#                 # filter(n->occursin(Regex(regex,"i"),n),
#                 # map(x->basename(only(DataFrames.metadata(x))[2]),
#                 # dfs)
#                 # )
#             end
    
#             cs= [:default		,
#                   :blues		,
#                   :bluesreds		,
#                   :darkrainbow		,
#                   :darktest		,
#                   :grays		,
#                   :greens		,
#                   :heat		,
#                   :lightrainbow		,
#                   :lighttest];
    
    
#             # function dfpall(dfs::Vector{DataFrame})
#             #     odf=DataFrame[]
#             #     odf = map(x->innerjoin(dfs[1],x,on="date",makeunique=true),dfs)
    
#             #     # for i in dfs; 
#             #     #     push!(odf,innerjoin(dfs[1],i,on="date",makeunique=true))
#             #     # end 
#             #     out = reduce(hcat, odf) 
#             # end
    
    
#             # z=fread("eva")
#             # dfpall(z)
#             #filterplot("ev",z)
    
#             # dfs = z
#             # for df in 2:length(dfs)
#             #     y = filter(x->!occursin("date",x),
#             #     names(dfs[df]))
#             #     s = map(y -> Symbol(y),y)
#             #     println(s,df,dfs[df])
#             # end
    
    
#             #xdf(df)
#             #fdf(df)
    
    
#             function p()
#                 return(pwd())
#             end
    
#             function yrsum(x::String)
#                 df = readdf(x)
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 #ti=DataFrames.metadata(df)|>only|>last|>basename
#                 df[!, :year] = year.(df[!,:date]);
#                 df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
#                 return(df_yearsum)
#             end
    
#             function yrsum(x::DataFrame)
#                     df = x
#                     y = filter(x->!occursin("date",x),names(df))
#                     s = map(y -> Symbol(y),y)
#                     #ti=DataFrames.metadata(df)|>only|>last|>basename
#                     df[!, :year] = year.(df[!,:date]);
#                     df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
#                     return(df_yearsum)
#             end
    
#             function yrmean(x::String)
#                     df = readdf(x)
#                     y = filter(x->!occursin("date",x),names(df))
#                     s = map(y -> Symbol(y),y)
#                     #ti=DataFrames.metadata(df)|>only|>last|>basename
#                     df[!, :year] = year.(df[!,:date]);
#                     df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
#                     return(df_yearsum)
#             end
    
#             function yrmean(x::DataFrame)
#                 df = x
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 #ti=DataFrames.metadata(df)|>only|>last|>basename
#                 df[!, :year] = year.(df[!,:date]);
#                 df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
#                 return(df_yearsum)
#             end
    
#             function bardf(x::String)
#                 "with String"
#                 df = readdf(x)
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
#                 ti=DataFrames.metadata(df)|>only|>last|>basename
#                 df[!, :year] = year.(df[!,:date]);
#                 df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
#                 @df df_yearsum Plots.plot(:year,
#                     cols(s),
#                     legend = :topright, 
#                     title=ti,
#                     seriestype=:bar)
#             end
    
#             function bardf(x::Regex)
#                 "with regex, and new metadata extraction"
#                 df = readdf(x)
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
#                 ti=DataFrames.metadata(df)|>only|>last|>basename
#                 df[!, :year] = year.(df[!,:date]);
#                 df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
#                 @df df_yearsum Plots.plot(:year,
#                     cols(s),
#                     legend = :topright, 
#                     title=ti,
#                     seriestype=:bar)
#             end
    
#             function bardf(x::DataFrame)
#         "with DataFrame input"
#             df = x
#             y = filter(x->!occursin("date",x),names(df))
#             s = map(y -> Symbol(y),y)
#             ti=DataFrames.metadata(df)|>only|>last|>basename
#             df[!, :year] = year.(df[!,:date]);
#             df_yearsum = combine(groupby(df, :year), y .=> sum .=> y);
#             @df df_yearsum Plots.plot(:year,
#                 cols(s),
#                 legend = :topright, 
#                 title=ti,
#                 seriestype=:bar)
#             end
    
#             function bardfm(x::String)
#         "with String"
#         df = readdf(x)
#         y = filter(x->!occursin("date",x),names(df))
#         s = map(y -> Symbol(y),y)
#         ti=DataFrames.metadata(df)|>only|>last|>basename
#         df[!, :year] = year.(df[!,:date]);
#         df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
#         @df df_yearsum Plots.plot(:year,
#             cols(s),
#             legend = :topright, 
#             title=ti,
#             seriestype=:bar)
#             end
    
#             function bardfm(x::Regex)
#                     "with regex, and new metadata extraction"
#                     df = readdf(x)
#                     y = filter(x->!occursin("date",x),names(df))
#                     s = map(y -> Symbol(y),y)
#                     #ti=DataFrames.metadata(df)|>collect|>only|>last|>basename
#                     ti=DataFrames.metadata(df)|>only|>last|>basename
#                     df[!, :year] = year.(df[!,:date]);
#                     df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
#                     @df df_yearsum Plots.plot(:year,
#             cols(s),
#             legend = :topright, 
#             title=ti,
#             seriestype=:bar)
#             end
    
#             function bardfm(x::DataFrame)
#                 "with DataFrame input"
#                     df = x
#                     y = filter(x->!occursin("date",x),names(df))
#                     s = map(y -> Symbol(y),y)
#                     ti=DataFrames.metadata(df)|>only|>last|>basename
#                     df[!, :year] = year.(df[!,:date]);
#                     df_yearsum = combine(groupby(df, :year), y .=> mean .=> y);
#                     @df df_yearsum Plots.plot(:year,
#                         cols(s),
#                         legend = :topright, 
#                         title=ti,
#                         seriestype=:bar)
#             end
    
#             function cnt()
#                 return(length(readdir(pwd())))
#             end

    
#             function fdf(regex::AbstractString,dfs::Vector{DataFrame},f::Function)
#                 "selects first match and applies function..."
#                 df = dfs[map(n->occursin(Regex(regex,"i"),n),
#                 map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
#                 )] |> first
#                 f(df)
            
#                 #like: fdf("win",dfs,yrsum) 
#                 #like: fdf("win",dfs,describe) 
#                 #indexin(1:length(dfs),
#                 #map(x->basename(only(DataFrames.metadata(x))[2]),dfs))
#             end
    
    
#             function route_from_dir(dir::String)
#                 dirs = readdir(dir)
#                 routes::Vector{String} = []
#                 for directory in dirs
#                     if isfile("$dir/" * directory)
#                         push!(routes, "$dir/$directory")
#                     else
#                         if ~(directory in routes)
#                             newread = dir * "/$directory"
#                             newrs = route_from_dir(newread)
#                             [push!(routes, r) for r in newrs]
#                         end
#                     end
#                 end
#                 routes
#             end
    
#             ##https://chifi.dev/weve-been-writing-julia-wrong-speed-up-julia-with-annotations-9144845c3e24
#             # route_from_dir(".")
#             # ll()
    
#             mutable struct LinearRegression{A<:AbstractFloat, B<:AbstractFloat, P<:Function}
#                 a::A
#                 b::B
#                 predict::P
#                 function LinearRegression(x::Array,y::Array)
#                     # a = ((∑y)(∑x^2)-(∑x)(∑xy)) / (n(∑x^2) - (∑x)^2)
#                     # b = (x(∑xy) - (∑x)(∑y)) / n(∑x^2) - (∑x)^2
#                     if length(x) != length(y)
#                         throw(ArgumentError("The array shape does not match!"))
#                     end
#                     Σx::Float64 = sum(x)
#                     Σy::Float64 = sum(y)
#                     xy::Array = x .* y
#                     Σxy::Float64 = sum(xy)
#                     x2::Array{Float64} = x .^ 2
#                     Σx2::Float64 = sum(x2)
#                     n::Int64 = length(x)
#                     # Calculate a
#                     a::Float64 = (((Σy) * (Σx2)) - ((Σx * (Σxy)))) / ((n * (Σx2))-(Σx^2))
#                     b::Float64 = ((n*(Σxy)) - (Σx * Σy)) / ((n * (Σx2)) - (Σx ^ 2))
#                     predict(xt::Array) = (xt = [i = a + (b * i) for i in xt]::Array)
#                     return new{Float64, Float64, Function}(a::Float64, b::Float64, predict::Function)
#                 end
#             end
    
#             # x = randn(50000000)
#             # y = randn(50000000)
#             # @time LinearRegression(x, y).predict(y)
    
    
#             function dfyrs(df::DataFrame;)
#                 ti = DataFrames.metadata(df)|>only|>last|>basename
#                 fact,logy = 1,0
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 df[!, :year] = year.(df[!,:date]);
#                 df = combine(groupby(df, :year), y .=> sum .=> y);
#                 #df = df[!,Not("date")]
#                 ln = Symbol.(filter(x->!occursin("year",x),names(df)))
#                 nrows=size(df)[2]-1
#                 if nrows == 1
#                     ln = only(ln)
#                     fig = 
#                     PlotlyJS.plot(
#                     PlotlyJS.scatter(x=df.year, y=df[!,ln],
#                     name=ln,type="bar")
#                     );
#                     PlotlyJS.relayout!(fig,
#                         height=600*fact,width=900*fact,
#                         title_text="Series of "*ti)
#                 else
#                     fig = PlotlyJS.make_subplots(
#                         shared_xaxes=true, 
#                         shared_yaxes=true    
#                         );
#                     for i in ln
#                         PlotlyJS.add_trace!(fig, 
#                         PlotlyJS.scatter(x=df.year, y=df[:,i],
#                         name=i));
#                     end
#                     if logy == true
#                         PlotlyJS.relayout!(fig,yaxis_type="log",
#                         height=600*fact,width=900*fact,
#                         title_text="Series of "*ti)
#                     else
#                         PlotlyJS.relayout!(fig,
#                         height=600*fact,width=900*fact,
#                         title_text="Series of "*ti)
#                     end
#                 end
#                 display(fig)
#             end
    
#             function dfpjs(df::DataFrame;)
#                 nrows=size(df)[2]-1 
#                 ti = DataFrames.metadata(df)|>only|>last|>basename
#                 fig = PlotlyJS.make_subplots(
#                     shared_xaxes=true, 
#                     shared_yaxes=true    
#                     );
#                 for i in 1:nrows;
#                     PlotlyJS.add_trace!(fig, 
#                     PlotlyJS.scatter(x=df.date, y=df[:,i],
#                     name=names(df)[i]));
#                 end
#                 fact,logy = 1,0
#                 if logy == true
#                     PlotlyJS.relayout!(fig,yaxis_type="log",
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 else
#                     PlotlyJS.relayout!(fig,
#                     height=600*fact,width=900*fact,
#                     title_text="Series of "*ti)
#                 end
#                 display(fig)
#             end
    
    
    
#             function monsum(x::String)
#                 df = readdf(x)
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 df[!, :month] = month.(df[!,:date]);
#                 df_monthsum = combine(groupby(df, :month), y .=> sum .=> y);
#                 return(df_monthsum)
#             end
    
#             function monsum(x::DataFrame)
#                 df = x
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 df[!, :month] = month.(df[!,:date]);
#                 df_monthsum = combine(groupby(df, :month), y .=> sum .=> y);
#                 return(df_monthsum)
#             end
    
#             function monmean(x::String)
#                 df = readdf(x)
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 df[!, :month] = month.(df[!,:date]);
#                 df_monthsum = combine(groupby(df, :month), y .=> mean .=> y);
#                 return(df_monthsum)
#             end
    
#             function monmean(x::DataFrame)
#                 df = x
#                 y = filter(x->!occursin("date",x),names(df))
#                 s = map(y -> Symbol(y),y)
#                 df[!, :month] = month.(df[!,:date]);
#                 df_monthsum = combine(groupby(df, :month), y .=> mean .=> y);
#                 return(df_monthsum)
#             end
    
#             function barp(x::DataFrame)
#                 "with DataFrame input"
#                     df = x
#                     ti=DataFrames.metadata(df)|>only|>last|>basename
#                     if any(x->occursin("year",x),names(df))
#                         ln = Symbol.(filter(x->!occursin("year",x),names(df)))
#                         @df df Plots.plot(:year,
#                             cols(ln),
#                             legend = :topright, 
#                             title=ti,
#                             seriestype=:bar) #color=:lightrainbow
#                     elseif any(x->occursin("month",x),names(df))
#                         ln = Symbol.(filter(x->!occursin("month",x),names(df)))
#                         @df df Plots.plot(:month,
#                             cols(ln),
#                             legend = :topright, 
#                             title=ti,
#                             seriestype=:bar)
#                     elseif (
#                         any(x->occursin("month",x),names(df)) & 
#                         any(x->occursin("year",x),names(df))            
#                         )
#                         ln = (filter(x->!occursin("month",x),names(df)))
#                         ln = Symbol.(filter(x->!occursin("year",x),ln))
#                         @df df Plots.plot(:month,
#                             cols(ln),
#                             legend = :topright, 
#                             title=ti,
#                             seriestype=:bar)
#                     else
#                         dfp(df)        
#                     end
#             end
    
#             function so_read(x::AbstractString)
#                 "--- reader with drop exept of first col ---"
#                 ms=["-9999","lin","log"]
#                 df::DataFrame = CSV.read(x,DataFrame,
#                 missingstring=ms,
#                 types = Float64,
#                 delim="\t",
#                 silencewarnings=true,
#                 normalizenames=true,
#                 drop=(i, nm) -> i == 4) |> dropmissing
#                 df.YY=map(x ->Int(x),df.YY);
#                 df.MM=map(x ->Int(x),df.MM);
#                 df.DD=map(x ->Int(x),df.DD);
#                 df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#                 df=df[:,Cols(4,end)] #only 1stlevel and date col.
#                 DataFrames.metadata!(df, "filename", x, style=:note);
#             end
#         end
#     end
    
#     println("you are here: ",p())

# end
