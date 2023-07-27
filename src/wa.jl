module wa
using DataFrames, CSV, Statistics, Dates, StatsPlots, Distributions
using DelimitedFiles, Grep , Printf
using Rasters
using PrettyTables
import NCDatasets
default(show = true)

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

#"dySeries"|>vgr

function vgpyo(snippet::AbstractString)
    owd=pwd()
    cd("C:/Users/Public/Documents/Python_Scripts")
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

# cdb()
# fdd()


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
function filterplot(regex::AbstractString,dfs::Vector{DataFrame})
    "selects first match and plots..."
    df = dfs[map(n->occursin(Regex(regex,"i"),n),
    map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
    )] |> first
    dfp(df)
    #indexin(1:length(dfs),
    #map(x->basename(only(DataFrames.metadata(x))[2]),dfs))
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

#dfpall(dxx)
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
        
#kge_read(pwd(),"out")
#println("try ",kge_read(pwd(),"out"))

# md = loadso(p,"so")
# md[end-1]|>vio

# z = loadalldfs(p)
# vv = listdfs(p)
# px = z[end-1]|>vio;
# plot!(title=vv[end-1])

#repl location
#C:\Users\chs72fw\.julia\logs\repl_history.jl

# #Rcalls
# using RCall
# jmtcars = reval("mtcars");
# gg = rimport("ggplot2")
# gg.ggplot(jmtcars,gg.aes(:mpg, :wt))+gg.geom_point()

# @rlibrary "data.table"
# # mode: r
# wa.dd <- function(x,flag=T) {
#     x <-
#       data.table::fread(
#         file = x,
        
#         header = T,
#         check.names = T,
#         skip = ifelse(test=flag, yes="YY",no=0),
#         na.strings = "-9999"
#       )
#     x = na.omit(x[, lapply(.SD, function(z)
#       as.numeric(as.character(z)))], cols = 1:4)
#     x[, t := do.call(paste, c(.SD, sep = "-")), .SDcols = 1:4]
#       x = x[, t := data.table::as.IDate(t)][, -c(1:4)]
#     data.table::setkey(x = x, 't')
#     }
# # time: 2023-03-30 11:32:45 W. Europe Summer Time
# # mode: r
#   oo=wa.dd("preci_1970.txt")
#   md = rcopy(R"oo")
#   describe(md)

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


#df=dfr(r"qgk")
#"D:/Wasim/regio/out/rc200/v5/station4/"|>cd
# df = readdf(r"qg")
# #df = dfr(r"qg")
# df = readdf(r"wind")
# df = dfr(r"ou")
# ftp(df)
# lplot(df)
# lplot(df[!,Cols(r"date|131|6")])
# ftp(df[!,Cols(r"date|135|136")])

# df[!,Cols(r"date|135|136")]|>qqpl
# #"D:\\Wasim\\streu\\out\\cw0"|>cd

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
# r"qi_"|>readras|>descr

# r"qi_"|>dfp
# using PlotlyJS
# plotlyjs()
# r"qi_"|>dfbarjs
# r"qi_"|>bardf

# v=r"qi_"|>globdf
# dfp(second(v))
#descr(r)
#lastdim
#d=(name(r.dims)[end])
#for i in range(1,(size(r)[end]));describe(r[d(i)]);end #nope
#describe(r[(r.dims)[end]=1])

# r = readras(r"tso")
# r[:,:,[2,3]]|>plot
# for i in range(1,(size(r)[end]));
#     sm=r[:,:,i];
#     describe(sm);end

function dfr(x::AbstractString)
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

function dfr(x::Regex)
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

function readalloutput()
    """
    reads timeseries and stores to Vector{DataFrame}
    reads NetCDFs and stores to Vector{Any}
    """
    cwd = "."
    dfs=loadalldfs(cwd)
    ncs=readallras(cwd)
    return(dfs,ncs)
end

# dfs,ncs = readalloutput()
# # filterplot("qg",dfs)
# # filterplot("rad",ncs)


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

# df=readf(r"qoutjl")
# df|>tpjs
# qpl(df)
#r"sch"i|>tpjs

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

#rsqgrep()

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

# function nsegrep()
#     """
#     ...sort try... but... hmm
#     When j=1, the modified NSeff is not inflated by the squared values 
#     of the differences, because the squares are replaced 
#     by absolute values.
#     """
#     path = pwd()
#     #match = Dict(SubString{String}, Vector{String})
#     match = []
#     files = glob(r"_output.txt|_outputjl") #non-recurse
#     #@printf("Searching for values > 0.3 in files matching pattern %s\n", path)
#     for file in filter(file -> endswith(file, "_output.txt"), files)
#         output = DelimitedFiles.readdlm(file, '\t', String)
#         DataFrame(output[4:end,:],:auto)


#         fn = first(split(file, "_qout"))
#         #append!(match,fn)
#         m = [fn,Grep.grep(r"mNSE", output)]
#         append!(match,m)
#     end
#         if !isempty(match)
#             split_elements = [split(string(item[1])) for item in match]
#             sort!(split_elements, by = x -> parse(Float64, x[end]);rev=true)
#             # Combine the sorted elements into a string
#             #combined_string = join([join(item, "     ") for item in split_elements], "\n")
#             combined_string = [join(item, "     ") for item in split_elements]
#             for line in combined_string
#                 line = strip(line)  # remove leading and trailing whitespace
#                 line = join(split(line), " ")  ##remove inner whitespaces
#                 printstyled(rpad("$fn:", 55), lpad("$line\n", 10), color = :green)
#             end
#         end
#     #end
# end

# z=glob(r"_output.txt|_outputjl")|>first
# m=Grep.grep(r"mNSE", DelimitedFiles.readdlm(z, '\t', String))
# m=only(m)
# parse(Float64, split(m)[end])



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

#nsegrep()
# #filterplot("ar",dfs)
# xdd = mall(dfs)
# dfp(xdd)

# function filterplot(regex::AbstractString,dfs::Vector{Any})
#     "selects first match and plots..."
#     df = dfs[map(n->occursin(Regex(regex,"i"),n),
#     map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
#     )] |> first
#     dfp(df)
# end
# waread(file_path|>first)

# v=[]
# for z in file_path;
#     df = waread(z);
#     dropmissing!(df)
#     observed  = df[:,1]
#     simulated = df[:,2]
#     kge_value = kge1(simulated,observed)
#     nse_value = nse(simulated,observed)
#     #nm = basename(file_path)
#     nm="tst"
#     println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
#     printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
#     push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value))
#     v = DataFrame(v)    
# end

# PrettyTables.pretty_table(v)
# xd = file_path|>lastbefore|>readdf
# ftp(xd)
# tpjs(xd)

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

#@ncrm #works

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

# function findindf(df::DataFrame,col,x)
#     #filter(row -> contains(col,x), df)
#     filter(row -> occursin(Regex(x,"i"),col), df)
# end
# #filter(row -> contains(row.name,"_qout"), ds)
# findindf(df,"name","E")
# x="Et"
# col="name"
# filter(row -> occursin(Regex(x,"i"),df[!,col]), df)

function findindf(df::DataFrame, x)
    """
    like Grep.grep("x",df)
    """
    filter(row -> any(occursin(Regex(x, "i"), 
        string(value)) for value in row), 
            eachrow(df))
end

# findindf(df,"Et")
# findindf(df,"15")
# findindf(df,"-")
# Grep.grep(r"Et",df.name)
# Grep.grep(r"Et",df)
# Grep.grep("Et",df)


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

end ##end of module


#xdd = wa.xread(pt)
#wa.dfp(pt)
#include(raw"C:\Users\Public\Documents\Python_Scripts\julia\wa.jl")
#wa.aplot(xdd)