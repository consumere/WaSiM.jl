#functions
# using DataFrames, CSV, Statistics, Dates, PlotlyJS
#using Printf, DataFrames, CSV, Statistics, Dates, Rasters, Distributions, StatsPlots, PlotlyJS
# no JS plots
using Printf, DataFrames, CSV, Statistics, Dates, Rasters, Distributions, StatsPlots
#using Query
function setup()
#    include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/functions.jl")
    thisfile="/mnt/c/Users/Public/Documents/Python_Scripts/julia/functions_v1.9.jl"
    include(thisfile)
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

list(x) = Any[i for i ∈ x]

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

function ddense(df::DataFrame)
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df density(cols(s), legend = :topright)
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

function dfp(df::DataFrame;)
    #ti = DataFrames.metadata(df)|>only|>last|>basename 
    ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
    end
    if (any(x->occursin("year",x),names(df)))
        s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
        @df df Plots.plot(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    @df df Plots.plot(:date,cols(s),
    #legend = :bottomright, 
    legend = :outertopright,
    title=ti)
    end
end

# function f(x, y, z);x + y + z;end
# args = (1, 2, 3)
# f(args...) #6 splatting

function dfp(df::String)
    df=readdf(df)
    o = DataFrames.metadata(df)|>collect
    ti = basename(o[1][2])
    if (any(x->occursin("year",x),names(df)))
        s = Symbol.(filter(x->!occursin("year",x),names(df)))
        @df df Plots.plot(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot(:date,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
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
        @df df Plots.plot(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot(:date,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
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
        @df df Plots.plot(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot(:date,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
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
        @df df Plots.plot!(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    @df df Plots.plot!(:date,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    end
end

function dfp!(df::String)
    df=readdf(df)
    o = DataFrames.metadata(df)|>collect
    ti = basename(o[1][2])
    if (any(x->occursin("year",x),names(df)))
        s = Symbol.(filter(x->!occursin("year",x),names(df)))
        @df df Plots.plot!(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot!(:date,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
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
        @df df Plots.plot!(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot!(:date,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
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
        @df df Plots.plot!(:year,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot!(:date,cols(s),
#legend = :bottomright,
legend = :outertopright, title=ti)
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

function waread(x::AbstractString)
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
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
    metadata!(df, "filename", x, style=:note);
end

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
    ms=["-9999","lin","log","-9999.0"]
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
# DataFrames.combine(groupby(df,:year),:=>sum)
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
    ti = DataFrames.metadata(df)|>only|>last|>basename
    s = names(df)[1:2]
	t2 = string.(ti,"\n",s[1],"|",s[2],ti)
    StatsPlots.plot( 
	qqplot(df[!,1],df[!,2], qqline = :fit), 
	qqplot(Cauchy,df[!,2]), 
	qqnorm(df[!,2], qqline = :R),
    title = t2)
end

qqp=qpl

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

function vio(mm::Regex)
    """
    vioplot wasim timeseries 
    """
    df = readdf(mm)
    str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    o = DataFrames.metadata(df)|>collect
    ti = basename(o[1][2])
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
    @df df Plots.plot(:date,cols(s),legend = :outertopright, title=ti)
    #legend = :bottomright,
    end
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
# dfm = DataFrames.combine(groupby(df, :month), y .=> sum .=> y);
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

##unlist df
#[eachrow(observed)...]
#vec(Matrix(observed))

function nse(df::DataFrame)
    observed, simulated = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
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

function kge_fread()
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

function kge_df(ext::Regex)
    """
    should be recursive, with ext funcs 
    """
    m=rglob(ext)
    files=filter(x->!occursin(r"xml|txt|mon|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",x),m)
    v = []
    for file in files
        if isfile(file)
            # dd = CSV.read(file,DataFrame,missingstring="-9999",delim="\t")
            # observed  = dd[:,5]
            # simulated = dd[:,6]
            dd  =   readdf(file)
            observed  = dd[:,end-2]
            simulated = dd[:,end-1]
            kge_value = kge2(observed, simulated)
            nse_value = nse(observed, simulated)
            nm = basename(file)
            println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
            printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
            push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
            v = DataFrame(v)
        end
    end
    return(v)
end


function kge2(observed, simulated)
    r = cor(observed, simulated)
    α = std(simulated) / std(observed)
    β = mean(simulated) / mean(observed)
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

#ds = kgerec("qout")
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

function rp(x::Regex)
    v = rglob(x)
    file = v[broadcast(x->endswith(x,"nc"),v)]|>first;
    ts=read(Raster(file,missingval=0))
    contourf(ts; c=cgrad(:thermal),size=(1200, 800))
end

function rp(file::AbstractString, lyr::Int)
    ts=read(Raster(file,missingval=0))
    x = ts[t=lyr]
    contourf(x; c=cgrad(:thermal),size=(1200, 800))
end

function rp(x::Raster)
    contourf(x; c=cgrad(:matter),size=(1200, 800),
    xlabel="",ylabel="")
end

mask_trim(raster, poly) = Rasters.trim(Rasters.mask(raster; with=poly); pad=10)

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

function rplot!(file::Regex)
	#v::Vector{String} = readdir();
    #file=r"win"
    #v = occursin(file,filter(x->endswith(x,"nc"),readdir()))
    z = first(filter(x->endswith(x,"nc") & occursin(file,x),readdir()))
    # v = v[broadcast(x->endswith(x,"nc"),v)];
    # z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
    # filter(x -> occursin(Regex(file,"i"),x), readdir())
    xr = read(Raster(z;
        crs=EPSG(25832),
        missingval=0))
	Plots.plot!(xr;
    c=cgrad(:thermal),
    xlabel="",
    ylabel="",
    size=(1200*.66, 800*.66))
end

function rpm(x::Regex;msk::Float64,gt=false)
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

function cpl(file::Raster)
    #x=Rasters.rebuild(file;missingval=-9999)
    x=Rasters.rebuild(file;missingval=0)
    #x=x[t=1]
    Plots.contourf(x; 
    c=cgrad(:matter),
    xlabel="",
    ylabel="")
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
    #describe(x)
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
    cwd = pwd()
    dirs = readdir(cwd)
    printstyled("$cwd\n",color=:magenta)
    for dir in dirs
        if isdir(dir)
            @printf("%-5s | ","$dir");
        end
    end
end

function fd(cwd::AbstractString)  
    dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
    for dir in dirs
        if isdir(dir)
            @printf("%-8s |","$dir");
        end
    end
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

dfl = lplot

function dfl!(regex::AbstractString,dfs::Vector{DataFrame})
    "selects first match and plots..."
    df = dfs[map(n->occursin(Regex(regex,"i"),n),
         map(x->basename(only(DataFrames.metadata(x))[2]),
         dfs))] |> first
         ln = Symbol.(filter(x->!occursin("date",x),names(df)))
         nm = propertynames(df)[1:end-1];
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
        ti = raw""
    end
         @df df Plots.plot!(:date,cols(ln),yaxis=:log,title=ti)  
end

function dfl!(df::DataFrame)
    ti = try
        DataFrames.metadata(df)|>only|>last|>basename
    catch
        @warn "No basename in metadata!"
    ti = raw""
    end
    ln = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot!(:date,cols(ln),yaxis=:log,title=ti)     
end

function dfl!(df::String)
    df=readdf(df)
    ti = try
        DataFrames.metadata(df)|>only|>last|>basename
    catch
        @warn "No basename in metadata!"
    ti = raw""
    end
    ln = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot!(:date,cols(ln),yaxis=:log,title=o)     
end

function dfl!(x::Regex)
    df=readdf(x)
    ti = try
        DataFrames.metadata(df)|>only|>last|>basename
    catch
        @warn "No basename in metadata!"
    ti = raw""
    end
    ln = Symbol.(filter(x->!occursin("date",x),names(df)))
    @df df Plots.plot!(:date,cols(ln),yaxis=:log,title=o)     
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

function getnames(rg::String,dfs::Vector)
	nms = [];
	for i in dfs;	
		x=collect(DataFrames.metadata(i))[1][2]|>basename
		push!(nms, x)
	end
	return(filter(x->occursin(Regex(rg,"i"),x),nms))
end

function getnames(dfs::Vector)
	nms = [];
	for i in dfs;	
		x=collect(DataFrames.metadata(i))[1][2]|>basename
		push!(nms, x)
	end
	return(nms)
end
function getnames(ncs::Vector{Raster})
	x=map(x->name(x),ncs)
	return(x)
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
readalldfs = loadalldfs

function alldf(path::Regex)
    """
    like loadalldfs but with txt match
    """
    v = readdir();
    v = v[broadcast(x->!endswith(x,"nc"),v)];
    files = v[(broadcast(x->occursin(path,x),v))];
    dfs::Vector{DataFrame} = []
    for file in files
        if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
           file_path = file
	   println("reading ",file_path,"...")
	   p1 = readdf(file_path)
	   push!(dfs, p1)
        end
    end
    return(dfs)
end

function readalldf(pattern::AbstractString)
    """
    reads recursivley and stores to Vector{DataFrame}
    """
    v::Vector{String} = recursive_glob_prfx(pwd(),pattern)
    dfs::Vector{DataFrame} = map(x->readdf(x),v)
    return(dfs)
    pts = map(x->dirname(x),v)
    nm=getnames(dfs)
    # map((x, y) -> joinpath(x,y), v, nm)
    # broadcast((x, y) -> joinpath(x,y), v, nm)
    # #(+).([1, 2, 3], [4, 5, 6])
    # (joinpath).(v, nm)
    # #d = (joinpath).(v, nm)
    #map(x->x*" stored from",nm)
    d = DataFrame(hcat(v,nm), :auto)
    rename!(d,1=>"path",2=>"filename")
    printstyled(d , color=:green    )
    #printstyled("$d\n" , color=:green    )
end


function loadalldfs(path::AbstractString)
    files = readdir(path)
    dfs = DataFrame[]
    #nms = []
    for file in files #&& occursin(Regex(prefix),file)
        if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|year|yrly|nc|png|svg|zip|tar",file))
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

# function regand(v::String,xv::String)
#     needle=join(xv,"+.*");
#     z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
# return(z)
# end

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



function hometeo()
    cd("/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/")
    println("you are here: ",pwd())
end

function homreg()
    cd("/mnt/d/Wasim/regio/out/");
    println("you are here: ",pwd())
    fd()
end

function hombr()
    cd("/mnt/d/Wasim/Tanalys/DEM/brend_fab/out");
    println("you are here: ",pwd())
    fd()
end

function writewa(file::AbstractString, df::DataFrame)
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

wawrite = writewa

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



function vgr(snippet::AbstractString)
    owd=pwd()
    cd("/mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions")
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

function vgjl(snippet::AbstractString)
    owd=pwd()
    cd("/mnt/c/Users/Public/Documents/Python_Scripts/julia")
    files = filter(file -> endswith(file, ".jl"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if Base.contains(line,snippet)
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
    owd=pwd()
    cd("/mnt/c/Users/Public/Documents/Python_Scripts")
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


# function vgctl(snippet::AbstractString)
#     """
#     hint. vgctl("set \$TS")
#     """
#     owd=pwd()
#     nwd="/mnt/d/Wasim/regio/control/"
#     nwd2="/mnt/d/temp/saale/control/"
#     cd(nwd)
#     println("greps from *ctl from  \n$nwd \n$nwd2 ...")
#     files = filter(file -> endswith(file, ".ctl"), readdir())
#     for file in files
#         open(file) do f
#             counter = 0 # Zähler initialisieren
#             for line in eachline(f)
#                 counter += 1 # Zähler erhöhen
#                 if contains(line,snippet)
#                     printstyled("$counter:$nwd",color=:light_red) 
#                     printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
#                     printstyled("$line\n",color=:green,bold=true) 
#                 end
#             end
#         end
#     end
#     cd(nwd2)
#     files = filter(file -> endswith(file, ".ctl"), readdir())
#     for file in files
#         open(file) do f
#             counter = 0 # Zähler initialisieren
#             for line in eachline(f)
#                 counter += 1 # Zähler erhöhen
#                 if contains(line,snippet)
#                     printstyled("$counter:$nwd2",color=:light_red) 
#                     printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
#                     printstyled("$line\n",color=:green,bold=true) 
#                 end
#             end
#         end
#     end
#     cd(owd)
# end

#vgctl("rcm200")

function vgctl(snippet::AbstractString)
    """
    hint. vgctl("set \$TS")
    """
    owd = pwd()
    nwd_list = [
        "/mnt/d/Wasim/regio/control/",
        "/mnt/d/temp/saale/control/",
        "/mnt/d/Wasim/Tanalys/DEM/brend_fab/control/",
        "/mnt/d/Wasim/regio/control/",
    ]
    
    for nwd in nwd_list
        cd(nwd)
        println("greps from *ctl from $nwd...")
        files = filter(f -> endswith(f, ".ctl"), readdir())
        for file in files
            counter = 0 # Zähler initialisieren
            open(file) do f
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line, snippet)
                        printstyled("$counter:$nwd", color=:light_red) 
                        printstyled("$file:\t", color=:light_magenta, underline=true, blink=false, bold=true) 
                        printstyled("$line\n", color=:green, bold=true) 
                    end
                end
            end
        end
    end
    
    cd(owd)
end



function vgrepl(snippet::AbstractString)
    """
    greps from repl_history
    """
    #owd=pwd()
    #cd("C:/Users/Public/Documents/Python_Scripts")
    file = raw"/home/ubu/.julia/logs/repl_history.jl"
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

#vgrepl("Pkg")

#"mask"|>vgpy
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

# #using GeoArrays
# function gplot(r::AbstractString)
#     geoarray = GeoArrays.GeoArray(ArchGDAL.readraster(lk))
#     geoarray|>plot
# end

function rglob(prefix::AbstractString)
    rootdir::String = "."
    results::Vector{String} = []
    #results = [] #Vector{Any}
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
    rootdir::String = "."
    results::Vector{String} = []
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
        legend = :outertopright,
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
    @df df Plots.plot!(:date,cols(s),
    legend = :outertopright)
    #legend = :bottomright)
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

function dfpall(x::Regex)
    "reads, reduces + merges by date and plots"
    files = rglob(x)
    dfs = DataFrame[]
    for file in files
        if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
           file_path = file
	   println("reading ",file_path,"...")
	   p1 = readdf(file_path)
       
       for x in 1:size(p1,2)-1
        rename!(p1,x=>basename(file_path)*names(p1)[x])
       end
       
       #rename(p1,Not(Cols("date")))
	   push!(dfs, p1)
        end
    end
    df = reduce((left, right) -> 
      innerjoin(left, right, on = :date,makeunique=true), 
      dfs)
    ##to preserve column order and names (date at last position)
    df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
    y = filter(x->!occursin("date",x), names(df))
    s = map(y -> Symbol(y),y)
    @df df Plots.plot(:date,
            cols(s),
            legend = :outertopright)
            #legend = false)
end

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

#r"sb"|>dfpall
# r"wind"|>ldfpall
# r"qg"|>dfpall

function dfpall(files::Vector{Any})
    "reads, reduces + merges by date and plots"
    #files
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
      innerjoin(left, right, on = :date,makeunique=true), 
      dfs)
    ##to preserve column order and names (date at last position)
    df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
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
    ##to preserve column order and names (date at last position)
    df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
    y = filter(x->!occursin("date",x), names(df))
    s = map(y -> Symbol(y),y)
    @df df Plots.plot(:date,
            cols(s),
            #yaxis = :log,
            legend = :bottom)
end


function mall(files::Vector{String})
    "reads, reduces + merges by date"
    #files
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
      innerjoin(left, right, on = :date,makeunique=true), 
      dfs)
    ##to preserve column order and names (date at last position)
    df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
    return(df)
end

#v = r"qg"|>rglob
#vm=mall(v)
#bardfm(vm)

function mall(files::Vector{DataFrame})
    "reduces + merges by date"
    df = reduce((left, right) -> 
      innerjoin(left, right, on = :date,makeunique=true), 
      files)
    ##to preserve column order and names (date at last position)
    df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
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


# function p()
#     return(pwd())
# end

function yrsum(x::String)
    df = readdf(x)
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    #ti=DataFrames.metadata(df)|>only|>last|>basename
    df[!, :year] = year.(df[!,:date]);
    df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
    return(df_yearsum)
end

function yrsum(x::DataFrame)
    df = x
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    #ti=DataFrames.metadata(df)|>only|>last|>basename
    df[!, :year] = year.(df[!,:date]);
    df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
    return(df_yearsum)
end

function yrmean(x::String)
    df = readdf(x)
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    #ti=DataFrames.metadata(df)|>only|>last|>basename
    df[!, :year] = year.(df[!,:date]);
    df_yearsum = DataFrames.combine(groupby(df, :year), y .=> mean .=> y);
    return(df_yearsum)
end

function yrmean(x::DataFrame)
    df = x
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    #ti=DataFrames.metadata(df)|>only|>last|>basename
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
#     ti=DataFrames.metadata(df)|>only|>last|>basename
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
#     y = filter(x->!occursin("date|year",x),names(df))
#     s = map(y -> Symbol(y),y)
#     ti=DataFrames.metadata(df)|>only|>last|>basename
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
#         y = filter(x->!occursin("date|year",x),names(df))
#         s = map(y -> Symbol(y),y)
#         ti=DataFrames.metadata(df)|>only|>last|>basename
#         df[!, :year] = year.(df[!,:date]);
#         df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
#         @df df_yearsum Plots.plot(:year,
#             cols(s),
#             legend = :topright, 
#             title=ti,
#             seriestype=:bar)
# end

function bardf(x::String)
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
    y = filter(x->!occursin("date|year",x),names(df))
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

#r"sb05"|>dfbarjs

# df = readdf(r"sb05")
#df = df[!,Cols(r"date","_5","_139")]
# nrows=size(df)[2]-1
# fig = PlotlyJS.make_subplots(shared_xaxes=true, shared_yaxes=true)

# function contour3()
# trace = PlotlyJS.contour(x=df[:,2],y=df[:,3],z=df[:,4])
# layout = PlotlyJS.Layout(title="Setting the X and Y Coordinates in a Contour Plot")
# PlotlyJS.plot(trace, layout)
# #contour(x=x, y=y, z=z)


# for i in 1:nrows;
#     PlotlyJS.add_trace!(fig, 
#     PlotlyJS.contour(x=df.date, y=df[:,i],
#     name=names(df)[i]))
# end

    
function dfcntjs(df::Regex;)
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
        PlotlyJS.contour(x=df.date, y=df[:,i],
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

# r"sb05"|>dfcntjs


function monsum(x::String)
    df = readdf(x)
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

function rplot(x::Regex;missval="")
    """
    rplot(x::Regex)
    reads first match of regex wasim nc
    """
    rootdir="."
    results = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (occursin(x,filename)) && occursin(r".nc",filename) && 
                (!occursin(r"txt|yrly|png|svg|grd",filename)) 
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    file = first(results)
    pt = normpath(file)
    printstyled("now plotting $pt ...\n",color=:green)
    #missval=0;
    if isempty(missval) | (missval==0)
        missval=0
        @warn "missingval set to zero!"
    end
    xr = read(Raster(file;crs=EPSG(25832),missingval=missval))
        Plots.plot(xr;
        c=cgrad(:thermal),
        xlabel="",
        ylabel="",
        size=(1200*.66, 800*.66))   
end

# rplot(r"qd";missval=-0.1)
# rplot(r"rain")
# rplot(r"rain";missval=-1)


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


function dfsjs(dfs::Vector{DataFrame};save="",fact=.7,logy=false)
    "plots and adds"
    df = dfs[1]
    nrows=size(df)[2]-1
    ti = try
        DataFrames.metadata(df)|>only|>last|>basename
    catch
        @warn "No basename in metadata!"
        ti = raw""
    end
    p = PlotlyJS.make_subplots(
        shared_xaxes=true, 
        shared_yaxes=true    
        );
    for i in 1:nrows;
        PlotlyJS.add_trace!(p, 
        PlotlyJS.scatter(x=df.date, y=df[:,i],
        name=ti*"_"*(names(df)[i])));
    end
    for i in 2:length(dfs)
        nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
        println("adding $nm")
        #s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
        #@df dfs[i] Plots.plot!(:date,cols(s),label="$nm")
        df = dfs[i]
        nrows=size(df)[2]-1
        for i in 1:nrows;
            PlotlyJS.add_trace!(p, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=nm*"_"*(names(df)[i])));
        end
        #addtraces!(p)
    end
    #fact,logy = .60,0 #hardcoded
    if logy == true
        PlotlyJS.relayout!(p,
        template="seaborn",
        yaxis_type="log",
        height=600*fact,width=900*fact)
        #,title_text="Series of "*ti)
    else
        PlotlyJS.relayout!(p,
        template="seaborn",
        height=600*fact,width=900*fact)
        #,title_text="Series of "*ti)
    end
    display(p)
    if !isempty(save) 
        PlotlyJS.savefig(p,save*".png")
        printstyled("$save saved as $save*.png! \n",color=:green)
    end
end


# pwd()
# vv = rglob(regand("rgex","44"))
# dfpjs("./w_readgrids/rgexggb.v44.2021")
# vv = filter(x->!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",x),vv)
# dfs = readall(vv)
# dfsjs(dfs)
# dfsjs(dfs;logy=true)
#dfsjs(dfs;logy=true,save="test")
# vgjl("theme")
# vgjl("template")
# vgjl("climex")

function jread(jsonfile::AbstractString)
    """
    reads GeoJSON, returns DataFrame
    """
    #using GeoJSON
    geojson_file=jsonfile
    jsonbytes = read(geojson_file) # read the geojson file as bytes
    fc = GeoJSON.read(jsonbytes)
    cr = DataFrame(hcat([(x) for x in fc.geometry]...),:auto)|>permutedims
    nm = DataFrame(hcat([(x) for x in fc.Names]...),:auto)|>permutedims
    ht = DataFrame(hcat([(x) for x in fc.ht]...),:auto)|>permutedims
    df = hcat(nm,ht,cr,makeunique=true)
    rename!(df,["Location","Height","X","Y"])
    metadata!(df, "filename", jsonfile, style=:note);
    return(df)
end

function agjson(jsonfile::AbstractString)
    """
    reads json and transforms points...
    but AG transformation is wrong :(
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
# pt="/mnt/d/ClimateExplorer/precipitation/pre_ogr.geojson"
# we=agjson(pt)


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

function ll()
    readdir()
end

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
# @edit pwd()
# ##from docs:
# readdir(; join::Bool=false, sort::Bool=true) =
#     readdir(join ? pwd() : ".", join=join, sort=sort)

# ll(;reg=true) |> third|>rp
# ll(;reg=true) |>x->x[end-4]|>facets

# ll(;reg=true) |>y -> filter(x->occursin(r"temp",x),y)

# ll(;reg=true) |>y -> filter(x->occursin(r"temp",x),y)|>only|>facets

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
    plot(td;
    title = "Tdiff[mm] of "*ti,
    c=cgrad(:matter))
    ##for overwriting force
    #write("tdiff.nc",td;force=true)
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
    hier müsste dann push!(p,...) and @layout plot kommen...
    """
    path = pwd()
    files = readdir(path)
	for file in files
           i = joinpath(path, file)
           if isfile(i) && occursin(Regex(ext,"i"),file) && endswith(file, ".nc")
            r=read(Raster(i,missingval=0,mappedcrs=EPSG(25832)));
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
            fact=.75
            if (r.dims[end]|>length == 1)
                @warn("only one layer available...")
                p=Plots.plot(r;
                # xlabel="",
                # ylabel="",
                #c=cgrad(:thermal),
                c=cgrad(:matter),
                size=(1200*fact, 800*fact));
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
                size=(1200*fact, 800*fact));
                display(p)
           #savefig(p,outname)
           #println(basename(outname)," saved!");
        end
   end
end

# x="/mnt/d/Wasim/Goldbach/revision/fab150/pestpp/v5/pestout/old/windfab150.2018.nc"
# x="/mnt/d/Wasim/Goldbach/revision/fab150/pestpp/v5/pestout/old/tsoilfab150_stack.2018.nc"
# r=read(Raster(x,missingval=0,mappedcrs=EPSG(25832)))
# r.dims[end]|>length == 1
# facets(x)
# vgjl("cgrad")
# vgjl("layout")
function facets(ext::Raster)
    """
    like stackplot, but for interactive view
    """
    r=ext
    fact=.75
    if (r.dims[end]|>length == 1)
        @warn("only one layer available...")
        p=Plots.plot(r;
        c=cgrad(:matter),
        size=(1200*fact, 800*fact));
        display(p)
    else
        @warn("subsetting first layer...")
        ee = Int(r.dims[3][end])
        rn = r[t=2:ee];    #subset till end
        p=Plots.plot(rn;
        xlabel="",
        ylabel="",
        c=cgrad(:thermal),
        size=(1200*fact, 800*fact));
        display(p)
   end
end


function facets(ext::Regex)
    """
    like stackplot, but for interactive view
    """
    grids = filter(x -> occursin(ext,x), readdir())
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
            fact=.75
            if (r.dims[end]|>length == 1)
                @warn("only one layer available...")
                p=Plots.plot(r;
                c=cgrad(:matter),
                size=(1200*fact, 800*fact));
                display(p)
            else
                @warn("subsetting first layer...")
                ee = Int(r.dims[3][end])
                rn = r[t=2:ee];    #subset till end
                p=Plots.plot(rn;
                xlabel="",
                ylabel="",
                c=cgrad(:thermal),
                size=(1200*fact, 800*fact));
                display(p)
        end
   end
end

dfread = readdf

##
second(x) = x[2]
third(x) = x[3]
lastbefore(x) = x[end-1]
# geht auch:
# v |> x -> x[2]
# v |> x -> x[Not(2)]
# rglob(r"pr")|>lastbefore     
#rglob(r"pr")|>lastbefore|>dfread|>describe   
# rglob(r"pr")|>lastbefore|>yrmean

function reorder!(df::DataFrame)
    """
    date to last position
    """
    df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
#    return(df)
end


function dpr(x)
    df=readdf(x)
    v = (names(df[!,1:2]))
    a = reshape(v, 1, 2)
    Plots.plot(df.date,[df[!,2], df[!,1]], 
    #label=["Simulated" "Observed"], xlabel="Date", ylabel="Value")
    label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    annotate!(last(df.date), 0.85*maximum(df[!,1]),
    text("R² = $r2", 10, :black, :right))
    kge = round(kge2(df[!,1], df[!,2]), digits=2)
    annotate!(
        #:topright,
        last(df.date), 0.95*maximum(df[!,1]),
    text("KGE = $kge", 10, :black, :right))
end


function dpr!(x)
    df=readdf(x)
    v = (names(df[!,1:2]))
    a = reshape(v, 1, 2)
    Plots.plot!(df.date,[df[!,2], df[!,1]], 
    label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    annotate!(last(df.date), 0.85*maximum(df[!,1]),
    text("R² = $r2", 10, :black, :right))
    kge = round(kge2(df[!,1], df[!,2]), digits=2)
    annotate!(
        #:topright,
        last(df.date), 0.95*maximum(df[!,1]),
    text("KGE = $kge", 10, :black, :right))
end

function qall_num()
    files = rglob("qgko")
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

function qall()
    files = rglob("qgko")
    for file in files
        # Load the file into a DataFrame
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
            #dout.basin = map(x -> parse(Int, x), dout.basin)
            #dout.basin = map(x -> parse(string, x), dout.basin)
            #dout.basin = map(x -> rename(x => "B_"*x), dout.basin)
            dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
            return(dout)
        catch
            @warn("error! ")
            # Skip files that can't be loaded as a DataFrame
            continue
        end
    end
end


function theplot(x::AbstractString)
    df = DataFrame(CSV.File(x))
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    ndf = df[!,Not(1:4)]
    rename!(ndf, [:Observed, :Simulated, :Date])
    overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    r2 = overall_pearson_r^2
    nse_score = nse(ndf[!, :Observed], ndf[!, :Simulated])
    kge_score = kge2(ndf[!, :Observed], ndf[!, :Simulated])
    #ti = "Time Series of $(uppercase(first(split(basename(x), '-'))))"
    ti = first(split(basename(x),"_"))
    #subs = "Pearson r: $(round(overall_pearson_r, digits=2))\nPearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    #p = plot(title=[ti, subs], ylabel="[unit/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    p = plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
    plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    annotate!(
    #nrow(ndf), 0.95*maximum(ndf.Observed),
    :bottomright,
    text("$subs", 10, :black, :right))
    return p
end

ftp(z::AbstractString) = theplot(first(filter(x->occursin(Regex(z,"i"),x),filter(x->endswith(x,"qout"),readdir()))))

function theplot(x::DataFrame)
    ndf = x
    rename!(ndf, [:Observed, :Simulated, :Date])
    overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    r2 = overall_pearson_r^2
    nse_score = nse(ndf[!, :Observed], ndf[!, :Simulated])
    kge_score = kge2(ndf[!, :Observed], ndf[!, :Simulated])
    ti = collect(DataFrames.metadata(df))[1][2] |>basename
    subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    p = plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
    plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    annotate!(
    :bottomright,
    text("$subs", 10, :black, :right))
    return p
end


function linp(x::AbstractString)
    df = DataFrame(CSV.File(x))
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    ndf = df[!,Not(1:4)]
    rename!(ndf, [:Observed, :Simulated, :Date])
    overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    r2 = overall_pearson_r^2
    nse_score = nse(ndf[!, :Observed], ndf[!, :Simulated])
    kge_score = kge2(ndf[!, :Observed], ndf[!, :Simulated])
    ti = first(split(basename(x),"_"))
    subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    p = plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", legend=:topleft)
    plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
    plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    annotate!(
    :bottomright,
    text("$subs", 10, :black, :right))
    return p
end


function linp(x::DataFrame)
    ndf = copy(x)
    rename!(ndf, [:Observed, :Simulated, :Date])
    overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    r2 = overall_pearson_r^2
    nse_score = nse(ndf[!, :Observed], ndf[!, :Simulated])
    kge_score = kge2(ndf[!, :Observed], ndf[!, :Simulated])
    ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
    end
    subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    p = plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", legend=:topleft)
    plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
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


function dfr(x::AbstractString)
    """
    --- reader with fewer constrains ---
    no |> dropmissing 
    df[!,Cols(r"^Col|date")] |>dfp  
    """
    ms=["-9999","lin","log"]
    df::DataFrame = CSV.read(x,    DataFrame,    missingstring=ms,
    delim="\t",    normalizenames=true,
    drop=(i, nm) -> i == 4) 
    df.YY=map(x ->Int(x),df.YY);
    df.MM=map(x ->Int(x),df.MM);
    df.DD=map(x ->Int(x),df.DD);
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,4:end]
    df=df[!,Not(Cols(r"^Column"))] #drops names starting with Column, usually all missing
    DataFrames.metadata!(df, "filename", x, style=:note);
end


function dfr(file::Regex)
    """
    takes first match 
    --- reader with fewer constrains ---
    no |> dropmissing 
    df[!,Cols(r"^Col|date")] |>dfp  
    """
    ms=["-9999","lin","log"]
    # x = Regex(x,"i")|>glob|>first
    # x = first(glob(x))

    z = first(filter(x->!endswith(x,"nc") & occursin(file,x),readdir()))

    
    df::DataFrame = CSV.read(z,    DataFrame,    missingstring=ms,
    delim="\t",    normalizenames=true,
    drop=(i, nm) -> i == 4) 
    df.YY=map(x ->Int(x),df.YY);
    df.MM=map(x ->Int(x),df.MM);
    df.DD=map(x ->Int(x),df.DD);
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,4:end]
    df=df[!,Not(Cols(r"^Column"))] #drops names starting with Column, usually all missing
    DataFrames.metadata!(df, "filename", x, style=:note);
end


function dfp(files::Vector{Any})
    "reads, reduces + merges by date and plots"
    #files
    dfs = DataFrame[]
    for file in files
        if isfile(file) && (!occursin(r"xml|year|fzt|ftz|log|ini|wq|.yrly|.nc|.png|.svg",file))
           file_path = file
	   println("reading ",file_path,"...")
	   p1 = dfr(file_path)
       ##renamer
       for x in 1:size(p1,2)-1
         rename!(p1,x=>basename(file_path)*names(p1)[x])
        end
	   push!(dfs, p1)
        end
    end
    df = reduce((left, right) -> 
      innerjoin(left, right, on = :date,makeunique=true), 
      dfs)
    ##to preserve column order and names (date at last position)
    df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
    y = filter(x->!occursin("date",x), names(df))
    s = map(y -> Symbol(y),y)
    @df df Plots.plot(:date,
            cols(s),
            #yaxis = :log,
            #legend = :bottom)
            legend = false)
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
# #names(ncs[1])

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

#rp3("/mnt/d/temp/saale/output/thulba/v0/windsaale.2017.nc")
#rp3("/mnt/d/temp/saale/output/thulba/v0/sb05saale.2017.nc")


function globdf(x::Regex)
    """
    greps ts from current dir Regex
    """
    filter(file -> (occursin(x,file) & 
    (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
    ), readdir())
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
    df = CSV.read(x,    DataFrame,    
    missingstring=ms,
    delim="\t",    
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

println("you are here: ",pwd())
fd(pwd())
#set backend from gr() to 
#plotlyjs()


#r"qges"|>dfp
