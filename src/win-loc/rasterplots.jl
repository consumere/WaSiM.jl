# include("src/MyPackage.jl");
# using SnoopPrecompile, Preferences
# #set_preferences!(SnoopPrecompile, "skip_precompile" => ["PackageA", "PackageB"])
# SnoopPrecompile.verbose[] = true   # runs the block even if you're not precompiling, and print precompiled calls

#localARGS = ( if isempty(ARGS) ; ["arg1", "arg2"] ; else ARGS ; end )
#localARGS = ( if isempty(ARGS) ; ["arg1", "arg2"] ; else ARGS ; end )

using Base
mkdir("/mnt/d/wslconda/rainfarm/compcache")
pt="/mnt/d/wslconda/rainfarm/compcache"
cd(pt)
pwd()

#Base.compilecache("a3a2b9e3-a471-40c9-b274-f788e487c689")
id = Base.PkgId(Rasters)
#Base.compilecache("/home/ubu/.julia/compiled/v1.8/Rasters/3uKXQ_4rCK0.ji")
Base.compilecache(id)
Base.compilecache(Base.PkgId(Plots))

#lyr = isempty(ARGS) ? Int(1) : ARGS[2]
if length(ARGS[1]) == 0
	println("need args! <file><opt:layer>...")
    exit()
end

if length(ARGS[2]) == 0
	lyr=1
	println("skipping to layer 1...")
end

file=ARGS[1];
lyr=parse(Int,ARGS[2]);
#lyr = isempty(ARGS) ? Int(1) : ARGS[2]
file="/mnt/d/Wasim/regio/rcm/rcm.art_bfid"

# using GMT
# C = makecpt(0,5);
# imshow("forcing-moselle.nc?precip", proj=:guess, cmap=C, colorbar=true)

pt="/mnt/d/Wasim/regio/rcm/"
walkdir(pt)
#set_preferences!(SnoopPrecompile, "skip_precompile" => ["Rasters", "Plots"])
using Rasters,Plots
println("load",pwd()," ",file," ","an subset to Layer t=",lyr,"...")

#ts[t=lyr]|>plot
s=file
m=match(r".*[.]",s)
outfile = string(m.match,"png")

ts=Raster(file,missingval=-9999)
x = ts[t=lyr]
print("saving countour plot to",outfile,"...")

trim(x)|>plot
# x[Where(x)]

plotsize = (1600,800)
#p = contourf(x; dpi=300, size=(800, 400))
x=trim(x)
p = contourf(x; dpi=300, size=plotsize)

using Shapefile
# Load using Shapefile.jl
shapefile_name = "/mnt/d/Wasim/regio/rcm/ezg.shp"
ezg = Shapefile.Handle(shapefile_name)
#ezg|>plot
mask_trim(raster, poly) = trim(mask(raster; with=poly); pad=10)
ezg.shapes[end]
msk=mask_trim(ts, ezg.shapes[end])
plot(msk)

# default(show = true)
# #p = plot(x,size=plotsize)
# display(p)
#savefig(p, "t.png")
savefig(p, outfile)
print(outfile,"... saved! \n"); nothing


# Julia also provides * for string concatenation:
# julia> greet * ", " * whom * ".\n"
# "Hello, world.\n"


# if you only want the geometries and not the metadata in the DBF file
table = Shapefile.Table(shapefile_name)
geoms = Shapefile.shapes(table)
# whole columns can be retrieved by their name
table.rcm  # => Union{String, Missing}["Square with triangle missing", "Smaller triangle", missing]
# example function that iterates over the rows and gathers shapes that 
# meet specific criteria
function selectshapes(table)
    geoms = empty(Shapefile.shapes(table))
    for row in table
        if !ismissing(row.TestDouble) && row.TestDouble < 2000.0
            push!(geoms, Shapefile.shape(row))
        end
    end
    return geoms
end

# the metadata can be converted to other Tables such as DataFrame
using DataFrames
df = DataFrame(table)
#using GeoInterface
#GeoInterface.coordinates(Shapefile.shape(first(table)))

using ArchGDAL
const ag = ArchGDAL;
mys = ag.read(shapefile_name)

using GeoArrays
fn="/mnt/d/remo/cordex/eobs/out2.nc"
geoarray = GeoArrays.read(fn)
geoarray.crs
const ga = GeoArrays;
ga.coords(geoarray, [1,1])
#geoarray.EPSG = "4326"
epsg!(geoarray, 4326) 

plot(geoarray)

Plots.plotly(df[:,2])

using PlotlyJS
PlotlyJS.plot(cos, 0, 2π, mode="lines", Layout(title="cos(t)"))
PlotlyJS.plot(cos, 2π^2, 2π, mode="lines", Layout(title="okay"))


#large Data
const pj = PlotlyJS
using Random
N = 10000
pj.plot(scattergl(
    x=randn(N), y=randn(N), mode="markers",
    marker=attr(color=randn(N), colorscale="Viridis", line_width=1)
))

ts.data
ts[:,:,1]

df = DataFrame(ts[:,:,1])
rename!(df, ""=>"BK")
#columns(df)
#setproperty!(df[:,3],"BK")

#mx = parse.([Int],df.X)
#parse.(Int, df[!,"X"])
# mx = round.(df.X)
# my = round.(df.Y)
# my = df.BK
# mx=my
# pj.plot(df.X,my)
broadcast(x -> parse(string,x), df[!,"BK"] )
df

s = join(ENV)
for i in split(s,";");if occursin(r"JAVA",i) print(i); end;end  
for i in split(s);if occursin(r"QGIS",i) print(i); end;end  

s=ENV
#o_result = filter( p -> p[1], s)
#get(Dictionary_name, Key_name, Default Value)
get(s,"USER",0)
get(s,"PATH",1)
keys(s)
#values(Dictionary_name)
values(s)

# filter( p -> startswith("u"), s)
# filter( p -> contains("ubu"), s)
# filter( p ->  p[2], s)

for i in values(s);if occursin(r"QGIS",i) print(i); end;end  
for i in values(s); if occursin(r"QGIS",i) ; end;end  

using Grep
ENV |> grep(r"qgis"i)

filter(line -> occursin(r"qg"i, string(line)),ENV)

pat=r"QGIS"i
pat=r"xterm"
filter(el->occursin(pat, string(el)), ENV)



using PlotlyJS, Random
Random.seed!(42)
N = 100
random_x = range(0, stop=1, length=N)
random_y0 = randn(N) .+ 5
random_y1 = randn(N)
random_y2 = randn(N) .- 5

plot([
    scatter(x=random_x, y=random_y0, mode="markers", name="markers"),
    scatter(x=random_x, y=random_y1, mode="lines", name="lines"),
    scatter(x=random_x, y=random_y2, mode="markers+lines", name="markers+lines")
])

using DataFrames, CSV, Dates  
# df = DataFrame(CSV.File(path, header=skip, delim=sep,missingstring="-9999"));
# df = DataFrame(CSV.File(path, header=skip,missingstring="-9999"));
path="/mnt/d/Wasim/regio/out/r4/qgesrcm.v4.2017"
skip=1
df = CSV.read(path, DataFrame, header=skip,missingstring="-9999",)

dd = CSV.read(path, DataFrame, header=skip,missingstring=["-999","-9999", "--"],
comment="-")
df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), dd)

#comment=r"[a-z]")
#d = CSV.read(path, DataFrame,ignorerepeated=true,delim='\t')
#df = filter(d.YY => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), d)
#df = dd
#df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
#df = CSV.read(path, DataFrame, missingstring="-9999")
#df=df[in(start:stop).(year.(df.date)),Not(1:4)]
#df = df[Not(1:2),:]
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
#println(path,"\neval on:\n",(propertynames(df)[5:end-1]),"\n on years",start," to ",stop)
start,stop = 2015,2016
#df=df[in(start:stop).(year.(df.date)),Not(1:4)]
df=df[in(start:stop).(year.(df.date)),Not(1:4)]

plot([
    scatter(x=random_x, y=random_y0, mode="markers", name="markers"),
    scatter(x=random_x, y=random_y1, mode="lines", name="lines"),
    scatter(x=random_x, y=random_y2, mode="markers+lines", name="markers+lines")
])

for i in propertynames(df)[1:end-1];print(i,"\t");end

ps = []
#for i in propertynames(df)[1:end-1];
#df[:,1:end-1]
for i in 1:size(df)[2]-1;
    #print(i);end
    push!(ps,plot(
        [scatter(x=df.date, y=df[:,i], mode="lines", name="lines")]
        ))
    ;end
    
p = plot([scatter(x=df.date, y=df[:,1], 
mode="lines", name="Basin " * string(1))])

for i in 2:size(df)[2]-1;
    add_trace!(p, scatter(
        x=df.date, y=df[:,i], 
    name=i), row=1, col=1)
        ;end
    
p



using DataFrames, CSV, Dates  
#path="/mnt/d/Wasim/regio/out/r4/qgesrcm.v4.2017"
dir="/mnt/d/Wasim/regio/out/v8"
for (looproot, dir, filenames) in walkdir(dir)
    for i in filenames
        if occursin(r"^qo.*", i) println(i) end
    end
end


cd(dir)
skip=1
path="qout"
dd = CSV.read(path, DataFrame, #header=skip, 
missingstring=["-999","-9999", "--"],comment="-")
df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), dd)
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
df=df[:,Not(1:4)]
#println(path,"\neval on:\n",(propertynames(df)[5:end-1]),"\n on years",start," to ",stop)
#start,stop = 2015,2016
#df=df[in(start:stop).(year.(df.date)),Not(1:4)]
p = plot([scatter(x=df.date, y=df[:,1], 
mode="lines", name=string(propertynames(df)[1]))])
#name="Basin " * string(1))])
for i in 2:size(df)[2]-1;
    add_trace!(p, scatter(
        x=df.date, y=df[:,i], 
    name=string(propertynames(df)[i])))
    ;end
    #name="Basin " * string(i))) #, row=1, col=1)
relayout!(p, height=600, width=600, title_text=path)
p.plot

savefig(p,"myp.png")


# function topo_surface()
#     #z = Vector
# trace = surface(z=z)
#     layout = Layout(title="Mt. Bruno Elevation", autosize=false, width=500,
#                     height=500, margin=attr(l=65, r=50, b=65, t=90))
#     plot(trace, layout)
# end
# topo_surface()

function with_make_subplots1()
    p = make_subplots(rows=2, cols=1, shared_xaxes=true, 
    vertical_spacing=0.02)
    add_trace!(p, scatter(x=df.date, y=df[:,1], 
    mode="lines", name=string(propertynames(df)[1])), row=2, col=1)
    add_trace!(p, scatter(x=df.date, y=df[:,2], 
    mode="lines", name=string(propertynames(df)[2])), row=1, col=1)
    relayout!(p, title_text="Stacked "*path)
    p
end
with_make_subplots1()

function with_make_subplots2()
    p = make_subplots(rows=2, cols=1, shared_xaxes=true, shared_yaxes=true,
    vertical_spacing=0.08,
    subplot_titles=[string(propertynames(df)[1]) string(propertynames(df)[2]);]
    )
    add_trace!(p, bar(x=df.date, y=df[:,1],name=""), row=2, col=1)
    add_trace!(p, bar(x=df.date, y=df[:,2],name=""), row=1, col=1)
    relayout!(p, title_text="Stacked "*path)
    p
end
with_make_subplots2()


function wrap(path)
    dd = CSV.read(path, DataFrame, #header=skip, 
    missingstring=["-999","-9999", "--"],comment="-")
    df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), dd)
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
    df=df[:,Not(1:4)]
    p = make_subplots(rows=2, cols=1, shared_xaxes=true, shared_yaxes=true,
    vertical_spacing=0.08,
    subplot_titles=[string(propertynames(df)[1]) string(propertynames(df)[2]);]
    )
    add_trace!(p, bar(x=df.date, y=df[:,1],name=""), row=2, col=1)
    add_trace!(p, bar(x=df.date, y=df[:,2],name=""), row=1, col=1)
    relayout!(p, title_text="Stacked "*path)
    p
end
wrap("cloudrcm.v8.2017")

function wrapall(path)
    dd = CSV.read(path, DataFrame, #header=skip, 
    missingstring=["-999","-9999", "--"],comment="-")
    df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), dd)
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
    df=df[:,Not(1:4)]
    nrows=size(df)[2]-1
    st=[]
    for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
    p = make_subplots(rows=nrows, cols=1, 
    shared_xaxes=true, 
    shared_yaxes=false,
    vertical_spacing=0.08,
    #subplot_titles= st;
    )
    for i in 1:nrows;
            add_trace!(p, 
            bar(x=df.date, y=df[:,i],
            name=st[i]), 
            #name=""), 
            row=i, 
            col=1);
    end
    relayout!(p, title_text="Stacked "*path)
    p
end
#wrapall("albercm.v8.2017")
wrapall("wurzrcm.v8.2017")


path=fn

function pline(path::String,skip::Int)
    skip = isempty(skip) ? Int(1) : skip
    dd = CSV.read(path, DataFrame, header=skip, 
    missingstring=["-999","-9999", "--"],comment="-",delim="\t")
    df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), dd)
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
    df=df[:,Not(1:4)]
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
            scatter(x=df.date, y=df[:,i],
            name=st[i]),   row=i,     col=1);
    end
    relayout!(p, title_text="Series of "*basename(path))
    p
end

pline("wurzrcm.v8.2017")
fn="/mnt/d/Wasim/regio/out/v8/ts_avgrcm.v8.2017"
fn="/mnt/d/Wasim/regio/out/v8/rgexrcm.v8.2017"
pline(fn,6)
pline(fn,1)



##read all dfs by match

re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";
occursin(re,"123.")


path="/mnt/d/Wasim/regio/out/lowres/c5/loc2/qout_wlf"

dd = CSV.read(path, DataFrame, header=skip, 
           #drop = ,
            missingstring=["-999","-9999", "--"],comment="-",delim="\t")
            
            occursin(re,"123.")

rl=readlines(path)

xx = filter(el->occursin(re, (el)),rl )

for line in eachline(path)
    z=filter(el->occursin(re, string(el)),line)
    print(z)
end

for line in eachline(path)
    z=filter(el->occursin(re, string(el)),line)
    v=CSV.read(z, DataFrame, header=0)
end


df = DataFrame(rand(5, 5), :auto)


path="/mnt/d/Wasim/regio/out/lowres/c5/ts_avgrcm.c5.2017"
#file = CSV.File(path; dateformat="yyyy\tmm\tdd\thh",missingstring=["-999","-9999", "--"],comment="-",delim="\t")
file = CSV.File(path; dateformat="yyyy mm dd hh",missingstring=["-999","-9999", "--"],comment="-",delim="\t")
df = DataFrame(file)
df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
df=df[:,Not(1:4)]


path="/mnt/d/Wasim/regio/out/v8/tst"

df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan, isempty)), 
CSV.read(path,DataFrame,missingstring=["-999","-9999", " "],comment="-",delim="\t",ignorerepeated=true)
)

df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
df=df[:,Not(1:4)]
first(df)
# reading data without specific rows

#http://www.queryverse.org/Query.jl/stable/standalonequerycommands/#The-@dropna-command-1
using Query, DataFrames, CSV
df = CSV.read(path,DataFrame,header=1,missingstring=["-999","-9999"],comment="-",delim="\t",ignorerepeated=true)

df = CSV.File(path)
df = DataFrame(df)
q = df |> @dropna() |> DataFrame

dd = CSV.read(path,DataFrame,missingstring="-9999",
delim="\t",ignorerepeated=true)  |> @dropna() |> DataFrame

first(dd)


function baronly(df::DataFrame)
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
            bar(x=df.date, y=df[:,i],
            name=st[i]), 
            #name=""), 
            row=i, 
            col=1);
    end
    relayout!(p, title_text="Stacked "*path)
    p
end
baronly(df)

function fline(df::DataFrame)
    p = make_subplots(rows=2, cols=1, shared_xaxes=true, shared_yaxes=true,
    vertical_spacing=0.08,
    subplot_titles=[string(propertynames(df)[1]) string(propertynames(df)[2]);]
    )
    add_trace!(p, bar(x=df.date, y=df[:,1],name=""), row=2, col=1)
    add_trace!(p, bar(x=df.date, y=df[:,2],name=""), row=1, col=1)
    relayout!(p, title_text="Stacked "*path)
    p
end
fline(df)

function stackline(df::DataFrame)
p = plot([scatter(x=df.date, y=df[:,1], 
mode="lines", name=string(propertynames(df)[1]))])
for i in 2:size(df)[2]-1;
    add_trace!(p, scatter(x=df.date, y=df[:,i], name=string(propertynames(df)[i])));
end
relayout!(p, height=600, width=980, title_text=basename(path))
p.plot
end


function pline(path::AbstractString)
    df = CSV.read(path,DataFrame,missingstring="-9999",delim="\t",ignorerepeated=true)  |> @dropna() |> DataFrame
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
    df=df[:,Not(1:4)]
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
            scatter(x=df.date, y=df[:,i],
            name=st[i]),   row=i,     col=1);
    end
    relayout!(p,height=600*2,width=900*2,title_text="Series of "*basename(path))
    p
end
pline(path)



function aline(path::AbstractString)
    df = CSV.read(path,DataFrame,missingstring="-9999",delim="\t",ignorerepeated=true)  |> @dropna() |> DataFrame
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
    df=df[:,Not(1:4)]
    nrows=size(df)[2]-1
    st=[]
    for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
    p = plot([scatter(x=df.date, y=df[:,1],mode="lines", name=string(propertynames(df)[1]))])
    for i in 2:size(df)[2]-1;
        add_trace!(p, scatter(x=df.date, y=df[:,i], name=string(propertynames(df)[i])));
    end
    relayout!(p,
    yaxis_type="log",
    height=600*1,
    width=900*1,
    title_text="Series of "*basename(path))
    #Layout(xaxis_type="log")
    p
end
aline(path)



using PlotlyJS, CSV, DataFrames
df = dataset(DataFrame, "iris")
features = [:sepal_width, :sepal_length, :petal_width, :petal_length]
plot(df, dimensions=features, color=:species, kind="splom")

df = CSV.read(path,DataFrame,missingstring="-9999",
delim="\t",ignorerepeated=true)  |> @dropna() |> DataFrame
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
df=df[:,Not(1:4)]

# #features = [propertynames(df[:,1:end-1])]
# st=[]
# #for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
# for i in 1:size(df)[2]-1; push!(st,Symbol(propertynames(df)[i]));end
# #plot(df, dimensions=features, color=:date, kind="splom")
# plot(df, dimensions=st, kind="splom")

nrows=size(df)[2]-1
stackline(df)
##############################
using Plots
p = Plots.plot(1)
anim = @animate for x=0:0.1:5
    push!(p, 1, sin(x))
end
gif(anim)

##############################

function check_str(a)
    try
        parse(Float64,a)
        true
    catch
        false
    end
end

using BenchmarkTools
@btime check_str("60.0a") #264.600 μs
@btime occursin(re,"60.0a") #61.224 ns



using PlotlyJS, CSV, HTTP, DataFrames
# Read data from a csv
df = CSV.File(
    HTTP.get("https://raw.githubusercontent.com/plotly/datasets/master/api_docs/mt_bruno_elevation.csv").body
) |> DataFrame
z_data = Matrix{Float64}(df)

layout = Layout(
    title="Mt Bruno Elevation",
    autosize=false,
    scene_camera_eye=attr(x=1.87, y=0.88, z=-0.64),
    width=500, height=500,
    margin=attr(l=65, r=50, b=65, t=90)
)
plot(surface(
    z=z_data,
    contours_z=attr(
        show=true,
        usecolormap=true,
        highlightcolor="limegreen",
        project_z=true
    )
), layout)


# occursin(ext, string(file_path))
# pat=r"xterm"
# filter(el->occursin(pat, string(el)), ENV)
##################query################
function f_read(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    dfs = DataFrame[]
    for file in files
        file_path = joinpath(path, file)
        #if isfile(file_path) && endswith(file, ext)
        #if isfile(file_path) && contains(file, ext)
        if isfile(file_path) && occursin(ext, string(file))
            println(replace("reading  $file_path", "\\"  => "/")," now...")
            df = CSV.read(file_path,DataFrame,missingstring="-9999",comment="-",delim="\t")  |> @dropna() |> DataFrame
            #df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), dd)
            df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
            df=df[:,Not(1:4)]
#            println(replace("reading on $file_path", "\\"  => "/")," done!")
            push!(dfs, df)
        elseif isdir(file_path)
            dfs_in_subdir = f_read(file_path, ext)
        end
    end
    return(dfs)
end
directory_path = "/mnt/d/Wasim/regio/out/v8"
ext = "v8.2017"
#ext = "albe*v8"
#ext = "albe*v8"
f_read(directory_path,ext)


file_path = "/mnt/d/Wasim/regio/out/v8/albercm.v8.2017"
df = CSV.read(file_path,DataFrame,missingstring="-9999",comment="-",delim="\t")  |> @dropna() |> DataFrame
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
df=df[:,Not(1:4)]


function f(file)
    result = String[]
    for line in eachline(file)
        if startswith(line, "rout") || isempty(line)
            continue # skip the header and any empty lines
        elseif startswith(line, "lin")
            break # stop reading completely
        end
        push!(result, line)
    end
    return result
end

fn="/mnt/d/Wasim/regio/out/v8/qgkorcm.v8.2017"
z = f(fn)

first(z)
last(z)


for word in eachline(fn)
    if length(word) >= 24
        println(word)
    end
end


using Query, DataFrames, CSV, Dates, PlotlyJS

ms=["-9999","lin","log","LIN","LOG"]
missingstring=ms

df = CSV.read(fn,DataFrame;header=1,missingstring=ms) 
file = CSV.File(fn;delim="\t", missingstring=["-999", "NA"])

using CSV
data = """
code,age,score
0,21,3.42
1,42,6.55
-999,81,NA
-999,83,NA
"""
# by passing missingstring=["-999", "NA"], parsing will check each cell if it matches
# either string in order to set the value of the cell to `missing`
file = CSV.File(IOBuffer(data); missingstring=["-999", "NA"])

ms="-9999"
df = CSV.read(fn,DataFrame,missingstring=ms, #ignorerepeated=true,
delim="\t",comment="-",silencewarnings=true)  |> @dropna() |> DataFrame

df = CSV.File(fn,
missingstring=ms,
delim="\t",comment="-",silencewarnings=true)  |> DataFrame  |> @dropna() |> DataFrame

for word in CSV.File(fn)
    if isdigit(word)
        println(word)
    end
end



s = df[:,1] 
parse.(Float64, s)
parse.(Year, s )
d = parse.(Date, s )
d.year
#step=Year)

Date(s, dateformat="YYYY")


#https://github.com/JuliaData/SplitApplyCombine.jl
using SplitApplyCombine

# ]
# status --outdated

# using UpdateJulia
# update_julia()


using PkgCleanup;
PkgCleanup.artifacts()
using Pkg, Dates;
Pkg.gc(; collect_delay=Dates.Day(0))  

#As commented by Bogumił Kamiński, deleting the .julia/registries/General directory and running Pkg.update() worked.!!

Pkg.update()

#https://discourse.julialang.org/t/how-to-remove-rows-containing-missing-from-dataframe/12234/4

function pline(path::AbstractString)
    ms=["-9999","lin","log","LIN","LOG","--"] #comment="-",
    #df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",comment="-",ignorerepeated=true,silencewarnings=true,typemap=Dict(Int64=>String))  |> @dropna() |> DataFrame
    df = CSV.read(path,DataFrame,missingstring=ms,delim="\t",ignorerepeated=true,silencewarnings=true,typemap=Dict(String=>Int64))
    df = df[completecases(df), :]
    #df = filter( [2]=> x -> !any(f -> f(x), (ismissing)), df)
    #df = filter( [5]=> x -> isnumeric, df)
    #parse.(Date, df[:,1:4])
    #parse.(Date, string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH")
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
    df=df[:,Not(1:4)]
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
            scatter(x=df.date, y=df[:,i],
            name=st[i]),   row=i,     col=1);
    end
    relayout!(p,height=600*2,width=900*2,title_text="Series of "*basename(path))
    p
end


##################query################
function f_read(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    dfs = DataFrame[]
    for file in files
        file_path = joinpath(path, file)
        #if isfile(file_path) && endswith(file, ext)
        #if isfile(file_path) && contains(file, ext)
        if isfile(file_path) && occursin(ext, string(file))
            println(replace("reading  $file_path", "\\"  => "/")," now...")
            df = CSV.read(file_path,DataFrame,missingstring="-9999",comment="-",delim="\t")  |> @dropna() |> DataFrame
            #df = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), dd)
            df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
            df=df[:,Not(1:4)]
#            println(replace("reading on $file_path", "\\"  => "/")," done!")
            push!(dfs, df)
        elseif isdir(file_path)
            dfs_in_subdir = f_read(file_path, ext)
        end
    end
    return(dfs)
end
directory_path = "/mnt/d/Wasim/regio/out/v8"
ext = "v8.2017"
#ext = "albe*v8"
#ext = "albe*v8"
f_read(directory_path,ext)

path="/mnt/d/Wasim/regio/out/v8/qgkorcm.v8.2017"
path="/mnt/d/Wasim/regio/out/v8/wurzrcm.v8.2017"

ms=["-999","-9999","lin","log","LIN","LOG"]
df = CSV.read(path,DataFrame,missingstring=ms,
delim="\t",
comment="-",
#ignorerepeated=true,
silencewarnings=true,
ntasks=4,
downcast=true,
#types=Dict(1=>Int,2=>Int,3=>Int),
#types=Dict(1=>String,2=>String,3=>String),
#typemap=Dict(String=>Int64))
normalizenames=true,drop=(i, nm) -> i == 4)







#drop=(i, nm) -> i == 2)
#   •  drop: inverse of select; an AbstractVector of Integer, Symbol, String, or Bool, or a "drop" function of the form (i, name) ->
#drop::Bool; columns in the collection or for which the drop function returns true will ignored in the resulting CSV.File. Invalid
#values in drop are ignored.

#ArgumentError: `drop` keyword argument must be an `AbstractVector` of `Int`, `Symbol`, `String`, or `Bool`, 
#or a selector function of the form `(i, name) -> keep::Bool`

df = df[completecases(df), :]
df = dropmissing(df, :) #same
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
df=df[:,Not(1:4)]


df.date = Date.(string.(df.YY,df.MM,df.DD),"yyyymmdd");


#############mask
using Rasters, Plots, Dates
fn="/mnt/d/remo/cordex/eobs/rr_ens_spread_0.1deg_reg_v26.0e.nc"
shp="/mnt/d/Wasim/regio/rcm/ez-wgs.shp"

# Load and plot the file
#awap = read(Raster(fn, :rr; date=DateTime(2001,01,01))) #,"yyyy-mm-dd"
#awap = read(Raster(fn, :rr; date=Date(2001,01,01)))
awap = read(Raster(fn, :rr; date=DateTime(2001,12,1)))
awap = Raster(fn,key=:rr)
a = plot(awap; clims=(10, 45))
# Create a mask my resampling a worldclim file
wc = Raster(WorldClim{Climate}, :prec; month=1)
wc_mask = resample(wc; to=awap)
# Mask
awap_masked = mask(awap; with=wc_mask)
b = plot(awap_masked; clims=(10, 45))

using Shapefile
shp = "/mnt/d/Fernerkundungsdaten/glo30/eubasins.shp"
ezg = Shapefile.Handle(shp)
#EPSG(raster)
ezg|>plot
r = Raster("/mnt/d/Fernerkundungsdaten/glo30/sglo30_basin.tif",missingval=-32768)
#poly = ezg
mask_trim(raster, poly) = trim(mask(raster; with=poly); pad=10)
#ezg.shapes[end]
msk=mask_trim(r, ezg.shapes[5])
plot(msk)



using ArchGDAL; const ag = ArchGDAL
z = ag.readraster(fn)
dem = ArchGDAL.readraster("/mnt/d/Fernerkundungsdaten/glo30/sglo30_basin.tif")
#dem|>Plot.plot
g = ag.load()

ts=Raster(file,missingval=-9999)
x = ts[t=lyr]
print("saving countour plot to",outfile,"...")
#x = ts[t=1]
#x = x .> 0
x=trim(x)

#mean_height = map(x -> mean(skipmissing(replace_missing(mask(dem_model; to=x)))), 	vect.geom)

# rmask = ts .> 0
# rmask = x .> 0 
# #g = mask(x;rmask)
# rmask = x .<= 0
# g = x .* rmask
# #g = x .- rmask
# #g = x .> 1
# g = x .> rmask
# g |> Plots.plot

# w= replace_missing(x)
# w|>Plots.plot
#using Statistics
#d = zonal(mean, ts;of=:t)
#crs(x)
#mask(x,with=x>0)

#p = contourf(x; dpi=300, size=(800, 400))

file="/mnt/d/Wasim/regio/rcm/rcm.art_bfid"

# using GMT
# C = makecpt(0,5);
# imshow("forcing-moselle.nc?precip", proj=:guess, cmap=C, colorbar=true)

pt="/mnt/d/Wasim/regio/rcm/"
walkdir(pt)
s=file
m=match(r".*[.]",s)
outfile = string(m.match,"png")
ts=Raster(file,missingval=-9999)
lyr=1
x = ts[t=lyr]
#trim(x)|>plot
x=trim(x)
contourf(x;)
using Shapefile
shapefile_name = "/mnt/d/Wasim/regio/rcm/ezg.shp"
ezg = Shapefile.Handle(shapefile_name)
mask_trim(raster, poly) = trim(mask(raster; with=poly); pad=10)
ezg.shapes[end]
msk=mask_trim(ts, ezg.shapes[end])
plot(msk)

msk=mask_trim(r, ezg.shapes[1:end])
plot(msk)
contourf(msk)

#using NCDatasets
using Rasters
fn="/mnt/d/remo/cordex/eobs/rr_ens_spread_0.1deg_reg_v26.0e.nc"
shp="/mnt/d/Wasim/regio/rcm/ez-wgs.shp"
fn
#rra = Raster(fn,:rr)
using Dates
#awap = read(Raster(fn, :rr; time=Date(2001,12,1)))
w = Raster(fn, :rr;)

Pkg.resolve()

using PlotlyJS, CSV, DataFrames
md="/mnt/d/Wasim/ecad/metadata.txt"

df = CSV.read(md,DataFrame,ntasks=4,downcast=true)

dx = df|> dropmissing

# ms=["-999","-9999","lin","log","LIN","LOG"]
#     df = CSV.read(md,DataFrame,ntasks=4,downcast=true)
#     #missingstring="-9999", #also windows
#     missingstring=ms,
#     delim="\t",comment="-",
#     silencewarnings=false,
    
#     normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing


fn="/mnt/d/Wasim/ecad/sources.txt"
dd = CSV.read(fn,DataFrame)
dropmissing!(dd)

first(dd)
#names(dd) = parse.(String,(dd[1,:]))
#dd = dd[1:end,:]


using Query
const q = Query

dd = CSV.read(fn,DataFrame,skipto=25,header=24,stripwhitespace=true)
print(names(dd))
subset(groupby(dd, :STAID), :LAT => x -> maximum(x) < 50)


fn="/mnt/d/Wasim/ecad/QQ_SOUID210797.txt"
dd = CSV.read(fn,DataFrame,skipto=23,header=22,missingstring="-9999",stripwhitespace=true)|>dropmissing


##df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");
using Dates
dd.DATE = Date.(string.(dd.DATE),"yyyymmdd")

using Plots
dd[:,3:4]|>Plots.plot
Plots.bar(dd.DATE,dd.QQ;)
last(dd,3)


function f(file)
    result = String[]
    for line in eachline(file)
        if startswith(line, "LATITUDE") || isempty(line)
            continue # skip the header and any empty lines
        elseif startswith(line, "ELEVATION")
            break # stop reading completely
        end
        push!(result, line)
    end
    return result
end

cd("/mnt/d/Wasim/ecad/rr")
z=readdir()[2]
mm = f(z)
mm


ms="-999999"
#file = CSV.File(z; limit=20,missingstring=ms)
file = CSV.read(z,DataFrame; limit=10,missingstring=ms)


path=readdir()[2]
re="LON"
for line in eachline(path)
    z=filter(el->occursin(re, string(el)),line)
    v=CSV.read(z, DataFrame, header=0)
end


z=readdir()[5]
#m=(filter(line -> occursin(r"^[L]|^[ELEVATION]",line),readlines(open(z))))
m=(filter(line -> occursin(r"^L|^ELEV",line),readlines(open(z))))
#m=DataFrame(filter(line -> occursin(r"^L|^ELEV",line),readlines(open(z))))

df = DataFrame(A=Int[], B=String[])
#push!(df, Dict(:B => "daf", :A => 3))
#push!(df, Dict(:B => m, :A =>[1,2,3]))
push!(df, Dict(:B => m[1], :A =>1))

df = DataFrame(B=String[])
for i in 1:length(m)
    push!(df, Dict(:B => m[i]))
end
#eachrow(push!(df, Dict(:B => m)))

df = DataFrame(B=String[])
m=(filter(line -> occursin(r"^L|^ELEV",line),readlines(open(z))))
for i in 1:length(m)
    push!(df, Dict(:B => m[i]))
end
using Query
df |> @map({a=_.B * "  -> added string ", _.B}) |> DataFrame
#replace("Hello worldHello", r"Hello$"=>"")

df |> @map({before=_.B,after=replace(_.B,"+"=>"")}) |> DataFrame
df |> @map({before=_.B,after=replace(_.B,r"^[^+]*+|:"=>" ")}) |> DataFrame

#remove +only
df |> @map({before=_.B,after=replace(_.B,r"[([^+)]*"=>"")}) |> DataFrame

CSV.read(z,DataFrame; skipto=4, limit=10,missingstring=ms)

pt="/mnt/d/Wasim/ecad/rr/stations.txt"
CSV.read(pt,DataFrame;skipto=2000,limit=5,header=16)
df = CSV.read(pt,DataFrame;skipto=9900,limit=5,header=16,stripwhitespace=true)
propertynames(df)


df[:,:HGT]
df[:,4]

reg = r"^[^+]*+|:"
reg = r"[([^+)]*"
df |> @map({name=_.STATIONNAME,before=_.LAT,after=replace(_.LAT,reg=>"")}) |> DataFrame

##some other tuts
iris = CSV.read((joinpath(dirname(pathof(DataFrames)),
                                 "..", "docs", "src", "assets", "iris.csv")),DataFrame)

sort!(iris, [:Species, :PetalLength], rev=[true, false])
sort!(iris, [order(:Species, rev=true), :PetalLength])

@time df = CSV.read(pt,DataFrame;skipto=9900,limit=5000,header=16,stripwhitespace=true)
@time df = CSV.read(pt,DataFrame;skipto=9900,header=16,stripwhitespace=true)
@time df = CSV.read(pt,DataFrame;skipto=9900,header=16,stripwhitespace=true,ntasks=4)
@time df = CSV.read(pt,DataFrame;skipto=9900,header=16,stripwhitespace=true,ntasks=16)

propertynames(df)


using Plots; gr()
x = range(-3, 3, length=30)
Plots.surface(
  x, x, (x, y)->exp(-x^2 - y^2), c=:viridis, legend=:none,
  nx=50, ny=50, display_option=Plots.GR.OPTION_SHADED_MESH,  # <-- series[:extra_kwargs]
)

#set other backend
plotlyjs()


pt="/mnt/d/Wasim/ecad/ger-info.txt"	
fd = CSV.read(pt,DataFrame;skipto=3,header=2)

fn="umst_peca13074.dat"
dd = CSV.read(fn,DataFrame,skipto=23,header=22,missingstring="-9999",stripwhitespace=true)|>dropmissing
using Dates
dd.DATE = Date.(string.(dd.DATE),"yyyymmdd")

using Plots
dd[:,3:4]|>Plots.plot
Plots.bar(dd.DATE,dd.QQ;)
last(dd,3)

###so geht beides nicht aus geopandas ->julia... :( aber mit 
# using KML
# fp="/mnt/d/Wasim/ecad/rr/climexp230407.kml"
# file = KMLFile(fp)
#Plots.plot(file)
gd = GeoDataFrames.read(fp)
Plots.plot(gd.Name)

using Shapefile
path="/mnt/d/Wasim/ecad/rr/gdfcrds.shp"
table = Shapefile.Table(path)
# if you only want the geometries and not the metadata in the DBF file
geoms = Shapefile.shapes(table)

###geht
fp="/mnt/d/Wasim/ecad/rr/crds.json"
using GeoJSON, DataFrames
jsonbytes = read(fp);
fc = GeoJSON.read(jsonbytes)
using DataFrames
df = DataFrame(fc)


using ArchGDAL
x=ArchGDAL.read(fp)
d = DataFrames.DataFrame(ArchGDAL.getlayer(x, 0))
d
#https://github.com/JuliaGeo/Leaflet.jl
#using Leaflet, Blink, GADM
fnt="/mnt/d/Wasim/Goldbach/out/sglo_v45/ei_masked.v45.tif"
fnt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/precbr50.sum.nc"
ArchGDAL.imread(fnt)

spatialref = ArchGDAL.importEPSG(4326)
#Spatial Reference System: +proj=lcc +lat_0=45.3333333333333 + ... _defs
ArchGDAL.toPROJ4(spatialref)


using Rasters
t2=read(Raster(fnt))
using Plots
t2|>plot
t2[Band=5]|>plot

x = t2[Band=5]
y = x .* 0.1
y = x .+ 0.1
y|>contourf

fnt="/mnt/d/remo/cordex/wgs/utm/t2.nc"
fnt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/precbr50.sum.nc"

x = read(Raster(fnt,missingval=-9999))
x|>plot
y = x .+ 0.1
y|>contourf
x[:,1:50, 1]|>bar
x[20:125,5:50, 1]|>bar
x[50:125,1:50, 1]|>plot

cl=colormap("RdBu")
plot(x;c=cl)

#z = read(Raster(fnt,missingval=[0:10,-9999]))
#####so kann man values croppen! juhuuu
z = read(Raster(fnt,missingval=0))
plot(z[t=1];c=cl)

#set other backend
plotlyjs()
plot(z[t=1];c=cl)
#switch backend
gr()
x|>plot


fnt="/mnt/d/Wasim/Goldbach/revision/v5/albeggb.2021.nc"
z = read(Raster(fnt,missingval=0))
plot(z[t=1];c=cl)

pt="/mnt/d/Wasim/ecad/rr/dom_crds.tsv"
pt="/mnt/d/Wasim/ecad/rr/dom_crds.csv"
pt="/mnt/d/Wasim/ecad/rr/domcrds.csv"
pt="/mnt/d/Wasim/ecad/rr/domcrds.shp"
using GeoDataFrames,CSV,DataFrames
const gd = GeoDataFrames
#df = CSV.read(pt,DataFrame;comment="#",delim="\t")
df = CSV.read(pt,DataFrame;comment="#",select=1:4) #|>dropmissing
#df = gd.read(pt,options=["GEOM_POSSIBLE_NAMES=lon,lat"]) #,skiprows=1
#
###die layer auswahl << 0 >> oder name ist wichtig!
gdf = gd.read(pt, options=["GEOM_POSSIBLE_NAMES=points", "KEEP_GEOM_COLUMNS=NO"], 0)
propertynames(gdf)
gdf.geometry|>plot

#reprojection #https://juliapackages.com/p/geodataframes
using GeoFormatTypes; const GFT=GeoFormatTypes
df = gdf
rename!(df,:geometry => :geom)
df.geom = gd.reproject(df.geom, GFT.EPSG(4326), GFT.EPSG(25832))
df
propertynames(df)

##add columns to df
transform!(df, [:lon,:lat] => (*) => :crds) #will overwrite if exists
insertcols!(gdf, crds => gdf.lon + gdf.lat)

###get data in jl
fn="/mnt/d/Wasim/ecad/rr/rr_dom.txt"
ms="3.00e+33"
dd = CSV.read(fn,DataFrame,skipto=90000,
#header=1,
missingstring=ms,stripwhitespace=true,ntasks=4) #|>dropmissing


###newrastplot
fnt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/precbr50.sum.nc"
cl=colormap("Purples",logscale=true,1000)
z = read(Raster(fnt,missingval=0))
plot(z[t=1];c=cl)
plot(z[t=1];c=cgrad(:thermal, rev = true))
#see: https://docs.juliaplots.org/latest/generated/colorschemes/
plot(z[t=1];c=:oslo)

#
#import Pkg; Pkg.add("PythonPlot")
using Plots; pythonplot()
x = y = collect(range(-π, π; length = 100))
fn(x, y) = 3 * exp(-(3x^2 + y^2)/5) * (sin(x+2y))+0.1randn(1)[1]
surface(x, y, fn, c=:thermal, extra_kwargs=Dict(:subplot=>Dict("3d_colorbar_axis" => [0.9, 0.05, 0.05, 0.9])))


import Pkg; Pkg.add("UnicodePlots")
#https://docs.juliaplots.org/latest/backends/#[UnicodePlots](https://github.com/JuliaPlots/UnicodePlots.jl)
using Plots; unicodeplots()

extra_kwargs = Dict(:subplot=>(; border = :bold, blend = false))
p = plot(1:4, 1:4, c = :yellow; extra_kwargs)
plot!(p, 2:3, 2:3, c = :red)

using StatsPlots
using RDatasets 
school = RDatasets.dataset("mlmRev","Hsb82")
@df school density(:MAch, group = (:Sx, :Sector), legend = :topleft)
propertynames(school)


using DataFrames, CSV, Dates  
path="/mnt/d/Wasim/regio/out/r4/qgesrcm.v4.2017"
skip=1
dd = CSV.read(path, DataFrame, header=skip,missingstring=["-999","-9999", "--"],
comment="-")



@df dd density(:10, group = (:YY, :MM), legend = :topleft)

using CSV,DataFrames,Dates,StatsPlots; gr()
function ddense(path,skip::Int,start::Int,stop::Int)
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


###windows version:
function ddense(path,skip::Int,start::Int,stop::Int)
    df = CSV.read(path,DataFrame,skipto=skip,
    missingstring="-9999",delim="\t",comment="-",
    normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
    println(propertynames(df))
    @df df density(cols(start:stop), legend = :topleft)
end


#outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"

skip=1
skip = skip>=1 ? skip : 6
skip=0




function waread(filepath,skip::Int)
    ms=["-999","-9999","lin","log","LIN","LOG"]
    skip = skip<=1 ? 6 : 1
    df = CSV.read(filepath,DataFrame,skipto=skip,
    missingstring=ms,delim="\t",comment="-",silencewarnings=false,
    ntasks=4,downcast=true,normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
#    return(df)
end

#skip=6
df = CSV.read(filepath,DataFrame,skipto=skip,
    missingstring=ms,delim="\t",comment="-",silencewarnings=false,
    ntasks=4,downcast=true,normalizenames=true,drop=(i, nm) -> i == 4)

fnr="/mnt/d/Wasim/ecad/rr/tst.wa"
z = waread(fnr,5)
fnr="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/wasserkuppe.wind"
fnr="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/temp.2021"
ddense(fnr,8,1,3)


using RDatasets
iris = dataset("datasets", "iris")
@df iris andrewsplot(:Species, cols(1:4), legend = :topleft)

#AndrewsPlots are a way to visualize structure in high-dimensional data by depicting each row of an array or table as a line that varies with the values in columns. https://en.wikipedia.org/wiki/Andrews_plot
@df df andrewsplot(:tot_average, cols(1:4), legend = :topleft)
@df df andrewsplot(:date,cols(1:4),legend = :topleft)

using PkgCleanup
PkgCleanup.manifests()
PkgCleanup.artifacts()
using Pkg
using Dates
Pkg.gc(; collect_delay=Dates.Day(0))


ddense(v[end],4,2010,2011)

    # s = propertynames(i)[1:end-1]
    # @df i plot(:date,cols(s))
    # s = propertynames(y)[1:end-1]
    # @df y plot!(:date,cols(s))

    adf = vcat(dfs,cols = :date)
    adf = innerjoin(dfs,on = :date)
    
    
    adf = vcat(DataFrame(dfs[:RSM_77][:]))
    
    s = propertynames(dfs[1])[1:end-1]
    @df dfs[1] plot(:date,cols(s))
        for i in dfs[2:end];       
            s = propertynames(i)[1:end-1]
            @df i plot!(:date,cols(s))
        end
    
        
        for i in dfs[6:end];       
            s = propertynames(i)[1:end-1]
            @df i plot(:date,cols(s))
        end
    
    
v=ct("gw")
#z = v[Not(broadcast(x->occursin(r"nc$|xml|put",x),v))]
z = v[(broadcast(x->occursin(r"nc$",x),v))]
#z = v[(broadcast(x->occursin(r"^(?=[th]+$)(?=[nc]+$).*",x),v))]



# macro p_str(s) s end
# x = p"thrcm"
# y = p"nc$"
# z = Regex(join([x,y], "+"))
# l = v[(broadcast(x->occursin(z,x),v))]


# s="I'm searching for my funny words inside this text"
# occursin(".*(my)+.*(word)+|.*(my)+.*(word)+.*", s)
# occursin(r"/^(?=.*my)(?=.*word).*$/im", s)

#REGEX AND!!!
z = v[(broadcast(x->occursin(r".*(th)+.*(nc)+.*",x),v))]
a=[]
for i in z; r=Raster(i,missingval=0); push!(a,r);end
broadcast(x ->size(x),a)

#for i in a; plot(i);end

a[2]|>plot
a[4]|>plot

r=Raster(z[i],missingval=0)
plot(r)

x=a[4]+a[6]+a[8]
plot(x)
for i in z; r=Raster(i); push!(a,r);end

using Rasters, Plots, Dates, Shapefile, Downloads
using Rasters.LookupArrays
"/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/tar11"|>cd
shapefile_name="/mnt/d/Relief_DGMs/FABDEM/franken_fabdem_domain.shp"
shp = Shapefile.Handle(shapefile_name).shapes[1]

#rr = ct("raa01-sf_10000-1105261550.nc")
rr = Raster("raa01-sf_10000-1105261550.nc",crs=EPSG(25832), mappedcrs=EPSG(25832))
rc = crop(rr; to=shp)
rplot(rc)

lls = ct(".nc")
#/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/tar11:        22610.95 MB
lls = lls[1:end-1]
stk = Raster.RasterStack(filename=lls)

#typeof(lls[2])
stk = RasterStack(lls[8],layersfrom="Band1";crs=EPSG(25832), mappedcrs=EPSG(25832))
#stk = RasterStack(string.(lls[2:8]);crs=EPSG(25832), mappedcrs=EPSG(25832))
stk|>plot
#filepaths = ["file1.tif", "file2.tif", "file3.tif"]

r1 =Raster(lls[8],crs=EPSG(25832), mappedcrs=EPSG(25832))
r2 =Raster(lls[9],crs=EPSG(25832), mappedcrs=EPSG(25832))

#z=stack(r1,r2)
#z=[r1,r2]
z=vcat(r1,r2)
plot(z)


rc = crop(z; to=shp)
rc|>plot


lls[8:9]

# filelist=(
#  Band1="/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/tar11/raa01-sf_10000-1112312050.nc",
#  Band1="/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/tar11/raa01-sf_10000-1112312150.nc",
#  Band1="/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/tar11/raa01-sf_10000-1112312250.nc",
#  Band1="/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/tar11/raa01-sf_10000-1112312350.nc")
 
 stk = RasterStack(filelist;crs=EPSG(25832), mappedcrs=EPSG(25832))

 typeof(filelist)

 ol = Tuple(lls[8:9])
 
 #g(x) = string(x)[1] == '-' ? x : string("-", x)
 g(x) = string(x)[1] == '.' ? pwd() : string("-", x)
 g(".")
 fd()

 "/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2"|>cd
 z=Raster("fab.flk";missingval=-9999)
 z|>Plots.contourf  

 z|>

x=rebuild(z)
x=mean(z[Band(1)])


backends()
plotlyjs()
r=Raster("/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2/fab.slp")
# time: 2023-03-27 11:55:05 W. Europe Summer Time
# mode: julia
#z_data = Matrix{Float64}(r.data)
z_data = Matrix{Float64}(r.data[:,:,1])
layout = Layout(
    title="Slope of Fab",
    autosize=false,
    width=500,
    height=500,
    margin=attr(l=65, r=50, b=65, t=90)
)
PlotlyJS.plot(PlotlyJS.surface(z=z_data), layout)



r=Raster("/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2/fab.dhk")
file="/mnt/d/Wasim/regio/out/lowres/c8/qd__rcm_0300.sum.nc"
file="/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/c2/spin3/windfab.2002.nc"
r = read(Raster(file;crs=EPSG(25832),missingval=0))
#r = read(Raster(file,:qdrcm0300;crs=EPSG(25832),missingval=0))

name(r)
#z_data = Matrix{Float64}(r.data[:,:,1])
#typeof(r.data[:,:,1])
z_data = (r.data[:,:,1])
layout = Layout(
    title="Elevation",
    template="seaborn",
    autosize=false,
    scene_camera_eye=attr(x=1.87, y=0.88, z=-0.64),
    width=500, height=500,
    margin=attr(l=65, r=50, b=65, t=90),
    scene=attr(
            xaxis_nticks=20,
            zaxis_nticks=4,
            #camera_eye=attr(x=0, y=-1, z=0.5), #like zoomin
            aspectratio=attr(x=1, y=1, z=0.7)
        )
    )
PlotlyJS.plot(PlotlyJS.surface(
    z=z_data,
    contours_z=attr(
        show=true,
        usecolormap=true,
        highlightcolor="limegreen",
        project_z=true
    )
), layout)



flags = Dict(:tr => [1000,1000], :r => :near)
r2 = warp(r[t=1],flags)

using DimensionalData
DD=DimensionalData
A=r[t=1]
#dimres = map(d -> DD.basetypeof(d), dims(A, (Y,X)))
#A.dims

dims(A)
r2 = Rasters.resample(r[t=1],(1000,1000);)
#,:sum
(1000,1000)|>typeof
r[t=1]|>plot



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

file="/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/c2/spin3/precfab.sum.nc"
r = read(Raster(file;crs=EPSG(25832),missingval=0))

rasf = median_filter(r)
# Plot the original and filtered raster
plot(rasf)
plot(r)

r==rasf
describe(r)
describe(rasf)


using ArchGDAL
const AG = ArchGDAL
ds = AG.read(file)
# Holen Sie sich die erste Band als Raster
r = AG.getband(ds, 1)
# Definieren Sie den Ziel-CRS als WKT-String
target_crs = "PROJCS[\"WGS 84 / UTM zone 32N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",9],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32632\"]]"
# Definieren Sie die Ziel-Pixelgröße in Metern
target_pixelsize = (100.0, 100.0)
# Reprojizieren Sie das Raster in das Ziel-CRS und die Pixelgröße
#geht auch nicht!
r_reprojected = AG.reproject(r, target_crs, pixelsize=target_pixelsize)

AG.reproject|>methods


###AG -> Rasters
file="/mnt/d/Relief_DGMs/FABDEM/gb.tif"
xxv = AG.readraster(file)
my = Raster(xxv[:,:,1],(X,Y);mappedcrs=25832)
flags = Dict(:tr => [1000,1000], :r => :near)
plot(warp(my, flags)) #errors

my = Raster(file)
my.dims
plot(warp(my, flags))
#plot(my)
contourf(my)

##Rcall -> nctotiff
file="/mnt/d/Wasim/regio/out/lowres/c8/albercm.2003.tif"
file="/mnt/d/Wasim/regio/out/lowres/c8/albercm.2003.nc"
my = Raster(file)
plot(warp(my, flags))

my

using NetCDF
nc = NetCDF.open(file)
var = nc["albercm"]
data = var[:]
#NetCDF.close(file)
vr = Raster(var[:,:,1],(Y,X);missingval=-9999.0)
plot(vr)
plot(warp(vr, flags))


ln="/mnt/d/Wasim/regio/out/lowres/c8/qgesrcm.c8.2003"
dfp(ln)
df = readdf(ln)
denseplot(df)
vio(df)

pt="/mnt/d/Relief/DGM50_GDI/2022/vis/ufra/ufra_depth_20.tif"
r=Raster(pt)
r
#r[Band=:soil_class_id]
r[Band=4] |> plot
r

pt="/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2/bkid_brend.tif"
#rplot(pt,1)
r=Raster(pt)
z=vcat(r[Band(1)],r[Band(2)])
z|>plot
pt="/mnt/d/Relief/DGM50_GDI/2022/evi.tif"
r=Raster(pt)
"/mnt/d/Relief/DGM50_GDI/2022"|>cd
v=rglob(".nc")
pt="./vis/ufra/ufra_depth_30000.nc"
#r=Raster(pt,:sand;)
r=Raster(pt)
#shapefile_name="/mnt/d/Relief/DGM50_GDI/2022/reg_bez.shp"
shapefile_name="/mnt/d/Relief/DGM50_GDI/2022/reg_bez_utm.shp"
using Shapefile
shp = Shapefile.Handle(shapefile_name)
mask_trim(raster, poly) = trim(mask(raster; with=poly); pad=10)
shp.shapes[1:end]
shp
#reproject(shp.shapes[1],EPSG(25832))
#ERROR: ArgumentError: Cannot reproject from:  ...
#rl = Rasters.reproject(EPSG(4326),r)
msk=mask_trim(r, shp.shapes[2])
plot(msk)

#r|>describe

pt="/mnt/d/Relief/grid-ger/vg250-ew_12-31.utm32s.shape.ebenen/vg250-ew_ebenen_1231/"
cd(pt)
using Shapefile
shp = Shapefile.Handle("VG250_STA.shp")
gr()
shp.shapes[1]|>Plots.plot

Shapefile.write("ger.shp", shp.shapes[1])
v=getf("ger")
readlines(v[3])


using Shapefile, GeoJSON
# load a shapefile into GeoInterface compatible objects
#handle = Shapefile.Handle("path/to/a.shp")
# write them to a geojson file
GeoJSON.write("output.geojson", shp.shapes[1])


#Shapefile.write("tst.svg", shp.shapes[1])

using ArchGDAL, GeoFormatTypes
# load a shapefile as a dataset
dataset = ArchGDAL.read("VG250_STA.shp")
# get the first layer of the dataset
layer = ArchGDAL.getlayer(dataset, 0)
ArchGDAL.crs2transform
ArchGDAL.getcoord(layer)
# create a new CRS from EPSG code
#new_crs = ArchGDAL.SpatialRef(EPSG(4326))
new_crs = spatialref = ArchGDAL.importEPSG(4326)
# reproject the layer to the new CRS geht nicht.
new_layer = ArchGDAL.reproject(layer, ArchGDAL.toWKT(new_crs))


#######
fn="/mnt/d/Wasim/regio/rcm/rcm.art_bfid.filled.nc"
art = Raster(fn)

rplot(art,366)
describe(art)

@functionloc p()
@less p() #get out by pressing <<q>>
@assert p()
function oops(msg::AbstractString)
    println("\n*** Oops, $msg\n")
    # close database, close files, .....
    ccall(:jl_exit, Cvoid, (Int32,), 86) # Exit code 86 to the O/S
    ##or :  exit(86).
end

# Verwenden Sie die Funktion error, um eine Ausnahme auszulösen und die Ausführung zu unterbrechen. 
# Dies ist nützlich, wenn Sie den Fehler nicht erwartet haben oder keine alternative Berechnung durchführen können
# #@error p()
# Verwenden Sie die Makros @assert oder @error, um eine Bedingung zu überprüfen und eine Ausnahme auszulösen oder eine Fehlermeldung zu protokollieren. 

#climateplot on time with Rasters.jl
fn="/mnt/d/remo/cordex/eobs/qq_ens_spread_crop.nc"
cd("/mnt/d/remo/cordex/eobs/")
glob("qq")

fn=glob("qq")|>first
z=Raster(fn)
z[X=4,Y=4]|>plot


z[X=4,Y=4]|>typeof

fn="d:/remo/cordex/eobs/qq_ens_spread_0.1deg_reg_v26.0e.nc"

# Rasters.aggregate
# using Rasters, Statistics, Plots
# using Rasters: Center
# fn="/mnt/d/remo/cordex/eobs/qq_ens_spread_crop.nc"
# st = read(Raster(fn; date=DateTime(2000,12,15,00)))
# z = read(Raster(fn))
# #ag = Rasters.aggregate(Center(), z, (Y(20), X(20));)
# plot(ag)


##reproject rasters.
using Rasters
import ArchGDAL
r = Raster("/mnt/d/Wasim/regio/rcm200/v4/rcm.dhk") #,crs="epsg:25832")
flags = Dict(
"s_srs"=>"EPSG:25832",
"t_srs"=>"EPSG:4326",
"r"=>"cubicspline")
#-r resampling_method

rs = warp(r,flags)
#using Makie,CairoMakie
plot(rs)

##########reproject rasters pt2  
l="/mnt/d/Wasim/sinn/out/v3-f1/sb05sinn100.2016.nc"
m = towin(l)
r = Raster(m)
r = trim(r[t=1])
k = project(r)
k = rebuild(k, name=:sb05 )
#cmk.mkrheat(k)
Plots.plot(k)

##crop out
r_bounds = X(9.3 .. 10), Y(50 .. 50.5)
kn = k[r_bounds...]
Plots.plot(kn,title=name(kn))

#resample
r = Raster("d:/Wasim/regio/rcm200/v14/rcm.dhk")
flags = Dict(
"s_srs"=>"EPSG:25832",
"t_srs"=>"EPSG:4326",
"tr" => [0.01,0.01],
"r"=>"cubicspline")
#-r resampling_method
rs = warp(r,flags)
Plots.plot(rs)
rs = rs[r_bounds...]
rs = replace_missing(rs,0)
Plots.surface(rs,camera=(30,50),c=:thermal,
    title="Sinn")
Plots.surface(r,camera=(30,50),name(r))



#########
fn="D:/remo/qm/corgrids/proj/pre-cor3.nc"
shapefile_path="D:/Wasim/regio/rcm200/v14/catchment_v14.shp"
shp = gpd.read_file(shapefile_path)
shp = shp.to_crs(crs="EPSG:25832")
ds = xr.open_dataset(fn)
myvar = "pre"
ds[myvar].mean("time").plot()
shp.boundary.plot(ax=plt.gca(), color="red", 
    linewidth=1.75, linestyle="dotted")
ds[myvar].mean(["x","y"]).groupby("time.year").sum().plot()
ds.close()
ds = nothing
GC.gc()
vars()
###das ist auch sehr komisch.
"D:/Wasim/Tanalys/DEM/brend_fab/out/c9/c9-eobs4/"
"D:/Wasim/Tanalys/DEM/brend_fab/out/c10/s2/"
## es waren xc und yc float in der gridlist....
##nc scaling_factor checken!!!!
fn="D:/remo/cordex/eobs/v28/biascorr/sfcWind_cor.nc"
ds = xr.open_dataset(fn)
myvar = "sfcWind"
#Dimensions without coordinates: t, x, y
ds[myvar].mean("t").plot()
#ds.t as Array{Dates.DateTime,1} ?
@pyimport datetime
@pyimport numpy as np
pyslice = pybuiltin("slice")
pyint = pybuiltin("int")
#pyint(12354.654)
dt = [datetime.datetime(1990, 1, 1) + datetime.timedelta(days=pyint(i)) for i in ds.t]
ds.close()
fn="D:/remo/qm/corgrids/jlcor/pre-cor_extremes_raw.nc"
ds = xr.open_dataset(fn)
myvar = "pre"

shp = shp.to_crs(crs="EPSG:4326")
collect(ds.keys())
collect(ds.dims)
ds[myvar].mean("time").plot()
shp.boundary.plot(ax=plt.gca(), color="red", 
    linewidth=1.75, linestyle="dotted")
ds.close()


### jlcorgrids are wrong.
fn="D:/remo/qm/corgrids/proj/rh-cor.nc"
#ds = Raster(fn,key=:rh)
ds = xr.open_dataset(fn)
shp = shp.to_crs(crs="EPSG:25832")
ds["rh"].mean("time").plot()
shp.boundary.plot(ax=plt.gca(), color="red", 
    linewidth=1.75, linestyle="dotted")
ds.close()

fn="D:/remo/qm/corgrids/jlcor/rh-cor.nc"
ds = xr.open_dataset(fn)
shp = shp.to_crs(crs="EPSG:4326")
ds["rh"].mean("t").plot()
shp.boundary.plot(ax=plt.gca(), color="red", 
    linewidth=1.75, linestyle="dotted")


#aber rawgrids stimmen.
fn="D:/remo/qm/corgrids/jlcor/rh-cor_raw.nc"
ds = xr.open_dataset(fn)
ds["rh"].mean("time").plot()
shp.boundary.plot(ax=plt.gca(), color="red", 
    linewidth=1.75, linestyle="dotted")

fn=raw"D:\remo\qm\corgrids\tas\tas_cor_raw.nc"
ds = xr.open_dataset(fn)
myvar = collect(ds.keys())[end]
ds[myvar].mean("time").plot()
shp.boundary.plot(ax=plt.gca(), color="red", 
    linewidth=1.75, linestyle="dotted")

cdof(fn)
glob("sh")
cp("D:/remo/qm/corgrids/proj/mw_proj.sh","mw_tas.sh")
npplat()
glob("nc")

fn = "D:/remo/qm/corgrids/tas/tas_cor_raw.nc"
ds = xr.open_dataset(fn)
ds
using PyCall
pyproj = pyimport("pyproj")
# Rename the coordinates
ds = ds.rename(Dict("time" => "t", "longitude" => "x", "latitude" => "y"))
# Reorder the coordinates
ds = ds.transpose("t", "x", "y")
ds.to_netcdf("D:/remo/qm/corgrids/tas/tas_cor2.nc")
# Reproject to EPSG 25832
transformer = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:25832", always_xy=true)
x, y = transformer.transform(ds.coords["x"].values, ds.coords["y"].values)
ds.coords["x"].values = x
ds.coords["y"].values = y

#ds = ds.drop("lat").drop("lon")
#ds = ds.rename_dims(Dict("longitude"=>"x","latitude"=>"y"))
#DataArray.transpose(*dims, transpose_coords=True, missing_dims='raise')
#ds = ds.transpose("time","x","y")
#ds.to_netcdf("D:/remo/qm/corgrids/tas/tas_cor2.nc")
ds.close()
ds = nothing

@vv "ncrename"
#ncrename -v latitude,y -v longitude,x tas_cor_raw.nc tas_cor2.nc

##preci
fn = "D:/remo/qm/corgrids/pre/pre-cor_extremes.nc"
ds = xr.open_dataset(fn)
ds["pre"].mean("time").plot()
println(ds.coords)
ds = ds.rename(Dict("time" => "t", "longitude" => "x", "latitude" => "y"))
# Reorder the coordinates
ds = ds.transpose("t", "x", "y")
ds.to_netcdf("D:/remo/qm/corgrids/pre/pre-cor2.nc")
ds.close()
ds = nothing
GC.gc()
##humi
fn = "D:/remo/qm/corgrids/jlcor/rh-cor_raw.nc"
ds = xr.open_dataset(fn)
ds["rh"].mean("time").plot()
println(ds.coords)
ds = ds.rename(Dict("time" => "t", "lon" => "x", "lat" => "y"))
# Reorder the coordinates
ds = ds.transpose("t", "x", "y")
ds.to_netcdf("D:/remo/qm/corgrids/rh/rh-cor2.nc")
ds.close()
##wind
#fn = "D:/remo/qm/corgrids/wind/sfcWind-cor_raw.nc" #<-has errors!
fn = "D:/remo/qm/corgrids/wind/sfcWind-cor_extremes.nc"
ds = xr.open_dataset(fn)
ds["sfcWind"].mean("time").plot()
println(ds.coords)
ds = ds.rename(Dict("time" => "t", "longitude" => "x", "latitude" => "y"))
ds = ds.transpose("t", "x", "y")
ds.to_netcdf("D:/remo/qm/corgrids/wind/sfcWind-cor2.nc")
ds.close()
ds = nothing
##radiation
fn = "D:/remo/qm/corgrids/rsds/rsds-cor_raw.nc"
ds = xr.open_dataset(fn)
ds["rsds"].mean("time").plot()
println(ds.coords)
ds = ds.rename(Dict("time" => "t", "longitude" => "x", "latitude" => "y"))
ds = ds.transpose("t", "x", "y")
ds.to_netcdf("D:/remo/qm/corgrids/rsds/rsds-cor2.nc")

#error in date columns of air_humidity
#29.2.2016 12:0 h

x="D:/Wasim/Stationsdaten_GKD/meteo-gls/rad_2022/rad_clean.2022"
df=waread(x)
baryrsum(df)
stp(x)
x="D:/Wasim/regio/out/rc200/y12/f4/Unterweissenbrunn-qoutjl"
theplot(x)



xn="D:/Wasim/regio/s50/rcm.tpi"
@edit agheat(xn)
agheat(xn;cmap=:viridis
agheat(xn;cmap=:greys,anns=false)


"D:/Wasim/regio/s50/"|>cd
using ArchGDAL
# Open the dataset
dataset = ArchGDAL.readraster("rcm.dhk")
r = ArchGDAL.getband(dataset,1)
# Perform the TPI operation
tpi = ArchGDAL.gdaldem("TPI", r)
# Save the result to a new file
ArchGDAL.write(tpi, "rcm.tpi", driver="AAIGrid")


using ArchGDAL
# Define a function for the "TPI" operation
tpi_op() = gdal_operation("TPI")
# Convert the Raster to an AbstractDataset
ds = ArchGDAL.read("rcm.dhk")
# Call gdaldem with the AbstractDataset and the tpi_op function
ArchGDAL.gdaldem(ds, tpi_op())


#In the case of unsafe_gdaldem, it’s a lower-level function that interfaces directly with the GDAL library, which is written in C. The unsafe_ prefix indicates that this function doesn’t perform certain safety checks that would normally be done in Julia, so you need to be careful when using it.
# ArchGDAL.unsafe_gdaldem(
#     dataset::AbstractDataset,
#     processing::String,
#     options = String[]; # List of options (potentially including filename and open options). The accepted options are the ones of the gdaldem utility.
#     dest = "/vsimem/tmp",
#     colorfile)


using ArchGDAL
# Open the dataset
dataset = ArchGDAL.read("rcm.dhk")
# Perform the TPI operation
ArchGDAL.unsafe_gdaldem(dataset, "TPI", ["-of", "AAIGrid"], dest="rcm.tpi2")

using ArchGDAL

function gdaldem_tpi(input_file, output_file)
    dataset = ArchGDAL.read(input_file)
    ArchGDAL.unsafe_gdaldem(dataset, "TPI", ["-of", "AAIGrid"], dest=output_file)
end

@gl "dhm"
gdaldem_tpi("rcm.dhm", "rcm.tpi")
r=Raster("rcm.tpi2")
r=Raster("rcm.tpi")
Plots.plot(r;c=:lightrainbow,title="TPI",xlabel="",ylabel="")

"""
one of "hillshade", "slope", 
"aspect", "color-relief", "TRI", "TPI",
       "Roughness"
"""
function gdalop(input_file, output_file, operation::String="TPI")
    dataset = ArchGDAL.read(input_file)
    ArchGDAL.unsafe_gdaldem(dataset, operation, ["-of", "AAIGrid"], dest=output_file)
end


gdalop("rcm.dhm", "rcm.hsd","hillshade")
r=Raster("rcm.hsd")
Plots.plot(r;c=:greys,title="Hillshade",xlabel="",ylabel="")


#########extract from BIAScorr Rasters##############
"D:/Wasim/regio/rcm200/y7/x0/"|>cd
vgjl("surf")
agsurf("rcm.dhk")
r=Raster("rcm.dhk")
#ra=Rasters.aggregate(Rasters.Center(),r,(Y(125),X(125));)
ra=Rasters.aggregate(Rasters.Center(),r,2) #fact 2
ra=Rasters.aggregate(Rasters.Center(),r,20)
surface(ra;camera=(60, 60),xlab="",ylab="",zlab="Elevation")
ar = project(ra;misval=-9999.0f0)

import GeoDataFrames as gdf, GeoInterface
ge = gdf.read("D:/Wasim/regio/rcm/ezg_4326.json")

rglob("json","D:/Wasim/regio/rcm200/")

ge = gdf.read("D:/Wasim/regio/rcm200/y11\\catchment-y11-4326.geojson")
rcmoutline = reduce(gdf.union,ge.geometry)
#crds = GeoInterface.coordinates(rcmoutline)
#ar = rebuild(ar,missingval=-9999)
#rcmoutline = wa.reverse_coords(rcmoutline)
#arc = crop(ar;to=rcmoutline)
plot(arc)
arc = mask(ar;with=rcmoutline)
@vv " boolmask"
#zm = (arc .> float(0)) .& (ds .<= umask)
zm = (arc .> float(0))
zm = boolmask(zm;missingval=-9999)          
arc = Rasters.mask(arc; with=zm)

surface(arc;camera=(60, 60),xlab="",ylab="",
    zlab="[m]",title="DEM")


pww()
cd("D:/Wasim/regio/rcm200/y7/x1")
gdalop("rcm.dhm", "rcm.hsd","hillshade")
r=Raster("rcm.hsd")
Plots.plot(r;c=:greys,title="Hillshade",xlabel="",ylabel="")

r=Raster("rcm.dhk")
surface(r;camera=(60, 60),xlab="",ylab="",
zlab="dhk [m]",title="DEM 50m")  



