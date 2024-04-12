vgjl("wegfurth")
# zw Bischhofsheim, und Schönau a.d. Brend
#findindf(df, "23")        #ezg schweinhof

pt=raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022"
fl =  "wegfurth.wats"

x = joinpath(pt,fl)
using RCall
ct = raw"D:/Wasim/regio/rcm200/v8/catchment.shp"

@rput x
@rput ct
@rlibrary "terra"
R"
library(terra)
fin=sub(pattern='wats','geojson',x = x)
pt = terra::vect(fin)
ezg =terra::vect(ct)
plot(ezg)
plot(pt,cex=2,add=T)
"

gwf = waread(x)
gws = waread(r"gwst")
hgeo = waread(r"geo")

(r"geo")|>glob

findindf(df, "23")        #ezg schweinfhof
nd = mall([gwf, gws[!,Cols(23,"date")],hgeo])
#xd = innerjoin(gwf, gws[!,Cols(1,"date")], on = :date,makeunique=true)
xd = select(nd, 1:4)
names(xd)
describe(xd)
xd.gwm = xd[!,end] .- xd[!,1]
qplot(select(xd,[3,5])|>dropmissing)
select(xd,[2,3,5])|>dfl

pwd()
wawrite(xd,"sim-wegfurth.wa")

cd(raw"D:/Wasim/Tanalys/DEM/Grundwasser/gw2022")
v=["pfaffenhausen82a", "mühlfeldmu11", "stettens1", "heiligkreuzs8", "schönderlings6", "oberthulbab212", "rieneck164", "fellens4", "obersinn", "kothens5", "wegfurth", "unterelsbachs3"]
s = join(v,"|")   
z = glob(s)
z = filter(x->occursin(r"wats$",x),z)

om = mall(loadalldfs(z))
describe(om)
#odfr
wa.writedesc("gw-desc-merged.wats",om)
wa.writewa("gw-merged.wats",om)
dfp(om)

om = odfr("gw-merged.wats")
@vv "correlations[i]"
#@df om corrplot(cols(propertynames(om[!,Not(:date)])), grid = false)


cd("D:/Wasim/Tanalys/DEM/Grundwasser/gw2022/")
#ENV["PYTHON"]
vgjl("corr")
using PyCall
@pyimport pandas as pd


m = pd.read_csv("gw-merged.wats",delim_whitespace=true)
m.head()
m.columns
using PyPlot
#@pyimport matplotlib
#PyPlot.pyplot.switch_backend("TkAgg")  # or "Qt5Agg"
#m.iloc[:,1:5]
#m.drop([["YY", "MM", "DD", "HH"]],axis=1)
PyPlot.pygui(:qt5)
PyCall.pyexists("PyQt5")
pygui(true)         #set PyPlot to interactive mode
pd.plotting.scatter_matrix(m)
println(m.head())
println(m.columns)

s = m.iloc[:,4:]
pd.plotting.boxplot_frame(s) 


m.head()

column1 = convert(Array, m.fellens4)

function pydf_to_julia(py_df::PyObject)
    # Convert each column of the Python DataFrame to a Julia array
    col_names = py_df.columns  # Get the column names from the Python DataFrame
    col_arrays = [convert(Array, py_df[col]) for col in col_names]
    # Create a Julia DataFrame using the converted arrays and column names
    julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, col_arrays))    
    return julia_df
end


col_names = m.columns
col_arrays = [convert(Array, m[col]) for col in col_names]

function pydf(py_df::PyObject)
    # Convert each column of the Python DataFrame to a Julia array
    col_names = py_df.columns  # Get the column names from the Python DataFrame
    col_arrays = [convert(Array, py_df[col]) for col in col_names]
    jdf = DataFrame(col_arrays, :auto)
    #size(jdf)
    rename!(jdf, cn)
    return jdf
end



cn = convert(Array, col_names)
#col_arrays_t = permutedims(hcat(col_arrays...))
jdf = DataFrame(col_arrays, :auto)
names(jdf)
size(jdf)
nrow(jdf)
rename!(jdf, cn)
#dfp(jdf)



nd = pydf_to_julia(m)


py"""
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
        df.set_index('date', inplace=True)
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

pwd()

py"""
import pandas as pd
k=waread3("gw-merged.wats")
print(k.head())
"""

#@pyimport pandas as pd; pd.DataFrame(rand(5, 5))
py"""k.tail()"""

py"""k.reset_index(drop=False,inplace=True)"""
py"""k.tail()"""
#xo = py"""k.iloc[:,1:]"""
xo = py"""k"""
xo.head()
col_names = xo.columns
col_names = convert(Array, col_names)
xdf = wa.pydf(xo)
names(xdf)
dfp(xdf)

fn=raw"D:\Wasim\regio\out\rc200\x14\waba-input.wa"
# Python code to execute (concatenates two strings)
julia_str1 = "Hello"
julia_str2 = " World!"
py"""
result_str = $julia_str1 + $julia_str2
"""
# Access the result from Python in Julia
result_str = py"result_str"
println(result_str)
pyo = py"""waread3($fn).reset_index(drop=False)"""
pdf = wa.pydf(pyo)
names(pdf)
dfp(pdf)

fn = glob("wats")|>last
pyimport("sys").executable
using PyCall
# Import the pandas module
pd = pyimport("pandas")
# Access the version attribute
version = pd.__version__
println("Pandas version: ", version)


#pyo = py"""waread3($fn).reset_index(drop=False)"""



m = m.drop(["YY", "MM", "DD", "HH"], axis=1)
p = pd.plotting.scatter_matrix(m);
display(p)
show(p)

savefig(p,"scatter.png")
op()


@pywith pybuiltin("open")("file.txt","w") as f begin
    f.write("hello")
end
@pywith pybuiltin("open")("file.txt","r") as f begin
    print(f.read())
end 

#xd[!,1] .- hgeo[!,1]
dropmissing!(xd)
qplot(xd)
dpr(xd)
dfl(xd)
#@df xd[!,Not(:date)] heatmap
heatmap(
    select(xd,1),
    select(xd,2),
    )

first(xd,15)|>clipboard

x = 1:10
y = fill(NaN, 10, 100, 3)
for i = axes(y,3)
    y[:,:,i] = collect(1:2:20) .+ rand(10,100).*5 .* collect(1:2:20) .+ rand()*100
end

errorline(1:10, y[:,:,1], errorstyle=:ribbon, label="Ribbon")
errorline!(1:10, y[:,:,2], errorstyle=:stick, label="Stick", secondarycolor=:matched)
errorline!(1:10, y[:,:,3], errorstyle=:plume, label="Plume")

###########scattergeo#####################
using PlotlyJS
# Define the data for the scatter plot
data = scattergeo(
    lat = [40.712776, 51.507351, 35.689487],
    lon = [-74.005974, -0.127758, 139.691711],
    text = ["New York", "London", "Tokyo"],
    mode = "markers",
    marker_size = 10
)

# Create the layout for the plot
#PlotlyJS.relayout!(data,
lyout = PlotlyJS.Layout(
    title = "Scatter plot on a map",
    geo = attr(
        scope = "world",
        projection_type = "natural earth"
    )
)

PlotlyJS.plot(data,layout)

##########gwmap#######
raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022" |>cd
fn=glob("csv")
# dfs=[]
# for x in fn
#     dfgw = CSV.read(x,DataFrame,header=7,skipto=8,limit=10^5,decimal=',',select=1:2)
#     rename!(dfgw,1=>"date")
# end

import GeoDataFrames as GDF
import GeoFormatTypes as GFT

dfpt=[]
for x in fn
    dfhd = CSV.read(fn,DataFrame,header=2,limit=6,decimal=',',delim=";")
    crd = select(view(dfhd, 4:4, :),1:4)
    xd = DataFrame(geometry=AG.createpoint(parse(Float64, string(crd[1,2])),crd[1,end]), 
    name=basename(x))
    GDF.write(basename(x)*".geojson", xd; 
        layer_name=basename(x), 
        crs=GFT.EPSG(25832), driver="GeoJSON", 
        options=Dict("SPATIAL_INDEX"=>"YES"))
    push!(dfpt,xd)    
    println("$x done!")
end
latx()

dfpt

pts = reduce(hcat,dfpt)
pts = reduce((left, right) -> 
innerjoin(left, right, on = :geometry,makeunique=true), 
dfpt)


opt=[]
onames=[]
for x in fn
    dfhd = CSV.read(x,DataFrame,header=2,limit=6,decimal=',',delim=";")
    crd = select(view(dfhd, 4:4, :),1:4)
    k = parse(Float64, string(crd[1,2])),crd[1,end]
    if (k[1] < 50000 || k[2] < 50000)
        continue
    else
        pt = AG.createpoint(k)
        push!(opt,pt)    
        push!(onames,x)
        println("$x done!")
    end
end

opt
#opt = filter(x->(),opt)
ad = DataFrame(geometry=opt,name=map(x->replace(x,".csv"=>""),onames))
GDF.write("allpts.geojson", ad; 
        layer_name="Grundwasser", 
        crs=GFT.EPSG(25832), driver="GeoJSON", 
        options=Dict("SPATIAL_INDEX"=>"YES"))


#table = [(; geom=ad, name="test")]
ad.geometry
@vv ".geometry"

##YES reprojection

pts=[]
for geom in ad.geometry
    xc = AG.getx(geom,0)
    yc = AG.gety(geom,0)
    pt = AG.createpoint(xc,yc)
    pt = AG.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
    push!(pts,pt)
end

data = scattergeo(
    lat = map(x->AG.gety(x,0),pts),
    lon = map(x->AG.getx(x,0),pts),
    text = onames,
    mode = "markers",
    marker_size = 10
)

ndf = DataFrame(geometry=pts,name=map(x->replace(x,".csv"=>""),onames))
GDF.write("allpts_4326.geojson", ndf; 
        layer_name="Grundwasser", 
        crs=GFT.EPSG(4326), driver="GeoJSON", 
        options=Dict("SPATIAL_INDEX"=>"YES"))

# Create the layout for the plot
#PlotlyJS.relayout!(data,
lyout = PlotlyJS.Layout(
    title = "Groundwater Gauges",
    geo = attr(
        scope = "europe",
        projection_type = "natural earth"
    )
)

PlotlyJS.plot(data,lyout)


@which wa.#355#357

Main.wa.#355#357 = new_value


using ArchGDAL
using PlotlyJS
const AG = ArchGDAL
# using GeoJSON
# jsonbytes = read("allpts_4326.geojson") 
# fc = GeoJSON.read(jsonbytes)

#ad = GDF.read("allpts.geojson") 
ad = GDF.read("allpts_4326.geojson") 
#test without reprojection
data = scattergeo(
    lat = map(x->AG.gety(x,0),ad.geometry),
    lon = map(x->AG.getx(x,0),ad.geometry),
    text = ad.name,
    mode = "markers",
    marker_size = 8
)

lyout = PlotlyJS.Layout(
    title = "Groundwater Gauges",
    geo = attr(
        scope = "europe",
        projection_type = "natural earth"
    )
)

PlotlyJS.plot(data,lyout)        #no plot on 25832 :(



fn=glob("csv")
dfs=[]
for x in fn
    dfgw = CSV.read(x,DataFrame,header=7,skipto=8,limit=10^5,decimal=',',select=1:2)
    DataFrames.metadata!(dfgw, "filename", replace(x,".csv"=>""), style=:note);
    #rename!(dfgw,[1=>"date",2=>replace(x,".csv"=>"")])
    rename!(dfgw,1=>"date")
    push!(dfs,dfgw)
end

dfs = filter(x->nrow(x)>50,dfs)

x=dfs[5]

##write all out with map!
map(x->wawrite(x,DataFrames.metadata(x)|>only|>last|>z->(z*".wats")),dfs)
op()
wslpath()|>clipboard
#nds = loadalldfs(r"wats")


for x in dfs
    rename!(x,2=>DataFrames.metadata(x)|>only|>last)
end

filterplot("al",dfs)

dfp(dfs[1])
for x in dfs[30:end]
    display(dfp!(x))
end

adf = reduce((left, right) -> 
innerjoin(left, right, on = :date,makeunique=true), dfs)

wa.wawrite(adf,"gw-merged.wats")


## umlaut fixer
# dhg umlaut
# cp -v '/mnt/c/Users/Public/Documents/Python_Scripts/bashscripts/umlaut.pl' .
# find . -type f -name '*wats' | while read FILE ; do newfile="$(echo ${FILE} | ./umlaut.pl)" ; mv -v "${FILE}" "${newfile}";done
using PlotlyJS
plotlyjs()
#DF=readdf(r"*thul*wats")
#DF=wa.getdf("thul",dfs)
DF=last(dfs)

PlotlyJS.plot(DF, x=:date, z=Matrix(DF), 
    type="contour", 
    contours=attr(showlabels=true))


PlotlyJS.plot(DF, x=:date, 
    y=:weibersbrunns4,
    z=Matrix(DF), 
    type="contour", 
    contours=attr(showlabels=true))



##################
lk = raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022\gw-merged.wats"
#gwdat = waread(lk)
cd(dirname(lk))
fns=glob(r"wats$")
#dfs = loadalldfs(fns)
adf = reduce((left, right) -> 
    innerjoin(left, right, on = :date,makeunique=true), 
        map(x->readdf(x),fns[1:10]))

println(names(adf))
#wa.wawrite(adf,"gw-merged.wats")
dfp(adf)


###adjust headers...
cd(raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022")
#mkdir("bak")
fls=glob(r"wats$")
cpinto(fls,"bak")

df = waread(fls[1])

x=fls[end]

for x in fls
    df = waread(x)
    newname = split(DataFrames.metadata(df)|>only|>last|>basename, ".")[1]
    rename!(df,1=>newname)
    wawrite(df,x)
end

# function replace_string(filename::String)
#     basename = split(basename(filename), ".")[1]  # Extract the basename without the extension
#     content = read(filename, String)  # Read the content of the file as a string
#     new_content = replace(content, 
#     "Grundwasserstand [m ü. NN]" => basename)  # Replace the string with the basename
#     return new_content
# end

# file_list = ["file1.txt", "file2.txt", ..., "file46.txt"]  # Replace with your actual file list
# for file in file_list
#     modified_content = replace_string(file)
#     write(file, modified_content)
# end


#https://github.com/JuliaPy/PyCall.jl
__precompile__() # this module is safe to precompile
module MyModule
using PyCall

const scipy_opt = PyNULL()

function __init__()
    copy!(scipy_opt, pyimport_conda("scipy.optimize", "scipy"))
end
end

Main.MyModule.scipy_opt



import GeoDataFrames as GDF
import GeoFormatTypes as GFT
using CSV
using DataFrames
import ArchGDAL as AG
pt=raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022\Spielkasino_tmw.csv"

cd(dirname(pt))

dfhd = CSV.read(pt,DataFrame;header=2,limit=6,decimal=',',delim=";")
crd = select(view(dfhd, 4:4, :),1:4)
nam = split(basename(pt),"_")[1]
xd = DataFrame(geometry=AG.createpoint(parse(Float64, string(crd[1,2])),crd[1,end]), 
    name=nam)
    outname=nam*".geojson"
    GDF.write(outname,xd; 
        layer_name=nam, 
        #crs=GFT.EPSG(25832), 
        crs=GFT.EPSG(31468), 
        driver="GeoJSON", 
        geom_columns=[:geometry])
        
rm("Spielkasino.geojson")
#        options=Dict("SPATIAL_INDEX"=>"YES"))


#r"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022\Spielkasino.geojson"
pt=raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022\Fuchsstadt_Mu3.csv"
dfhd = CSV.read(pt,DataFrame;header=2,limit=6,decimal=',',delim="\t")
#crd = select(view(dfhd, 2:3, :),2)
dfhd= CSV.read(pt,DataFrame;header=3,limit=3,decimal=',',
    select=[2],delim="\t")

nam = split(basename(pt),"_")[1]
crd=view(dfhd,1:2,:)
xd = DataFrame(geometry=AG.createpoint(crd[1,1],crd[2,1]),name=nam)
#rm("Spielkasino.geojson")
outname=nam*".geojson"
rm(outname)
xd = hcat(xd[!,Not(Cols(r"geo"))],xd[:,Cols(r"geo")])
GDF.write(outname,xd; 
#        layer_name=nam, 
        #crs=GFT.EPSG(25832), 
        #crs=GFT.EPSG(31468), 
        driver="GeoJSON", 
        geom_columns=[:geometry])

setup()
plot(xd)

vgjl("wats")


fn=glob(r"Spielka|Fuchss")
fn = filter(x->!contains(x,"wats"),fn)
dfs=[]
for x in fn
    dfgw = CSV.read(x,DataFrame,header=7,skipto=8,limit=10^5,decimal=',',select=1:2)
    DataFrames.metadata!(dfgw, "filename", replace(x,".csv"=>""), style=:note);
    #rename!(dfgw,[1=>"date",2=>replace(x,".csv"=>"")])
    rename!(dfgw,1=>"date")
    push!(dfs,dfgw)
end
#dfs = filter(x->nrow(x)>50,dfs)
getnames(dfs)
first(dfs)
dfs[2]|>dfp
df=dfs[2]
rename!(df,2=>
DataFrames.metadata(df)|>only|>last|>z->(split(z,"_"))[1]
)
onam=DataFrames.metadata(df)|>only|>last|>z->(split(z,"_"))[1]
wawrite(df,onam*".wats")


fst=
ts = CSV.read(fn[1],DataFrame)
first(ts,16)|>println
dfgw = CSV.read(fn[1],DataFrame;
    header=2, delim=";",
    skipto=17,
    limit=10^5,decimal=',',select=1:2)
x=fn[1]
DataFrames.metadata!(dfgw, "filename", replace(x,".csv"=>""), style=:note);
rename!(dfgw,1=>"date")
select!(dfgw,[1,3])

onam=DataFrames.metadata(df)|>only|>last|>z->(split(z,"_"))[1]
rename!(dfgw,2=>onam)
dfgw[!, 2] = replace.(dfgw[!, 2], "," => ".")
dfgw[!, 2] = replace.(dfgw[!, 2], "---" => NaN)
dfgw[!, 2] = parse.(Float64,dfgw[!, 2] )
#dropmissing!(dfgw,2)
#parse.(Float64,dfgw[!, 2] )
#parse.(Date,dfgw[!, 1] )
dfgw[!, 1] = Date.(string.(dfgw[!, 1]),"dd.mm.yyyy")
#Date.(dfgw[!, 1])
dropmissing!(dfgw)
#dfgw[!, 2] = replace.(dfgw[!, 2], missing => NaN)
wawrite(dfgw,onam*".wats")
dfp(dfgw)
df = waread(r"Spielkasino.wats")
dfp(df)

op()
wslpath()|>clipboard
    
fn=raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\hiwi_saale\saale\Valentin\SWAT_Output_DWD\subout_ETmm"
xd = CSV.read(fn,DataFrame;)
et = waread(fn)
baryrsum(et)
f2 = raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\hiwi_saale\saale\Valentin\SWAT_Output_DWD\subout_PETmm"
pet = waread(f2)
baryrsum(pet)
td = innerjoin(pet,et,on=:date)
td = reorder_df(td)
td.Tdiff = td[!,1] .- td[!,2]
te(td)

ofl = replace(f2,"subout_PETmm"=>"Tdiff")
td = (reorder_df(td))
wawrite(td,ofl)

cd(raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\ANALYSIS")
fn=rglob(r"rad|rge")
#fn=rglob(r"qg")
#fn=dfonly(".")
ds = wa.qbb()
od = reduce((left, right) -> 
innerjoin(left, right, on = :basin,makeunique=true), 
        map(x->select(x,[1,4]),ds)               )
od = permutedims(od)
od = broadcast(x->parse.(Float64,x),od[Not(1),:])
filter(x->x>0,od[!,2])|>bar
yaxis!("COEFF.VAR.LIN")


raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\Goldbach_L\sglo"|>cd
kgegrepr()

# ds = wa.qbb()
# ds = filter(x->nrow(x)>0,ds)
# #filter(x->!ismissing(x)>0,ds)
# #filter(x->select(x,[1,4]),ds)
# sel = broadcast(x->select(x,[1,4]),ds)
# sel = filter(x->size(x,1)>0,sel)
# map(x->size(x,2),sel)
# map(x->isnumeric(x[!,2]),sel)

ds2 = qall(;recursive=true)
#select(sel[2],2)
sel = broadcast(x->select(x,[1,4]),ds2)
od = reduce((left, right) -> 
    innerjoin(left, right, on = :basin,makeunique=true), 
    sel             )
od = permutedims(od)
#od = broadcast(x->parse.(Float64,x),od[Not(1),:])
filter(x->x>-10,od[!,1])|>bar
yaxis!("COEFF.VAR.LIN")


raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\Goldbach_L\glo3"|>cd
kgegrepr()
cd("glo3_v4")
ls()
tree()
r"qg"|>glob
qgk()
glob("out")
sd = dfr(r"qges")
xc = ctl()|>first
xc = split(xc,"\"")|>first|>y->split(y," ")|>last
xc = replace(xc, "control"=>"D:/Wasim/Goldbach/control")
#infile = String(xc)
#wa.npp(infile)
infile=raw"D:\Wasim\Goldbach\control\glo3_v4_loc.ctl"
ofl = "route.txt"
routeg(infile, ofl)
sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/"
specfile=joinpath(sfpt,sfn)
sim = r"qges"|>glob|>first|>readdf
obs = readdf(specfile)
df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
df.n3 .= names(obs)[1]
rename!(df,1=>"sim",2=>"obs",3=>"name")
df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
df.name=map(x->replace(x,r"_>.*" => ""),df.name)
i=1
dm =innerjoin(
    sim[!,Cols("C"*string(i[1]),end)],
    obs,    
    on=:date)
onam = "Goldbach"*"-qoutjl"
DataFrames.metadata!(dm, "filename", onam, style=:note);
wawrite(dm,onam)
ftp(dm)
baryrsum(r"rg")
(r"qoutjl")|>dpr
(r"sb05")|>dfp
(r"qoutjl")|>dfp!
waba()
td=tdiff()
"D:/Wasim/streu/"|>cd
kgegrepr()
raw"out\coarse\v1"|>cd
kgegrep()
qgk()
glob("out")
glob("out")|>wa.second|>dpr
dm = glob("out")|>wa.second|>waread
ftp(dm)
savefig("qout-log.png")
"waba"|>glob
waba2()
#tsoilstrstack:default_cellsize = 226.898 ;
wslpath()|>clipboard
r = readras(r"sb1")
plot(r)
cdb()
nsegrepr()


x="D:/Wasim/Pegeldaten/fluesse-abfluss/schondra_graefendorf/24473055_beginn_bis_31.12.2021_tmw.csv"

hdr = CSV.read(x,DataFrame,limit=7,header=2)

ms=["-999","-9999","lin","log","LIN","LOG"]
md = CSV.read(x,DataFrame,header=11,
missingstring=ms,delim=";",decimal=',',    #skipto=11,
    limit=10^5,select=1:4)

plot(md.Datum, md.Mittelwert, ribbon=(md.Minimum,md.Maximum), fillalpha=0.2)

#Not(:Datum)

@df md plot(:Datum, cols(propertynames(md)[2:end]), 
        yaxis=:log)

fn="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/wasim_pest/app/input/sinn100.fzs"
fzplot2(fn)
fzplot(fn)
fn="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/wasim_pest/app/input/sinn100.fzs"
ezplot(dirname(fn))

"C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/wasim_pest/app/output/v0"|>cd
f=@gl "T_Lower_Boundary_Condition_000"
rpr(f)

ga = Rasters.read(Rasters.Raster(f))
ti=basename(f)
wa.rpr(ga)

x = ga
msk = float(4)
xs = (length(x.dims) > 2) ? x[t=Int(x.dims[3][end])+1]  : x
#plot(xs)
contourf(xs)
zm = (xs .> msk)
ti=Rasters.name(x)
M=Rasters.mask(xs; with=zm)|>trim;
contourf(M)
Plots.surface(
    M,
    legend=false,
    #legendtext="",
    grid=false,
    xlabel=" ",
    ylabel=" ",
    zlabel=" ",
    #zaxis = :ln,
    title=ti,
    camera = (20, 75, 22)        
    )




    using DataFrames
    using Plots
   
df = dfr("D:/Wasim/regio/out/rc200/x22/spin2/qgesrcm.x22.2017")



df = df[!,Cols(r"date|13")]

foreach(x->println(x),names(df))
s = Symbol.(filter(x -> !occursin(r"date|year", x), names(df)))
years = unique(Dates.year.(df[!, :date]))  # Use df[!, column] here

df
su = @rsubset df year.(:date)==years[1]
dfp(su)
typeof(years[1])
for yr in years[2:end]
    su = @rsubset df year.(:date)===yr
    #describe(su)|>println
    @df su Plots.plot!(:date,
        cols(s),legend = yr, title=false)
end
plot!()


function hydromon(df::DataFrame; leg = :outertopright)
    ti = try
        DataFrames.metadata(df) |> only |> last |> basename
    catch
        @warn "No basename in metadata!"
        raw""
    end
    s = Symbol.(filter(x -> !occursin(r"date|year|month"i, x), names(df)))
    years = unique(Dates.year.(df[!, :date]))  # Use df[!, column] here

    begin
        su = @rsubset df year.(:date)==years[1]
        su = monmean(su)
        @df su Plots.plot(:month,cols(s),label = years[1], 
            title=ti,legend=leg)
        for yr in years[2:end]
            su = @rsubset df year.(:date)===yr
            su = monmean(su)
            @df su Plots.plot!(:month,
                cols(s),label = yr)
        end
        month_abbr = ["Jan", "Feb", "Mär", "Apr", 
                    "Mai", "Jun", "Jul", "Aug", "Sep", 
                        "Okt", "Nov", "Dez"];
        #xticks!(0.5:11.5 , month_abbr)
        xticks!(1:12, month_abbr)
        plot!()
    end
end

hydromon(df)
dfm(df)
mbx(df)

df = dfr("D:/Wasim/regio/out/rc200/x22/spin2/qgesrcm.x22.2017")
hydromon(df)
hydromon(wa.selt(df,1))
hydromon(wa.selt(df,"C13"))

wa.selt(df,"C13")|>mbx

function hydro(df::DataFrame; leg = :outertopright, logy=false)
    ti = try
        DataFrames.metadata(df) |> only |> last |> basename
    catch
        @warn "No basename in metadata!"
        raw""
    end
    #df = wa.selt(df,4)
    s = Symbol.(filter(x -> !occursin(r"date|year|month"i, x), names(df)))
    years = unique(Dates.year.(df[!, :date]))
    #Scale of the axis. 
    #Choose from [:identity, :ln, :log2, :log10, :asinh, :sqrt].
    if logy 
        ylog = :log
    else
        ylog = :identity
    end

    begin
        su = @rsubset df year.(:date)==years[1]
        #unique(monthname.(su.date))
        tm_ticks = round.(su.date, Month(1)) |> unique

        @df su Plots.plot(:date,cols(s),
            yaxis = ylog,
            xticks=(tm_ticks, 
            #Dates.format.(tm_ticks, "uu/yyyy")), 
            Dates.format.(tm_ticks, "uu")), 
            xrot=45, xminorticks=true, 
            xlim=extrema(su.date),
            label = (years[1]),  #first
            title=ti,
            formatter=:plain,
            legend = leg
            );


        for yr in years[2:end]
            su = @rsubset df year.(:date)===yr
            #tm_ticks = round.(su.date, Month(1)) |> unique; 
            Plots.plot!(
                vec(Matrix(su[!,Not(:date)])))
            

            # @df su Plots.plot!(
            #     #:date,
            #     cols(s),
            #     yaxis=ylog,
            #     xticks=(tm_ticks, 
            #     #Dates.format.(tm_ticks, "uu/yyyy")), 
            #     Dates.format.(tm_ticks, "uu")), 
            #     xrot=45, xminorticks=true, 
            #     xlim=extrema(su.date),
            #     label = (yr), #first
            #     formatter=:plain)
        end
        plot!()
    end
end

df = dfr("D:/Wasim/regio/out/rc200/x22/spin2/qgesrcm.x22.2017")

hydro(wa.selt(df,2);logy=true)

hydro(df;logy=true)

wa.selt(df,4)|>dfl

dates = Date(2021, 1, 1):Day(1):Date(2023, 3, 28)
n = length(dates)
df = DataFrame(
  date=dates,
  col1=rand(n),
  col2=rand(n)
)
tm_ticks = round.(dates, Month(3)) |> unique; #quarterly ticks!
plot(dates, rand(n), 
    xticks=(tm_ticks, 
        Dates.format.(tm_ticks, "uu/yyyy")), 
        xrot=45, xminorticks=true, 
        xlim=extrema(dates))








        df = wa.selt(df,7)
        s = Symbol.(filter(x -> !occursin(r"date|year", x), names(df)))
        years = unique(Dates.year.(df[!, :date]))  # Use df[!, column] here
        ylog=:identity
        ti=""


        using DataFrames
        using Plots
        using Dates
        
function hydro(df::DataFrame; leg = :outertopright, logy = false)
            ti = try
                DataFrames.metadata(df).metadata[end].basename
            catch
                @warn "No basename in metadata!"
                raw""
            end
        
            s = filter(x -> !occursin(r"(?i)date|year|month", string(x)), names(df))
            years = unique(year.(df.date))
            years_str = string.(years)
        
            ylog = logy ? :log : :identity
            
            mn = [ monthabbr(x) for x in unique(month.(df.date)) ]
        
            hp1 = plot(#xticks = mn,
                       xrot = 45,
                       xminorticks = true,
                       yaxis = ylog,
                       #xlim = extrema(df.date),
                       title = ti,
                       formatter = :plain,
                       legend = leg)
        
            for yr in years
                su = filter(row -> year(row.date) == yr, df)
                hp1 = plot!(vec(Matrix(select(su, Not(:date)))),
                            label = yr,
                            formatter = :plain)
            end

            #lng = size(df,1) ./ length(years)
            #lng = size(df,1) ./ 12
            #monthday.(df.date)|>unique #das müsste besser gehen...
            st = 15:31:size(df,1)
            xticks!(st, mn)
            
        
            return hp1
end
df
hydro(df)

hydro(df;logy=true)



###if PGFPlotsX really needed, Latex has to be installed
##sudo apt install texlive #(280mb)

df=dfr(r"kothen+.*wats")
df=dfr(r"heili+.*wats")
df=dfr(r"fellen+.*wats")

dropmissing!(df)
median(tovec(df,1))
median(tovec(df,1))
mean(tovec(df,1))

df=dfr(r"wegf+.*wats")
plotly()
wa.hydro(df)
wa.hydromon(df)


nd = @rsubset df year.(:date)<=1989
dfp(nd)
wa.hydro(nd)

# Valid Operations
plotattr(:Plot)
plotattr(:Series)
plotattr(:Subplot)
plotattr(:Axis)
plotattr("size")

y = tovec(df,1)
Plots.scatter(y, thickness_scaling = 2)  # increases fontsizes and linewidth by factor 2
# good for presentations and posters
# If backend does not support this, use the function `scalefontsizes(2)` that scales
# the default fontsizes.
gr()
Plots.scatter(y, 
marker = (:hexagon, .05, 0.06, 
:green, stroke(.3, 0.2, :black, :dot)))

vgr("art-bfid")


v=["pfaffenhausen82a",  "oberthulbab212", "rieneck164", "fellens4", "obersinn", "kothens5", "wegfurth", "unterelsbachs3"]
s = join(v,"|")
z = glob(s)
z = filter(x->occursin(r"wats$",x),z)
om = mall(loadalldfs(z))
names(om)
wa.hydro(om;col=3)
#kothen:
# Grundwasserstand [m ü. NN]: 403,49
# Grundwasserstand unter Gelände [m]: 10,27
# Geländehöhe [m ü. NN]: 413,76

ko=selt(om,3)
ko.gw .= ko.kothens5 .- 413.76
wa.hydro(ko;col=3)


############dez 23#################################
@wasim
lk = raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022\gw-merged.wats"
cd(dirname(lk))
gwdat = waread(lk)
println(names(gwdat))
#makie
cmk.tsbar(gwdat)
my = cmk.yrmean(gwdat)
my.year
m = Matrix(my)
m[:,1] .= Int.(m[:,1])
f = Figure()
ax = Axis(f[1,1])
heatmap!(ax,m)
f

using CairoMakie

fig = Figure()
fig[1:2, 1:3] = [Axis(fig) for _ in 1:6]
supertitle = Label(fig[0, :], "Six plots", fontsize = 30)
sideinfo = Label(fig[1:2, 0], "This text is vertical", rotation = pi/2)
fig


############März 2024#################################
using Pkg
Pkg.gc(; collect_delay=Dates.Day(0))  #160mb

fn=raw"D:\Wasim\Pegeldaten\jossa.kml"
js = gdf.read(fn)
js.geometry

@time setup()       #24sec
cd(raw"E:\els")
#elsb = RasterStack(glob("tif"))
#elsb|>first|>plot
#mos = mosaic(first, elsb...)
#gdalbuildvrt oberelsbach.vrt *.tif
mos = Raster("oberelsbach.vrt")
res = resample(mos;method=:bilinear,res=100)
#res = project(res;src=EPSG(25832), dst=EPSG(25832))
contourf(res,levels=10,title="",xlabel="",ylabel="")
crs(res)
write("oberelsbach-100m.tif",res)
import GeoDataFrames as gdf
pt = "D:/Wasim/regio/rcm200/v13/cmtv13.shp"
pt = "D:/Wasim/regio/rcm200/v14/cmtv14.shp"
gp = gdf.read(pt)
ex = reduce(gdf.union,gp.geometry)
res = resample(mos;method=:bilinear,res=25)
plot(res)
plot!(gp.geometry,fillcolor=:transparent)
contourf(res)
surf(res)
# result_string = join(sort(ar)|>vec, "\n")
# cb(result_string)
# size(ar)
@edit fzplot(;dirpath="D:/Wasim/regio/rcm200/v14/")
rst.fzplot(;dirpath="D:/Wasim/regio/rcm200/v14/")
title!("FZS16")
cd("D:/Wasim/regio/rcm200/v14/")
savefig("fzs16.png")
rst.fzplot()
title!("FZS64")
savefig("fzs64.png")

rst.fzplot2("D:/Wasim/regio/rcm200/v14/rcm.fzs")


#############24.03.24#################################
cd(raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022")
fn=glob("csv")[1]
readlines(fn)[1:5]
@doc readdlm(fn,";")[1:5]

fn="D:/Wasim/regio/out/vt/vt3-f1/gwthrcm_1500.mit.nc"
@edit agheat(fn,roundto=4,step=50,cmap=:greys)
plotly()
rst.agheat(fn,roundto=4,step=50,cmap=:greys)



fn=raw"C:\Users\Public\Documents\web_tab.csv"
df = CSV.read(fn,DataFrame;            
    stripwhitespace=true,
    decimal=',',  
    select = 1:6,          
    types=Dict(
        1=>String,
        2=>String,
        3=>String,
        4=>String,
        #4=>Date,
        5=>Float64,
        6=>Float64
        ))
dropmissing!(df,4)
df[!, 4] = Date.(string.(df[!, 4]),"dd.mm.yyyy")
nn = replace.(names(df), "\uad" => "ss", "."=>"")
rename!(df,nn)
using DataFramesMeta
@rsubset df :2 != "BA"
cd(raw"D:\Wasim\Tanalys\DEM\Grundwasser\gw2022")
sf("json")
apt = gdf.read("allpts_4326.geojson")
df.name = lowercase.(df.Messstelle)
op()
sf("meta")
innerjoin(df,apt,on=:name)

gst = dfr("gwstrcm_saale_nn.wa")
names(gst)
df.name 

dfp("schoenderlings6.wats")

vgr("gw_idw")
vgjl("gw_idw")

# Get the names from both dataframes
#df_names = map(x -> length(x) >= 5 ? x[1:5] : x, df.name)
#df_names = map(x -> length(x) >= 3 ? x[1:nextind(x, 3)-1] : x, df.name)
df_names = first.(split.(df.name," "))
gst_names = names(gst)
# Initialize an empty array to store the matches
matches = []
# Loop over each name in df_names
for df_name in df_names
    # Loop over each name in gst_names
    for gst_name in gst_names
        # Check if df_name is a substring of gst_name
        if occursin(df_name, gst_name)
            # If it is, add it to the matches array
            push!(matches, df_name)
        end
    end
end
# Print the matches
println(matches)
size(matches)
df.nm = df_names
df2 = filter(row -> row.nm in matches, df)
#df2 = @rsubset gst names(gst) in matches
#df2 = select(gst,Cols(matches))

#Anti match
filter(!row -> row.nm in matches, df)

names(gst)
df2.Messstelle
filter(w -> occursin(r"fel",w.Messstelle), df)
#select!(df2,Not())
sort!(df2,ncol(df2)-3)
select(df2,Not(1:4))
names(gst)
#filter(w -> occursin(r"rie",w.Messstelle), df)
writedf(df2,"GW_metadata_saale.txt")
#npplat()

#df2.thk = select(df2,ncol(df2)-3) - select(df2,ncol(df2)-2)
df2.thk = df2[:, ncol(df2)-3] .- df2[:, ncol(df2)-2]

select(df2,Cols(r"Mess|leiter|Gelände"))

import GeoDataFrames as gdf, GeoInterface
#ge = gdf.read("D:/Wasim/regio/rcm/ezg_4326.json")
ge = gdf.read("D:/Wasim/regio/rcm200/y11\\catchment-y11-4326.geojson")
rcmoutline = reduce(gdf.union,ge.geometry)
#sx = gdf.within(rcmoutline,apt.geometry)
within_points = [ArchGDAL.within(point, rcmoutline) for point in apt.geometry]
sx = apt[within_points, :]
df3 = select(df2,Cols(r"Mess|leiter|nm"))
#rename!(df3,3=>"name")
#innerjoin(df3,sx,on=:name)
sx.nm = replace.(sx.name,
    "s1"=>"",
    "s3"=>"",
    "s4"=>"",
    "s8"=>"", 
    "mu11"=>"", 
    "furth"=>"furt", 
    "b212"=>"", r"[0-9]" => "")
df4 = innerjoin(df3,sx,on=:nm)
df4.lon = ArchGDAL.getx.(df4.geometry,0)
df4.lat = ArchGDAL.gety.(df4.geometry,0)
ss = select(df4,Not(3:5))
writedf(ss,"GW_Buntsandstein_saale.txt")

@doc wa.vef
