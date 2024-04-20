cd("/mnt/c/Users/chs72fw/.julia/dev/WaSiM/")
using Pkg
Pkg.activate(".")
Pkg.status()
Pkg.test()

force_recompile(package_name::String) = Base.compilecache(Base.identify_package(package_name))
force_recompile("NCDatasets")
Pkg.test()

#using Pkg
#Pkg.add(path="https://github.com/consumere/WaSiM.jl")

using WaSiM
using Pkg
Pkg.add.(["DataFrames", "CSV", "Statistics", "Dates", "StatsPlots", "Distributions", "DataFramesMeta", "DelimitedFiles", "Grep", "Printf", "PrettyTables", "Rasters", "NCDatasets", "ArchGDAL", "GeoInterface", "GeoDataFrames", "Shapefile", "InteractiveUtils", "Plots", "SHA"]);
v = Symbol.(["DataFrames", "CSV", "Statistics", "Dates", "StatsPlots", "Distributions", "DataFramesMeta", "DelimitedFiles", "Grep", "Printf", "PrettyTables", "Rasters", "NCDatasets", "ArchGDAL", "GeoInterface", "GeoDataFrames", "Shapefile", "InteractiveUtils", "Plots", "SHA"]);
for i in v
    @eval using $i
end
import WaSiM
function toMain()
    fnames = names(Main.WaSiM, all=true)
    for submodule in fnames
        @eval import Main.WaSiM.$submodule
    end
end
toMain()



# ##das müsste dann in dem Modul WaSiM stehen
# include("win/smallfuncs.jl")
# @cmk
#import Pkg; Pkg.add("CairoMakie")

import WaSiM
import WaSiM.du as d
d()
WaSiM.ls()
WaSiM.tree()


cd(src)
import WaSiM.rglob as g
#g("func-w")
g("raster")
g("py")
g("rca")
g("plotly")

# copy dev files to WaSiM.jl
src = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
dst = "/mnt/c/Users/chs72fw/.julia/dev/WaSiM/src"
files = ["./cairomakie.jl","./func-win.jl","./rasterfuncs.jl",
    "./pyjl.jl","./rcall.jl","./RCall_gof.jl","./wajs.jl"]
#"./win/smallfuncs.jl"
function copy_files_with_warnings(src::AbstractString, dst::AbstractString, files::Vector{<:AbstractString})
    for file in files
        #dst_file = joinpath(dst, basename(file))  # Remove leading "./" from file path
        
        src_file = joinpath(src,replace(file,"./"=>""))  
        dst_file = joinpath(dst, replace(file,"./"=>""))  # Remove leading "./" from file path
        
        if isfile(dst_file)
            printstyled("Warning: File $dst_file already exists and will be overwritten.\n",color=:red)
        end
        cp(src_file, dst_file; force=true)
    end
end
copy_files_with_warnings(src, dst, files)


#julia --threads auto -q --startup-file=no --project="." 
using Pkg; Pkg.add("JET")
__precompile__(false)
using JET
report_package("WaSiM")

cd("/mnt/c/Users/chs72fw/.julia/dev/WaSiM/")
using Pkg
Pkg.activate(".")
Pkg.status()
__precompile__(false)

#julia --threads auto -q --startup-file=no --project="."
using JET
report_package("WaSiM")

import WaSiM as wa
wa.lat()

using Pkg
Pkg.activate("/mnt/c/Users/chs72fw/.julia/dev/WaSiM")
Pkg.status()
import WaSiM as wa
ds=wa.findlog()

#add Leaflet, Blink
using Leaflet, Blink, Shapefile
#r = wa.readras("/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/c8/loc/sb1_fab_1000.mit.nc")
#layers = Leaflet.Layer.(r; color=:orange); 
#r = wa.agread("/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/c8/loc/sb1_fab_1000.mit.nc")
poly = Shapefile.Handle("D:/Wasim/regio/rcm200/v12/catchment.shp"|>towsl)
layers = Leaflet.Layer.(poly.shapes;); 
provider = Leaflet.Provider("OSM")
m = Leaflet.Map(; layers, provider, zoom=3, height=1000, center=[10.0, 50.0]);
w = Blink.Window(; body=m)

#Blink.AtomShell.electron()
#sudo apt-get install libasound2

z=Leaflet.Map(; layers, provider);
Blink.Window(; body=z) #nope.


pwd()
using Rasters, ArchGDAL, GeoDataFrames, DataFrames, Plots
g = GeoDataFrames.read("/mnt/d/Wasim/regio/rcm200/v6/catchment.shp")
import NCDatasets
r = Raster("/mnt/d/remo/qm/corgrids/jlcor/pre-REcor.nc";key=:pre)

out = ArchGDAL.reproject(g.geometry,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
r.metadata
#z = r[At("2012-01-01T12:00:00")]
z = r[Ti=365*10]

z = Rasters.mask(z;with=out)|>trim
plot(z)

using Dates
Pkg.add("CFTime")
using CFTime
z2 = r[Ti=Near(CFTime.DateTimeNoLeap(2016))]
#Plots.contourf(z2)
Rasters.mask(z2;with=out)|>trim|>plot

#df = r[X=Near(10),Y=Near(50)]
#df = r[X=Near(10.08),Y=Near(50.23)]
#ti = DateTime.(ti)
#dt = DateTimeNoLeap("21001231",dateformat"yyyymmdd");
# dt = parse(DateTimeNoLeap,"21001231",dateformat"yyyymmdd")
ti = Rasters.lookup(r,3)
dt=Date.(Dates.year.(ti), Dates.month.(ti), Dates.day.(ti))
#x,y = map(x->round(x ./2;digits=0),size(r)[1:2])
#size(r)
#r[X=Near(10),Y=Near(50)]
#df = DataFrame(r[X=Int(x),Y=Int(y)]',:auto)|>permutedims
#df = DataFrame(r[X=Near(10),Y=Near(50)]',:auto)|>permutedims

df = DataFrame(r[X=Near(10.08),Y=Near(50.23)]',:auto)|>permutedims
df = hcat(df,dt, makeunique=true)
rename!(df,1=>name(r),2=>:date)
r.metadata
DataFrames.metadata!(df, "filename", "pre-REcor.nc", style=:note);
df.pre = convert(Vector{Float64}, df.pre)
#using StatsPlots
import Pkg; 
Pkg.add("CairoMakie")
Pkg.add("Makie")
include("cairomakie.jl")
cmk.tsp(df)
cmk.tsbar(df)
cmk.cloudplot2(df)

using JET
__precompile__(false)
report_package("WaSiM")

cd("/mnt/d/remo/qm/corgrids/jlcor/")
describe(df)
ssup()
@wasim
#wa.wawrite(df,"pre_50_10.wa")
wa.wawrite(df,"pre_bad-kiss.wa")
#read obs.
sts = wa.stp("/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/pre_ce.txt")
sts = filter(x->occursin(r"Kiss"i,string(x.name)),sts)
ArchGDAL.reproject(sts.geometry,
    EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
###--> 10.08 50.23

ArchGDAL.reproject(sts.geometry,
    EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
    #snll= DataFrame(snll,:auto)

obs = wa.dfr("/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt")
ob = wa.selt(obs, r"Kiss")
pdf = wa.mall(df,ob)
wa.theme(:dao)
wa.hyeval(pdf;freq="Y",fun=sum,yscale=:identity)
Plots.xaxis!(rotation=45)
Plots.savefig("jlcor-pre_bad-kiss.svg")
wa.hyeval(pdf;yscale=:identity,freq="Q",fun=mean)
Plots.savefig("jlcor-bad-kiss-seasonal.svg")
wa.hyeval(pdf;yscale=:identity,freq="week",fun=mean)
Plots.savefig("jlcor-bad-kiss-week.svg")


using DataFrames, CSV, Statistics, Dates 
using GLM, StatsPlots, Distributions
using DataFramesMeta
using Grep, Printf
using PrettyTables
using Rasters
import DataFrames: combine, groupby, transform
import DelimitedFiles: readdlm
import NCDatasets
import ArchGDAL
import GeoInterface #for reverse_coords ->moved to rst
import GeoDataFrames
import Shapefile #for jlzonal
using SHA
import JSON #jsread
