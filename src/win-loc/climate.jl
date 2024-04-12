
#jl --threads auto -q
using Base.Threads
nthreads()
#import Pkg; Pkg.Registry.update(); Pkg.activate(); Pkg.instantiate();
#ENV["PYTHON"]="C:\\Users\\chs72fw\\Miniconda3\\python.exe"
using Conda
#Conda.add("cmocean") #<-musste noch dazu
using ClimateTools
#https://juliaclimate.github.io/ClimateTools.jl/stable/datasets/

#C = load(filename::String, vari::String; poly::Array, data_units::String, start_date::Tuple, end_date::Tuple, dimension::Bool=true)
fn="/mnt/d/remo/cordex/eobs/hu_ens_mean_crop.nc"
#  Temporal subsetting can be done by providing start_date and end-date Tuples of length 1 (year), length 3 (year, month, day) or 6 (hour, minute, second).
#C = load(fn,"hu";data_units="%",start_date=Tuple("2014"),end_date=Tuple("2015"))
#parsetuple(s::AbstractString) = Tuple(parse.(Float64, split(s, ',')))
#parsetuple("234,456,654")

fn="/mnt/d/remo/cordex/REMO11_2011_2015/prec_wgs_day_20110101-20151231.nc"
C = load(fn,"pr";)

#C = load(fn,"hu";data_units="%",start_date=Date(2014, 12, 31),end_date=Date(2015, 1, 31))
#parsetuple(p::String) = Tuple(Float64(x) for x in p)
#ClimateTools.regrid
am = annualmean(C)

# using Pkg; Pkg.add("PyCall")
# ENV["PYTHON"]=""
# Pkg.build("PyCall")
#add ClimatePlots
using ClimateTools
using ClimatePlots
ENV["MPLBACKEND"] = "TkAgg"
ENV["MPLBACKEND"] = "Qt5Agg"
ClimatePlots.contourf(am)
ClimatePlots.show()
#using Plots
mapclimgrid(C,"World")

ClimatePlots.plot(C)
#using PyPlot
#show()
#C = load(fn,"hu";data_units="%",start_date=parsetuple("2014"),end_date=parsetuple("2015"))
using Dates
t = Date(2014, 1, 31)

t=(2014, 1, 31)

#https://juliaclimate.github.io/ClimateTools.jl/stable/biascorrection/
#Quantile-quantile mapping (Themeßl et al. 2012, Grenier et al. 2015) is provided with ClimateTools.jl through the function qqmap.

qqmap(obs::ClimGrid, ref::ClimGrid, 
fut::ClimGrid; method::String="Additive", 
detrend::Bool=true, window::Int=15, 
rankn::Int=50, thresnan::Float64=0.1, 
keep_original::Bool=false, interp = Linear(), extrap = Flat())

using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)

function loadcdo()
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end
loadcdo()
cd("/mnt/d/remo/genRE")
fn=readdir()|>first
#ra = tt.load(fn,"pr")
lpt="/mnt/d/remo/genRE"
run(`cdo -v -settime,12:00 -daysum $(joinpath(lpt, "genRE_precipitation_hour_1.1.nc")) $(joinpath(lpt, "genRE_daily.nc"))`)

ssup()
@vv "yfirst"
@vv "extract"
@vv "pnts ="

fn="genRE_daily.nc"
run(`cdo griddes $fn`)

# using Conda
# Conda.add("pyproj")
using PyCall
@pyimport rioxarray as rx
@pyimport geopandas as gpd
@pyimport xarray as xr
cd("/mnt/d/remo/genRE")
pt="/mnt/d/Wasim/main_basin.geojson"
g = gpd.read_file(pt)
z = g.to_crs(crs="EPSG:32632")
#gpd.GeoDataFrame.to_json("main_basin_32632.geojson")
z.to_file("main_basin_32632.geojson", driver="GeoJSON")

function dsop(dataset_path, geojson_path;epsg=32632)
    # Open the dataset
    ds = xr.open_dataset(dataset_path)
    ds = ds.rio.write_crs(epsg) #ds, "EPSG:32632"
    println("$epsg assinged")

    # Read the GeoDataFrame from the GeoJSON file
    gdf = gpd.read_file(geojson_path)
    try
        # # Check if the GeoDataFrame and the dataset have the same CRS
        # if gdf.crs != ds.rio.crs
        #     # Reproject the GeoDataFrame to match the CRS of the dataset
        #     gdf = gdf.to_crs(ds.rio.crs)
        # end
        # Crop the dataset using the bounding box of the GeoDataFrame
        ds_cropped = ds.rio.clip(gdf.geometry)        
        # Reproject the dataset to EPSG:4326
        println("Reprojecting the dataset to EPSG:4326")
        ds_cropped = ds_cropped.reproject(crs="EPSG:4326")
        ds_cropped.to_netcdf("genRE_hrs_crop.nc")
        println("genRE_hrs_crop.nc written...")

        println("Resampling to daily ts")
        ds_daily = ds_cropped.resample(time="D").sum("time")
        # Save the resampled data to a new file
        ds_daily.to_netcdf("genRE_daily_4326.nc")
        println("genRE_daily_4326.nc written...")
        return ds_daily
    catch
        #throw(gdf)
        @error "Error occurred during cropping."
    end
end

# Example usage:
dataset_path = "genRE_precipitation_hour_1.1.nc"
###dataset_path = "genRE_daily.nc"
geojson_path = "main_basin_32632.geojson"
cropped_dataset = dsop(dataset_path, geojson_path)




ds = xr.open_dataset(dataset_path)
ds = ds.rio.write_crs("EPSG:32632")
ds.crs
gdf = gpd.read_file(geojson_path)
#subset ds by year > 1999
ds = ds.sel(time=ds.time.dt.year > 1999)

ds["pr"].mean()

using PyPlot
pygui(true)
ds["pr"].isel(time=1).plot()
ds["pr"].mean("x").mean("y").groupby("time.year").plot()

# ds["pr"]
# ds.sel(time=ds.time.dt.year == 2000).mean("x")

#.to_netcdf("genRE_precipitation_hour_2000.nc")

#@doc ds.rio.clip

cropping_geometries = [g for g in gdf.geometry]
ds_cropped = ds.rio.clip(cropping_geometries, crs=32632)
ds_cropped = ds.rio.clip(gdf.geometry, crs=32632)

ds_chunk = ds.isel(time=90:100)
ds_cropped_chunk = ds_chunk.rio.clip(gdf.geometry, crs=32632)

try
        # # Check if the GeoDataFrame and the dataset have the same CRS
        # if gdf.crs != ds.rio.crs
        #     # Reproject the GeoDataFrame to match the CRS of the dataset
        #     gdf = gdf.to_crs(ds.rio.crs)
        # end
        # Crop the dataset using the bounding box of the GeoDataFrame
                
        # Reproject the dataset to EPSG:4326
        println("Reprojecting the dataset to EPSG:4326")
        ds_cropped = ds_cropped.reproject(crs="EPSG:4326")
        ds_cropped.to_netcdf("genRE_hrs_crop.nc")
        println("genRE_hrs_crop.nc written...")

        println("Resampling to daily ts")
        ds_daily = ds_cropped.resample(time="D").sum("time")
        # Save the resampled data to a new file
        ds_daily.to_netcdf("genRE_daily_4326.nc")
        println("genRE_daily_4326.nc written...")
        return ds_daily



cd("/mnt/d/remo/genRE")
ssup()
@pj
pyjl.xrp()
@pyimport xarray as xr
ds = xr.open_dataset("daily_4326.nc")
ds.crs
a = ds.sel(time=ds.time.dt.year < 2016).resample(time="Y",skipna=true).sum()
a = a.rename_vars(Dict("daily_4326"=>"pr"))
a["pr"].mean("longitude").mean("latitude").plot()

da = ds.sel(time=ds.time.dt.year < 2016).rename_vars(Dict("daily_4326"=>"pr","longitude" => "x", "latitude" => "y"))
@pyimport pandas as pd
vgjl("M8[")
#add 1 hour to time
da["time"] = da["time"] + pd.Timedelta("1 hour")




using PyPlot
a = ds["daily_4326"].mean("longitude").mean("latitude").groupby("time.year")

mpt = "/mnt/c/users/chs72fw/documents/EFRE_GIS/Hydrologie/mhm/ubu/"
fn = "test_domain/ConfigFile.log"
da = CSV.read(joinpath(mpt,fn),DataFrame,header=false)
pwd()
@edit op()
op()
pwc()




