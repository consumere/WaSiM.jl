"d:/remo/genRE"|>cd
using Pkg
"D:/remo/qm/qm/"|>Pkg.activate
Pkg.update()

Pkg.instantiate() #to install all recorded dependencies.
pt="C:/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
include(pt)
Pkg.status()
pwc()
using OhMyREPL

#ncrename -O -v pr,pre genRE_main.nc
# @time_imports @pj
# m=tt.load("genRE_main.nc","pre")

using PyCall
@pyimport rioxarray as rx
@pyimport geopandas as gpd
@pyimport xarray as xr
cd("D:/remo/genRE")
pfix()
#dataset_path = "genRE_precipitation_hour_1.1.nc"
dataset_path = "genRE_daily_32632.nc"
###dataset_path = "genRE_daily.nc"
geojson_path = "main_basin_32632.geojson"
epsg=32632
ds = xr.open_dataset(dataset_path)
ds = ds.rio.write_crs(epsg) #
ds = ds.drop_vars(["lat", "lon"])
for var in ds.data_vars
    ds[var].attrs = delete!(ds[var].attrs, "grid_mapping")
end
gdf = gpd.read_file(geojson_path)
ds = ds.rio.clip(gdf.geometry) #no memerr.

#ds_daily = ds.resample(time="D").sum("time")
glob("geojson")

dl = ds.transpose("time", "y", "x").rio.reproject(dst_crs="EPSG:25832")
# clip it.
#gdf = gdf.to_crs("EPSG:25832")
# Save the dataset to a NetCDF file
dl.pr.mean()

glob("main")
dl.to_netcdf("genRE_main_25832.nc")
pfix()



##mem errors. ->
ncf = "pre-REcor_raw.nc"
ds = xr.open_dataset(ncf)
    if "string1" in keys(ds.dims)
        ds = ds.drop_dims("string1")
    end
    if "Regular_longitude_latitude" in ds.variables
        ds = ds.drop_vars("Regular_longitude_latitude")
    end
    # Rename dimensions
    ds = ds.rename_dims(Dict("time" => "t"))
    # Reorder dimensions
    ds = ds.transpose("t", "x", "y")
    # Update attributes
    ds.x.attrs["units"] = "m East of reference point"
    ds.y.attrs["units"] = "m North of reference point"
    #cvar = string(cvar)
    ds.pre.attrs["missing_value"] = -9999.
    outfile = replace(ncf, "_raw.nc" => ".nc")
    #outpt = joinpath(pwd(),outfile)
    ds.dims
    ds.keys()
    outpt = outfile

    ds.crs.attrs
    ds = ds.drop_vars(["crs"])
    ds.pre.attrs
    ds.to_netcdf(outpt)
    println("$ncf done... \n saved to $outpt")
    ds.close()
    dumpcmd = `wsl ncdump -h $outfile`
    run(dumpcmd)



# cdo -settaxis,1996-05-30,12:00:00,1day -setcalendar,365_day genRE_main.nc  gtmp.nc

cdo -settaxis,1996-05-30,12:00:00,1day -setcalendar,proleptic_gregorian genRE_main.nc gtmp2.nc


function loadcdo()
    pt="C:/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end
loadcdo()

using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall


"D:/remo/genRE/"|>cd
m=tt.load("gtmp2.nc","pre")
cm.plot(annualsum(m))
##tt.griddata
    cvar = Symbol("pre")
    # obspt = "D:/remo/cordex/eobs/v28/sfcWind/sfcWind-obs.nc"
    # obsras = tt.load(obspt,string(cvar))
    obsras = temporalsubset(m,(1997,01,01),(2015,01,01)) #2016-01-20T12:00:00
    cm.plot(annualsum(obsras))
    cm.contourf(obsras)
    
    simlnk="D:/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
    sim = tt.load(simlnk,string(cvar))
    @doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    # simh = temporalsubset(sim,(1997,01,01),(2015,01,01))
    # simp = temporalsubset(sim,(1997,01,01),(2100,01,01))
    # simh = tt.griddata(simh,obsras) #Time: 0:01:23
    # simp = tt.griddata(simp,obsras) #ETA: 0:07:57

    #tt.write(simh,"simh-jl.nc")
    ###tt.write(simp,"simp-jl.nc") #errors mem-error.
    sim = tt.griddata(sim,obsras) #ETA: 0:11:34 #7122.90 MB
    pwd()
    tt.write(sim,"sim-jl.nc") #NO mem-error.
    #out = qqmap(obsras,simh,simp) #takes 1h
    out = qqmap(obsras,
        temporalsubset(sim,(1997,01,01),(2015,01,01)),
        temporalsubset(sim,(1997,01,01),(2100,12,30)))
    
    pwd()
    #tt.write(out,"$cvar-REcor_raw.nc") #4808.5367 MB
    tt.write(out,"$cvar-REcor_365d.nc") #MB


    max_obs = annualsum(temporalsubset(obsras,(1990,01,01),(2020,12,31)))
    #cm.plot(max_obs)
    max_modelinterp = annualsum(sim)
    max_modelqqmap = annualsum(out)
    begin
        cm.plot(max_obs, label="genRE")
        cm.plot(max_modelinterp, label="REMO - interpolated")
        cm.plot(max_modelqqmap, label="REMO - bias corrected", 
        titlestr = "Effect of bias correction on annual precipitation sums",
        filename = "$cvar-sum.png")
    end


cd("D:/remo/genRE")
fn="D:/remo/qm/prec/QDM_adjust.nc"
#ERROR: NetCDF error: Variable 'time_bnds' not found in file D:/remo/qm/prec/QDM_adjust.nc 
r = Raster(fn) #,key=:scen


