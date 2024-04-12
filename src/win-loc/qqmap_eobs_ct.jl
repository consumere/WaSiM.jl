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

##tt.regrid test. ==> griddata
begin 
    cvar = Symbol("tas")
    obspt = "/mnt/d/remo/cordex/eobs/v28/tas/tas_obs.nc"
    temp = tt.load(obspt,string(cvar))
    simlnk="/mnt/d/remo/cordex/wgs/proj_tas_hist+rcp85_wgs84.nc"
    sim = tt.load(simlnk,string(cvar))
    @doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    #"1990-01-01","2099-12-31"
    simh = temporalsubset(sim,(1990,01,01),(2020,01,01))
    simp = temporalsubset(sim,(1990,01,01),(2099,01,01))
    simh = tt.griddata(simh,temp)
    simp = tt.griddata(simp,temp) #7min
    out = qqmap(temp,simh,simp;method="Additive")


    max_obs = annualmean(temporalsubset(temp,(1990,01,01),(2020,12,31)))
    max_modelinterp = annualmean(simp)
    max_modelqqmap = annualmean(out)
    begin
        cm.plot(max_obs, label="EOBS")
        cm.plot(max_modelinterp, label="REMO - interpolated")
        cm.plot(max_modelqqmap, label="REMO - bias corrected", 
        titlestr = "Effect of bias correction on annual mean temperature values",
        filename = "$cvar-mean.png")
    end
    lat()
    op()
    #cm.plot(max_obs, label="EOBS")

    cm.plot(annualmean(obsras))

    pwd()
    #write out
    tt.write(out,"tas_cor_raw.nc")
    @doc tt.write
    @edit tt.write(out,"tas_cor.nc")
    latsymbol, lonsymbol = ClimateTools.getsymbols(out)
    # x, y, timevec = ClimateTools.getdims(C)
    # longrid, latgrid = ClimateTools.getgrids(C)
    @doc tt.renameDim()
    import NCDatasets
    const nc = NCDatasets

    versioninfo()
    using Pkg
    Pkg.installed()["NCDatasets"]

    pkgversion("aster")
    pkgversion("aster ")
    pkgversion("Frames")


    @doc NCDataset
    #ds = nc.NCDataset("tas_cor_raw.nc","c") #overwrites... not good.
    # #das geht, aber reorder nicht. daher lieber xarray.. 
    ds = nc.NCDataset("tas_cor_raw.nc","a") #append
    tt.renameDim(ds,latsymbol,"y")
    tt.renameDim(ds,lonsymbol,"x")
    tt.renameDim(ds,"time" , "t")
    #@doc tt.replace_missing!
    #todo .drop_vars(["Regular_longitude_latitude", "lat", "lon"])
    #reorder float tas(t, y, x) ;
    nc.close(ds)
    # pwc()
    pwin()

    #ValueError("cannot rename 'latitude' because it is not 
        #found in the dimensions of this dataset ('x', 'y', 'time')")

    xr
    tt.write(out,"tas_cor_raw.nc")
    ncf = "tas_cor_raw.nc"
    ds = xr.open_dataset(ncf)
    if "string1" in keys(ds.dims)
        ds = ds.drop_dims("string1")
    end
    if "Regular_longitude_latitude" in ds.variables
        ds = ds.drop_vars("Regular_longitude_latitude")
    end
    # Rename dimensions
    ds = ds.rename_dims(Dict("latitude" => "y", 
        "longitude" => "x", "time" => "t"))
    # Reorder dimensions
    ds = ds.transpose("t", "x", "y")
    # Update attributes
    ds.x.attrs["units"] = "m East of reference point"
    ds.y.attrs["units"] = "m North of reference point"
    ds.tas.attrs["missing_value"] = -9999.
    outfile = replace(ncf, "_raw.nc" => ".nc")
    outpt = joinpath(pwd(),outfile)
    ds.dims
    ds.keys()
    ds.to_netcdf(outpt)
    println("$ncf done... \n saved to $outpt")
    ds.close()
    # subprocess = pyimport("subprocess")
    # subprocess.call(["ncdump", "-h", outpt])
    dumpcmd = `ncdump -h $outpt`
    run(dumpcmd)
    pwin()
end


##tt.griddata => wind
cd("/mnt/d/remo/qm/corgrids/wind")
begin 
    cvar = Symbol("sfcWind")
    obspt = "/mnt/d/remo/cordex/eobs/v28/sfcWind/sfcWind-obs.nc"
    obsras = tt.load(obspt,string(cvar))
    simlnk="/mnt/d/remo/cordex/wgs/proj_sfcWind_hist+rcp85_wgs84.nc"
    sim = tt.load(simlnk,string(cvar))
    #@doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    #"1990-01-01","2099-12-31"
    simh = temporalsubset(sim,(1990,01,01),(2020,01,01))
    simp = temporalsubset(sim,(1990,01,01),(2099,01,01))
    simh = tt.griddata(simh,obsras) #2min
    simp = tt.griddata(simp,obsras) #7:48 min
    
    cm.plot(annualmean(obsras))
    out = qqmap(obsras,simh,simp)
    
    max_obs = annualmean(temporalsubset(obsras,(1990,01,01),(2020,12,31)))
    #cm.plot(max_obs)
    max_modelinterp = annualmean(simp)
    max_modelqqmap = annualmean(out)
    begin
        cm.plot(max_obs, label="EOBS")
        cm.plot(max_modelinterp, label="REMO - interpolated")
        cm.plot(max_modelqqmap, label="REMO - bias corrected", 
        titlestr = "Effect of bias correction on annual mean surface wind values",
        filename = "$cvar-mean.png")
    end
    lat()
    op()
    pwd()
    tt.write(out,"$cvar-cor_raw.nc")
    tt.write(simh,"simh-jl.nc") #errors
    
    using NCDatasets
    # Write data to NetCDF file
    NCDatasets.Dataset("simh-jl.nc", "c") do ds
        for name in keys(simh.data_vars)
            var = simh[name]
            defVar(ds, name, var.dtype, dimensions(var))
            ds[name][:] = var.values
        end
    end

    ncf = "$cvar-cor_raw.nc"
    ds = xr.open_dataset(ncf)
    if "string1" in keys(ds.dims)
        ds = ds.drop_dims("string1")
    end
    if "Regular_longitude_latitude" in ds.variables
        ds = ds.drop_vars("Regular_longitude_latitude")
    end
    # Rename dimensions
    ds = ds.rename_dims(Dict("latitude" => "y", 
        "longitude" => "x", "time" => "t"))
    # Reorder dimensions
    ds = ds.transpose("t", "x", "y")
    # Update attributes
    ds.x.attrs["units"] = "m East of reference point"
    ds.y.attrs["units"] = "m North of reference point"
    #cvar = string(cvar)
    ds.sfcWind.attrs["missing_value"] = -9999.
    outfile = replace(ncf, "_raw.nc" => ".nc")
    outpt = joinpath(pwd(),outfile)
    ds.dims
    ds.keys()
    ds.to_netcdf(outpt)
    println("$ncf done... \n saved to $outpt")
    ds.close()
    # subprocess = pyimport("subprocess")
    # subprocess.call(["ncdump", "-h", outpt])
    dumpcmd = `ncdump -h $outpt`
    run(dumpcmd)
    pwin()
end


##tt.griddata => rh --errors on qqmap
# need: group='time.dayofyear'
begin 
    cvar = Symbol("rh")
    obspt = "/mnt/d/remo/cordex/eobs/v28/rh/rh-obs1.nc"
    obsras = tt.load(obspt,string(cvar))
    simlnk="/mnt/d/remo/cordex/wgs/proj_rh_hist+rcp85_wgs84.nc"
    sim = tt.load(simlnk,string(cvar))
    #@doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    #"1990-01-01","2099-12-31"
    simh = temporalsubset(sim,(1990,01,01),(2020,01,01))
    simp = temporalsubset(sim,(1990,01,01),(2099,01,01))
    simh = tt.griddata(simh,obsras) #2min
    simp = tt.griddata(simp,obsras) #7:48 min
    

    obspt = "/mnt/d/remo/cordex/eobs/v28/rh/rh-obs.nc"
    obsras = tt.load(obspt,string(cvar))
    cm.plot(annualmean(obsras))
    
    out = qqmap(obsras,simh,simp)
    #cm.plot(out)
    
    #out is empty. for rh-obs1.
    #cmethods from 2020
    out = load("/mnt/d/remo/cordex/eobs/v28/rh/rh_cor.nc",string(cvar))


    @pyimport cmethods
    eobspt="/mnt/d/remo/cordex/eobs/v28/rh/"
    obsh = xr.open_dataset(joinpath(eobspt,"rh-obs1.nc"))
    #tt.write(simh,"simh.nc")
    sim = xr.open_dataset(simlnk)
    #simp = xr.open_dataset("simp.nc")
    simp = sim
    start_date_h = "1990-01-01"
    end_date_h = "2020-01-01"
    #slice = pyimport("slice")

    py"$sim.head()"
    py"simh = $sim.sel(time=slice($start_date_h, $end_date_h))"
    #py"simh = $sim.sel(time=$start_date_h:$end_date_h)"
    myvar = "rh"
    # x0 = py"$obsh[$myvar][:,0,0]"
    # x1 = py"$sim.sel(time=slice('1950-01-01','2000-01-01'))[$myvar][:,0,0]"
    # x2 = py"$sim.sel(time=slice('1980-01-01','2099-12-31'))[$myvar][:,0,0]"
    # ls_result = cmethods.CMethods.linear_scaling(
    #     obs = x0,
    #     simh = x1,
    #     simp = x2,
    #     kind = "*" #*> not available for linear_scaling.
    # )
    #cm.plot(ls_result)

    obsh = xr.open_dataset(joinpath(eobspt,"rh-obs1.nc"))
    simh = py"$sim.sel(time=slice('1950-01-01','2000-01-01'))"
    simp = py"$sim.sel(time=slice('1980-01-01','2099-12-31'))"
    obsh.keys()
    obsh = obsh.rename(Dict("longitude"=>"lon",
    "latitude"=>"lat"))
    obsh.dims
    #x = obsh.reorder_levels(dim_order=["time","lat","lon"])
    simh = simh.rename(Dict("longitude"=>"lon",
    "latitude"=>"lat"))
    simp = simp.rename(Dict("longitude"=>"lon",
    "latitude"=>"lat"))

    # 3d = 2 spatial and 1 time dimension
    qdm_result = cmethods.CMethods.adjust_3d(
        method = "quantile_delta_mapping",
        #method = "quantile_mapping",
        obs = obsh[myvar],
        simh = simh[myvar],
        simp = simp[myvar],
        n_quaniles = 1000,
        detrended = true,
        kind = "*"
    )
    #[myvar]
    qdm_result.mean("time").plot()
    #ValueError("Dimensions {'lat', 'lon'} do not exist. Expected one or more of ('time', 'latitude', 'longitude')
    qdm_result.to_netcdf("rh_cor.nc")

    
    max_obs = annualmean(temporalsubset(obsras,(1990,01,01),(2020,12,31)))
    #cm.plot(max_obs)
    max_modelinterp = annualmean(simp)
    max_modelqqmap = annualmean(out)
    begin
        cm.plot(max_obs, label="EOBS")
        cm.plot(max_modelinterp, label="REMO - interpolated")
        cm.plot(max_modelqqmap, label="REMO - bias corrected", 
        titlestr = "Effect of bias correction on annual mean temperature values",
        filename = "$cvar-mean.png")
    end
    lat()
    op()

    #ValueError("cannot rename 'latitude' because it is not 
        #found in the dimensions of this dataset ('x', 'y', 'time')")

    xr
    tt.write(out,"tas_cor_raw.nc")
    ncf = "tas_cor_raw.nc"
    ds = xr.open_dataset(ncf)
    if "string1" in keys(ds.dims)
        ds = ds.drop_dims("string1")
    end
    if "Regular_longitude_latitude" in ds.variables
        ds = ds.drop_vars("Regular_longitude_latitude")
    end
    # Rename dimensions
    ds = ds.rename_dims(Dict("latitude" => "y", 
        "longitude" => "x", "time" => "t"))
    # Reorder dimensions
    ds = ds.transpose("t", "x", "y")
    # Update attributes
    ds.x.attrs["units"] = "m East of reference point"
    ds.y.attrs["units"] = "m North of reference point"
    ds.tas.attrs["missing_value"] = -9999.
    outfile = replace(ncf, "_raw.nc" => ".nc")
    outpt = joinpath(pwd(),outfile)
    ds.dims
    ds.keys()
    ds.to_netcdf(outpt)
    println("$ncf done... \n saved to $outpt")
    ds.close()
    # subprocess = pyimport("subprocess")
    # subprocess.call(["ncdump", "-h", outpt])
    dumpcmd = `ncdump -h $outpt`
    run(dumpcmd)
    pwin()
end


cd("/mnt/d/remo/cordex/eobs/v28/")
fdi()
rglob("rsds")
##tt.griddata => rsds
cd("/mnt/d/remo/qm/corgrids/rsds")
#begin 
    cvar = Symbol("rsds")
    obspt = "/mnt/d/remo/cordex/eobs/v28/rsds/rsds_eobs.nc"
    obsras = tt.load(obspt,string(cvar))
    simlnk="/mnt/d/remo/cordex/wgs/proj_rsds_hist+rcp85_wgs84.nc"
    sim = tt.load(simlnk,string(cvar))
    #@doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    #"1990-01-01","2099-12-31"
    # simh = temporalsubset(sim,(1990,01,01),(2020,01,01))
    # simp = temporalsubset(sim,(1990,01,01),(2099,01,01))
    # simh = tt.griddata(simh,obsras) #2:50 min
    # simp = tt.griddata(simp,obsras) #10:48 min

    simh = temporalsubset(sim,(1950,01,01),(2015,01,01))
    simp = temporalsubset(sim,(2015,01,01),(2099,01,01))
    simh = tt.griddata(simh,obsras) 
    simp = tt.griddata(simp,obsras) 

    out = qqmap(obsras,simh,simp)
    
    cm.plot(annualmean(obsras))
    #max_obs = annualmean(temporalsubset(obsras,(1990,01,01),(2020,12,31)))
    #cm.plot(max_obs)
    max_obs = annualmean(obsras)
    max_modelinterp = annualmean(simp)
    max_modelqqmap = annualmean(out)
    begin
        cm.plot(max_obs, label="EOBS")
        cm.plot(max_modelinterp, label="REMO - interpolated")
        cm.plot(max_modelqqmap, label="REMO - bias corrected", 
        titlestr = "Effect of bias correction on annual radiation values",
        filename = "$cvar-mean-1950.png")
    end

    lat()
    op()
    
    pwd()
    ls()
    tt.write(out,"$cvar-cor_raw.nc")
       
    ##redim for wasim
    ncf = "$cvar-cor_raw.nc"
    ds = xr.open_dataset(ncf)
    if "string1" in keys(ds.dims)
        ds = ds.drop_dims("string1")
    end
    if "Regular_longitude_latitude" in ds.variables
        ds = ds.drop_vars("Regular_longitude_latitude")
    end
    # Rename dimensions
    ds = ds.rename_dims(Dict("latitude" => "y", 
        "longitude" => "x", "time" => "t"))
    # Reorder dimensions
    ds = ds.transpose("t", "x", "y")
    # Update attributes
    ds.x.attrs["units"] = "m East of reference point"
    ds.y.attrs["units"] = "m North of reference point"
    #cvar = string(cvar)
    ds.rsds.attrs["missing_value"] = -9999.
    #outfile = replace(ncf, "_raw.nc" => ".nc")
    outfile = replace(ncf, "_raw.nc" => ".nc")
    cd("/mnt/d/remo/cordex/eobs/v28/rsds")
    outpt = joinpath(pwd(),outfile)
    ds.dims
    ds.keys()
    ds.to_netcdf(outpt)
    println("$ncf done... \n saved to $outpt")
    ds.close()
    # subprocess = pyimport("subprocess")
    # subprocess.call(["ncdump", "-h", outpt])
    run(`ncdump -h $outpt`)
    pwin()
end


lo = xr.open_dataset("rsds-cor.nc")
ls = xr.open_dataset("rsds_cor.nc")

# lo["rsds"].mean("x").mean("y").plot()
# ls["rsds"].mean("x").mean("y").plot()

begin
    lo["rsds"].mean("x").mean("y").groupby("time.year").mean().plot(label="bias corrected 1950 - 2099")
    ls["rsds"].mean("longitude").mean("latitude").groupby("time.year").mean().plot(label="old method")
end



@doc qqmap
#obsvec::Array{N,1}, refvec::Array{N,1}, futvec::Array{N,1}
ov = temporalsubset(obsras,(1990,01,01),(2020,12,31))
rv = temporalsubset(obsras,(1950,01,01),(2020,12,31))
#fv = temporalsubset(sim,(1980,01,01),(2099,12,31)) #griddata necessary
fv = simp
#out = qqmap(obsras,simh,simp)
size(rv[1],1)
size(fv[1],1)
#bias = qqmap(ov,rv,simp;method="Additive") #weil 2xobsgrids dublikate vorhanden. julia syn-error.


###xclim test
using PyCall
cd("/mnt/d/remo/cordex/eobs/v28/tas")
#using Conda
#Conda.add("xclim")
@pyimport xclim as xc
@pyimport xarray as xr
@pyimport matplotlib.pyplot as plt

#ds = xr.open_dataset("tas_cor.nc")
#tg = xc.atmos.tg_mean(ds.tas, freq="MS")
#tg.mean("longitude").mean("latitude").mean("time").plot()
#plt.show()
@pyimport xclim.sdba as sdba
@pyimport nc_time_axis
#@pyimport builtins.slice as slice
xbuiltin = pyimport("builtins")
slice = xbuiltin["slice"]
ls()
ref = xr.open_dataset("tas_obs.nc").sel(time=slice("2015-01-01","2017-01-01"))
hist = xr.open_dataset("simh.nc").sel(time=slice("2015-01-01","2017-01-01"))
sim = xr.open_dataset("simp.nc") #.tas

# hist = hist.sel(time=slice("2015-01-01","2017-01-01"))
# sim = sim.sel(time=slice("2015-01-01","2017-01-01"))
sim.keys()

QM = sdba.EmpiricalQuantileMapping.train(
    ref.tas, hist.tas, nquantiles=15, group="time", kind="+"
)

scen = sdba.EmpiricalQuantileMapping.adjust(QM, extrapolation="constant", interp="nearest")

 ref["tas"].groupby("time.dayofyear").mean("longitude").mean("latitude").plot(label="Reference")
hist["tas"].groupby("time.dayofyear").mean("longitude").mean("latitude").plot(label="Model - biased")
plt.show()

# sim.sel(time=slice("2000", "2015")).groupby("time.dayofyear").mean().plot(
#     label="Model - adjusted - 2000-15", linestyle="--"
# )
# scen.sel(time=slice("2015", "2030")).groupby("time.dayofyear").mean().plot(
#     label="Model - adjusted - 2015-30", linestyle="--"
# )
# plt.legend()

tasjl = xr.open_dataset("tas_cor.nc")   
taspy = xr.open_dataset("simh_adj20002020.nc")
tasjl["tas"].groupby("time.dayofyear").mean("longitude").mean("latitude").plot(label="JL corrected")
taspy.groupby("time.dayofyear").mean("longitude").mean("latitude").plot(label="PY corrected")
plt.show()


function loadcdo()
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end
loadcdo()

cd("/mnt/d/remo/qm/corgrids/jlcor/")
fn="sfcWind-cor.nc"
win = tt.load(fn,"sfcWind")
cm.contourf(win)
using PyCall
@pyimport xclim as xc
@pyimport xarray as xr
@pyimport matplotlib.pyplot as plt
fn="f_sfcWind.nc"
ds = xr.open_dataset(fn)

run(`cdo vardes $fn`)
#run(`cdo vardes pre-cor.nc`)

######genRE qqmap
cd("/mnt/d/remo/genRE/")
using Pkg
Pkg.activate("/mnt/d/remo/qm/qm/")
Pkg.status()
function loadcdo()
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end

loadcdo()
readdir()
@time_imports ssup()
df=dfr("prec_ext.wa")
@time_imports @pyjl
pyjl.pyplot_df(df)
pyjl.pybox(df)
glob(r"winlist")

fn="pr-dly.winlist"
readlines(fn)[1:10]|>x->replace.(x, "\t" => " ")
#readlines(fn)[5:10] .|> x -> replace(x, "\\t" => "\t")

pwc()

dl = tt.load("daily_4326.nc","pr")
using PyCall
@pyimport xclim as xc
@pyimport xarray as xr
@pyimport matplotlib.pyplot as plt
#fn="daily_4326.nc"
fn="pr-dly.nc"
ds = xr.open_dataset(fn)
fn="genRE_precipitation_hour_1.1.nc"
fn="genRE_daily_32632.nc"

dl = tt.load(fn,"pr")

@pyimport xclim.indices as xci
p = xr.open_dataset(fn).pr
p.attrs["units"] = "mm/day"
#pday_seasonality = xci.precip_seasonality(p)
p_month = xci.precip_accumulation(p, freq="30D")
p_month.attrs["units"] = "mm/month"
monseas = xci.precip_seasonality(p_month)


p_weekly = xci.precip_accumulation(p, freq="7D")
#put units need to be a rate
p_weekly.attrs["units"] = "mm/week"
pweek_seasonality = xci.precip_seasonality(p_weekly)

#bash

#cp -v daily_4326.nc dly.nc
ncatted -O -a grid_mapping_name,global,a,c,"EPSG:4326" dly.nc

fn = "dly.nc"
fn = "/mnt/d/remo/cordex/eobs/v28/pre/simh.nc"
#dl = tt.load(fn,"pr")
dl = tt.load(fn,"pre")
ds = tt.monthsum(dl)
cm.plot(ds, label="EOBS Monthly sum")


macro css_str(x::String)
    kvalues::Vector{SubString} = split(x, ";")
    Dict{Symbol, String}(begin
         splt = split(kval, ":")
         Symbol(splt[1]) => string(splt[2])
    end for kval in kvalues)::Dict{Symbol, String}
end

css"color:blue; background-color:yellow"
css"border-radius:5px;color:red;width:500px"


glob("genRE")
#dl = tt.load("genRE_precipitation_hour_1.1.nc","pr")

using PyCall
@pyimport xclim as xc
@pyimport xarray as xr
@pyimport matplotlib.pyplot as plt

ds = xr.open_dataset("genRE_precipitation_hour_1.1.nc")

cd("/mnt/d/remo/qm/corgrids/wind")
cvar = Symbol("sfcWind")
ds = tt.load("sfcWind-cor.nc",string(cvar))
ls()
ds = tt.load("sfcWind-cor_raw.nc",string(cvar))
cm.contourf(ds)
cm.rc("font", family="serif", serif=["cmr10"])
# Set use_mathtext property to True
plt.gcf().axes.formatter.use_mathtext = true
#cm.gcf().axes.formatter.use_mathtext = true
cm.plot(annualmax(ds))

#vaporpressure(specifichumidity::ClimGrid, sealevelpressure::ClimGrid, orography::ClimGrid,
#  daily_temperature::ClimGrid)
#Returns the vapor pressure (vp) (Pa) estimated with the specific humidity (q), the sea level
#  pressure (psl) (Pa), the orography (orog) (m) and the daily mean temperature (tas) (K).
# vp = \frac{q * sp}{q+0.622}
cvar = Symbol("sfcWind")
obspt = "/mnt/d/remo/cordex/eobs/v28/sfcWind/sfcWind-obs.nc"
obsras = tt.load(obspt,string(cvar))
simlnk="/mnt/d/remo/cordex/wgs/proj_sfcWind_hist+rcp85_wgs84.nc"
sim = tt.load(simlnk,string(cvar))
#@doc griddata
#Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
#where A and B are ClimGrid.
#"1990-01-01","2099-12-31"
simh = temporalsubset(sim,(1970,01,01),(2020,01,01))
simp = temporalsubset(sim,(1970,01,01),(2100,01,01))    #(2099,12,31)) #DateTime("2099-12-31T12:00:00")]
simh = tt.griddata(simh,obsras) #3 min
simp = tt.griddata(simp,obsras) #7 min

#64.745055 seconds (18.51 M allocations: 22.994 GiB, 1.87% gc time, 48.00% compilation time)
@time out = qqmap(obsras,simh,simp)
tt.write(out,"wind_jlcor.nc")

#max_obs = annualmean(temporalsubset(obsras,(1990,01,01),(2020,12,31)))
max_obs = annualmax(obsras)
max_modelinterp = annualmax(simp)
max_modelqqmap = annualmax(out)
begin
    cm.plot(max_obs, label="EOBS")
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum surface wind values",
    filename = "$cvar-max.png")
end


simh.dataunits
obsras.dataunits
#:time, [DateT

#tt.write(simh,"simh.nc") #errors.

#ok = tt.groupby(out,"time.dayofyear").mean()

#Combine the empirial Quantile-Quantile mapping (see qqmap) and Generalized Pareto Distribution bias-correction methods.
outextr = tt.biascorrect_extremes(obsras,simh,simp)
contourf(outextr)
contourf(out)
tt.write(outextr,"sfcWind-cor_extremes.nc")
begin
    cm.plot(max_obs, label="EOBS")
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(annualmax(outextr), label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum extremes surface wind values",
    filename = "$cvar-max-extremes.png")
end

###########now extremes biascor for precipitation #####################obsras = 
using Pkg
Pkg.activate("/mnt/d/remo/qm/qm/")
function loadcdo()
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end

loadcdo()
ssup()
cvar = Symbol("pre")
cd("/mnt/d/remo/qm/corgrids/pre")
#obspt = "/mnt/d/remo/genRE/genRE_precipitation_hour_1.1.nc"
#obspt = "/mnt/d/remo/genRE/genRE_main.nc"
#fn = "/mnt/d/remo/genRE/pre-daily.nc"
#fn = "/mnt/d/remo/genRE/genRE_daily_32632.nc"
#fn = "/mnt/d/remo/genRE/genRE_daily.nc"
obspt = "/mnt/d/remo/cordex/eobs/v28/pre/pre_rcm_obs.nc"
obsout = "obs_365d.nc"
run(pipeline(`cdo -v -setcalendar,365_day $obspt $obsout`))
simlnk="/mnt/d/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
run(pipeline(`cdo -v -setcalendar,365_day -seldate,1970-01-01,2020-01-01 -remapbil,$obsout $simlnk simh1.nc`))
run(pipeline(`cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,$obsout $simlnk simp1.nc`))

#obsras = tt.load(obspt,"pre") #string(cvar)
obsras = tt.load(obsout,"pre") 
#"1996-05-30T00:00:00"..."2016-01-16T00:00:00"
#obsras = temporalsubset(obsras,(1997,01,01),(2016,01,01))
contourf(obsras)
#sim = tt.load(simlnk,string(cvar))
simh = tt.load("simh1.nc",string(cvar))
simp = tt.load("simp1.nc",string(cvar))
# simh = temporalsubset(sim,(1997,01,01),(2016,01,01))
# simp = temporalsubset(sim,(1997,01,01),(2100,01,01))
# simh = temporalsubset(sim,(1970,01,01),(2020,01,01))
# simp = temporalsubset(sim,(1970,01,01),(2100,01,01)) 
#simp = temporalsubset(sim,(1997,01,01),(2050,01,01))
simh = tt.griddata(simh,obsras) #2 min
simp = tt.griddata(simp,obsras) #10 min

#out = qqmap(obsras,simh,simp) #julia breaks.
#@time outextr = tt.biascorrect_extremes(obsras,simh,simp) #also breaks on genRE.
#73.015514 seconds (49.13 M allocations: 37.786 GiB, 33.58% gc time, 80.00% compilation time)

#41.886294 seconds (46.73 M allocations: 37.221 GiB, 4.92% gc time, 76.21% compilation time)
@time outextr = tt.biascorrect_extremes(obsras,simh,simp) #
contourf(outextr)
max_obs = annualmax(obsras)
max_modelinterp = annualmax(simp)
max_modelqqmap = annualmax(outextr)
begin
    cm.plot(max_obs, label="EOBS")
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum precipitation values",
    filename = "$cvar-max2.png")
end
cp("pre-cor_extremes.nc","pre-cor_extremes_366d.nc")
tt.write(outextr,"pre-cor_extremes.nc")
pwc()
mns = tt.monthsum(obsras)
simn = tt.monthsum(simp)
outn = tt.monthsum(outextr)

begin
    cm.plot(mns, label="EOBS")
    cm.plot(simn, label="REMO - interpolated")
    cm.plot(outn, label="REMO - bias corrected")
    #titlestr = "Effect of bias correction on monthly sum precipitation values",
    #filename = "$cvar-monthsum.png")
end

bia = simn - outn
cm.plot(bia, label="REMO - bias")
cm.plot(annualsum(bia), label="REMO - bias")
cm.plot(annualsum(outn), label="bias corrected")

# <line level="4" description= "       routing_model 28.2.2012 12. h (time of meteo data)"/>
# <line level="4" description= "         temperature 29.2.2012 12. h "/>
# <line level="0" description= "--;error in date columns of precipitation;;first Meteo data set (this is the set from which data is get:;29.2.2012 12:0 h;;but in precipitation:;1.3.2012 12:0 h;;--;"/>
# <line level="0" description= "Exception in Modellapplication::run() caught!"/>
pwc()
#/mnt/d/remo/qm/corgrids/pre
###########make wasim gridlist#######################
cd("/mnt/d/remo/qm/corgrids/pre")
a="pre-cor_extremes.nc"
open("ts2", "w") do out
    #ts_content = String(read(pipeline(`wsl cdo showtimestamp $a`)))
    ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n") #s/ /\n/g
    #ts_content = replace(ts_content, r"\s+" => "\t")
    ts_content = replace(ts_content, r"-|T|:" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "\n" => "", count=1)
    #ts_content = replace(ts_content, r"\s+" => "\t")
    write(out, ts_content)
end

mv("pre-cor_extremes.winlist","pre-cor_extremes_366d.nclist")
function process_files(match::Regex,suffix=".winlist")
    files = filter(x -> occursin(match,x) && endswith(x, ".nc"), readdir())
    # if any(x -> occursin(match, x) && endswith(x, suffix), readdir())
    #     println("Corresponding files already exist. Exiting now!")
    #     return
    # end
    for var in files
        gridfile = string(replace(var,".nc" => ".gridfile"))
        global onam = string(replace(var,".nc" => suffix))

        open(gridfile, "w") do out
            write(out, String(read(pipeline(`cdo griddes $var`))))
        end
        

        open(onam, "w") do file
            println(file, "GRID LIST of $var")
            println(file, "YY\tMM\tDD\tHH\t9999")
            print(file, "YY\tMM\tDD\tHH\t")
            println(file, read_grid_value(gridfile, "xfirst"))
            print(file, "YY\tMM\tDD\tHH\t")
            println(file, read_grid_value(gridfile, "yfirst"))
            println(file, "YY\tMM\tDD\tHH\tListe")

            ts2 = "ts2"  # Replace this with the actual variable or data source
            lines = readlines(ts2)
            for (i, line) in enumerate(lines)
                global param = split(var, '-')[1]
                line = replace(line, "00" => string("\tD:/remo/qm/corgrids/jlcor/",var,"<",param,">",i-1))
                println(file, replace(line, r"\s+" => "\t"))
            end
        end

        println("$param written to $onam !")
        println("$onam done!")
    end

    println("Lists for WaSiM Input are written to the current directory. All done!")
end

function read_grid_value(gridfile, pattern)
    lines = readlines(gridfile)
    for line in lines
        if contains(line,pattern)
            match_obj = split(line,"=")|>last
            return isnothing(match_obj) ? "N/A" : strip(match_obj)
        end
    end
    return "N/A"
end
process_files(r"pre-cor_extremes+.*nc")
z=lat()
npp(z)
##fix for winlist
#cut -f1-5 tmp 
perl -i -F'\t' -lane 'print join("\t", @F[0..4])' pre-cor_extremes.winlist
#PS D:\remo\qm\corgrids\jlcor> cp -v "D:\remo\qm\corgrids\pre\pre-cor_extremes.nc" .
#PS D:\remo\qm\corgrids\jlcor> cp -v "D:\remo\qm\corgrids\pre\pre-cor_extremes.winlist" . 
cd("/mnt/d/remo/qm/corgrids/jlcor")

a="pre-cor_extremes.nc"
cp(a,"pre-cor_extremes_366d.nc")

# ncpdq -x -v lon,lat --rdr=time,x,y tmp.nc pre-cor_extremes.nc
# ncpdq -C -x -v lon,lat --rdr=time,x,y tmp.nc pre-cor_extremes.nc
# ncrename -O -d time,t pre-cor_extremes.nc
# vlat
# ncvarlst $lv
# #ncrename -O -d latitude,y -d longitude,x pre-cor_extremes.nc
# ncrename -O -v latitude,y -v longitude,x pre-cor_extremes.nc
# ncrename -O -v time,t pre-cor_extremes.nc
# ncvarlst $lv

# cp -v pre-cor_extremes.nc tmp2.nc
# ncpdq -x -v Regular_longitude_latitude --rdr=t,x,y tmp2.nc pre-cor_extremes.nc
# ncatted -O -a _FillValue,,o,f,-9999 pre-cor_extremes.nc
# #ncpdq -x -v Regular_longitude_latitude -v latitude -v longitude --rdr=t,x,y tmp2.nc pre-cor_extremes.nc

###fix for wasim --->only working method....!!!!
cd("/mnt/d/remo/qm/corgrids/jlcor")
cp -v "/mnt/d/remo/qm/corgrids/pre/pre-cor_extremes.nc" pre-cor_extremes_raw.nc
cp -v "/mnt/d/remo/qm/corgrids/pre/pre-cor_extremes.winlist" .
cp -v "/mnt/d/remo/qm/corgrids/pre/pre-cor_extremes.gridfile" .

using PyCall
###xr = pyimport("xarray")
@pyimport xarray as xr
ncf = "pre-cor_extremes_raw.nc"
ds = xr.open_dataset(ncf)
if "string1" in keys(ds.dims)
        ds = ds.drop_dims("string1")
end
if "Regular_longitude_latitude" in ds.variables
        ds = ds.drop_vars("Regular_longitude_latitude")
end
# Rename dimensions
ds = ds.rename_dims(Dict("time" => "t"))
ds = ds.rename_dims(Dict("longitude" => "x", "latitude" => "y"))
# Reorder dimensions
ds = ds.transpose("t", "x", "y")
# Update attributes
ds.x.attrs["units"] = "m East of reference point"
ds.y.attrs["units"] = "m North of reference point"
#cvar = string(cvar)
ds.pre.attrs["missing_value"] = -9999.
outfile = replace(ncf, "_raw.nc" => ".nc")
ds.dims
ds.keys()
outpt = outfile
#ds.crs.attrs
#ds = ds.drop_vars(["crs"])
ds.pre.attrs
ds.to_netcdf(outpt)
println("$ncf done... \n saved to $outpt")
ds.close()
dumpcmd = `ncdump -h $outfile`
run(dumpcmd)


cd("/mnt/d/remo/qm/corgrids/jlcor")
perl -pe 's;D:/remo/qm/corgrids/jlcor;stateini;g' pre-cor_extremes.winlist > pre-cor_extremes.nc_list
#cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs_fix2.nc $simlnk simp.nc
#"
f="d:/Wasim/Tanalys/DEM/brend_fab/in4/fab.tlowbdry.nc"
r=Raster(f)
r = r[t=1]
of="d:/Wasim/Tanalys/DEM/brend_fab/in4/lonlat.nc"
write(of,r)

