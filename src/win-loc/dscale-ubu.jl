"/mnt/d/remo/qm/qm/"|>cd

#before activating the env, one could try...
#pt="/mnt/c/Users/chs72fw/.julia/dev/WaSiM/src/wa.jl"
#@time include(pt)

## in ubuntu, its the base environment, so this is not necessary
using Pkg
Pkg.activate(".")
Pkg.status()
Pkg.update()
using Base.Threads
nthreads()
#ENV["PROJ_LIB"] = "c:/OSGeo4W64/bin"
#ENV["PROJ_LIB"]
#using Conda
#Conda.ROOTENV
#### daily mean temperature TG######
cdb()
fdi()
wrkdir="/mnt/d/remo/qm/tas"
cd(wrkdir)

xr = pyimport("xarray")
pd = pyimport("pandas")
np = pyimport("numpy")


obslnk="/mnt/d/remo/cordex/eobs/tg_ens_spread_crop.nc"
simlnk="/mnt/d/remo/cordex/wgs/proj_tas_hist+rcp85_wgs84.nc"

sim = xr.open_dataset(simlnk)
obs = xr.open_dataset(obslnk)

using PyCall
pybuiltin = pyimport("builtins")
slice = pybuiltin["slice"]
# xr.pandas.slice
# xd = obs.where(time>"1970-01-01")
obsh = obs.sel(time=slice("1970-01-01","1990-01-01"))
obsp = obs.sel(time=slice("1991-01-01","2017-01-01"))

obsh["tg"].mean("longitude").mean("latitude").groupby("time.year").mean().plot()
obsh.to_netcdf("tg-obsh.nc")
obsp.to_netcdf("tg-obsp.nc")
pwd()

using ClimatePlots
using ClimateTools
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
PyCall.pygui()
mpl=pyimport("matplotlib")
pygui(true)
#first regrid...
tas = load("/mnt/d/remo/cordex/wgs/proj_tas_hist+rcp85_wgs84.nc","tas")
nconly(".")
#tg="/mnt/d/remo/cordex/eobs/tg_enc_nck_crop.nc"
tg="/mnt/d/remo/cordex/eobs/tg_ens_spread_crop.nc"
tg = load(tg,"tg")
cm.plot(annualmean(tg))

obsh = load("/mnt/d/remo/qm/tas/tg-obsh.nc","tg")
obsp = load("/mnt/d/remo/qm/tas/tg-obsp.nc","tg")
#out = qqmap(obs,simh,simp;method="Additive")
#out = qqmap(tg,obsh,obsp;method="Additive") #errors!

cm.plot(annualmean(tg),label = "tg")
cm.plot(annualmean(out) ,label = "bias")


#TX is needed!
#function temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)
Date(1970,1,1)|>typeof
#obsh = obs.sel(time=slice('1970-01-01','1990-01-01'))
#obsp = obs.sel(time=slice('1991-01-01','2017-01-01'))
#ob = tt.temporalsubset(tg, Date(1970,1,1), Date(1991,12,31))

sy=(1970,01,01)
ey=(1991,01,01)
fy=(2017,01,01)
ne=(2100,01,01)
#ob = tt.temporalsubset(tg, (1970,1,1),(1991,12,31))
ob = tt.temporalsubset(tg, sy,ey)
osp = temporalsubset(tg, ey,fy)
cm.plot(annualmean(ob))

rf = tt.temporalsubset(tas,sy,ey)
fut = tt.temporalsubset(tas,ey,ne)
cm.plot(annualmean(fut))

cm.plot(ob)
cm.plot(annualsum(ob))
cm.plot(annualsum(rf))


# #doesnt work
# pt=raw"C:\Users\Public\Documents\Python_Scripts\julia\win\ClimateTools_functions.jl"
# include(pt)
# using GeoStats
#xd = regrid(ob,rf)

##so i need to save the data first and regrid it via cdo...
#NCDatasets.renameDim is a Function.
#renameDim(ds::NCDataset, oldname::Union{AbstractString, Symbol}, newname::Union{AbstractString, Symbol})
#renameDim(ob,:tg,:tas)
pwd()
# ClimateTools.write(ob,"obsh.nc") #also errors :(
# tt.renameDim(ob,:tg=>:tas)

obs = load("/mnt/d/remo/qm/tas/obs.nc","tas") 
simh=load("/mnt/d/remo/qm/tas/simh.nc","tas")
simp=load("/mnt/d/remo/qm/tas/simp.nc","tas")
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

out = qqmap(obs,simh,simp;method="Additive")
cm.plot(annualmean(out))
cm.plot(annualmax(out))
cm.plot(annualmin(out))
cm.plot(annualmean(out),filename="tas-corr.png")

ClimateTools.write(out,"tas-cor.nc")

doy = dayofyear(out)
#rm("tas_cmethods.nc")
#tcm=load("tas_cm.nc","tas")
#ClimatePlots.contourf(out,region="EU",filename="tst.png");
#cm.plot(out,filename="tst.png");

#now mean wind speed FG
# 20	 E-OBS comes as an ensemble dataset and is available on a 0.1 and 0.25 degree 
# 21	 regular grid for the elements 
# 22	 daily mean temperature TG, 
# 23	 daily minimum temperature TN, 
# 24	 daily maximum temperature TX, 
# 25	 daily precipitation sum RR, 
# 26	 daily averaged sea level pressure PP, 
# 27	 daily averaged relative humidity HU, 
# 28	 daily mean wind speed FG 
# 29	 daily mean global radiation QQ.
ssup()
wrkdir="/mnt/d/remo/qm/wind"
cd(wrkdir)
ept = "/mnt/d/remo/cordex/eobs/"
#bash
obslnk="/mnt/d/remo/cordex/eobs/fg_ens_spread_crop.nc"
simlnk="/mnt/d/remo/cordex/wgs/proj_sfcWind_hist+rcp85_wgs84.nc"
cdo -v chname,fg,sfcWind -setcalendar,365_day -seldate,1970-01-01,2017-12-31 $obslnk obs.nc
cdo -v -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapbil,obs.nc $simlnk simh.nc
cdo -v -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs.nc $simlnk simp.nc
mrz

######
using ClimatePlots
using ClimateTools
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
PyCall.pygui()
mpl=pyimport("matplotlib")
pygui(true)       #wrks in ubu only.
pwd()

myvar="sfcWind"

obs = load(joinpath(wrkdir,"obs.nc") ,myvar) 
simh =load(joinpath(wrkdir,"simh.nc"),myvar)
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

#NCDatasets.cfvariable how to manually override the missing_value attribute.

out = qqmap(obs,simh,simp;method="Additive")
cm.plot(annualmean(out))
cm.plot(annualmax(out))
cm.plot(annualmin(out))
cm.plot(annualmean(out),filename="$myvar-corr.png")

ClimateTools.write(out,"$myvar-cor.nc")

vlat
xrcrds $lv 9.88 49.789 

using DataFrames
df = wread("sfcWind-cor.wa_ctlp",skip=1)
plot(df.date,df.x5)
dy=yrsum(df)

contourf(out)
obs=load("../prec/obsp.nc","pre")
contourf(obs)
pr=load("../prec/pre_delta_up2.nc","pre")
contourf(pr;cm="turbo")
#contourf(pr;region="Gr")

# Csub = resample(pr, "JJA")
# contourf(Csub;cm="turbo")

#now daily mean global radiation QQ.
wrkdir="/mnt/d/remo/qm/rsds"
cd(wrkdir)
ept = "/mnt/d/remo/cordex/eobs/"
#bash
"/mnt/d/remo/qm/rsds"
# obslnk="/mnt/d/remo/cordex/eobs/qq_ens_spread_crop.nc"
obslnk="/mnt/d/remo/cordex/eobs/qq_ens_spread_0.1deg_reg_v27.0e.nc"
simlnk="/mnt/d/remo/cordex/wgs/proj_rsds_hist+rcp85_wgs84.nc"
cdo -v chname,qq,rsds -setcalendar,365_day -seldate,1970-01-01,2017-12-31 $obslnk obs.nc
cdo -v -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapbil,obs.nc $simlnk simh.nc
cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs.nc $simlnk simp.nc

# mrz

myvar="rsds"
obs = load(joinpath(wrkdir,"obs.nc") ,myvar) 
obs = load(joinpath(wrkdir,"obs_fix.nc") ,myvar)
simh =load(joinpath(wrkdir,"simh.nc"),myvar)
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

contourf(obs)
contourf(simp)
out = qqmap(obs,simh,simp;)
cm.plot(annualsum(out))
cm.plot(annualmax(out))
contourf(out)
cm.plot(annualsum(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

ClimateTools.write(out,"$myvar-cor.nc")



wrkdir="/mnt/d/remo/qm/prec"
cd(wrkdir)
myvar="pre"
obs = load(joinpath(wrkdir,"obs.nc") ,myvar) 
simh =load(joinpath(wrkdir,"simh.nc"),myvar)
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

out = qqmap(obs,simh,simp;method="Multiplicative",thresnan=0.001)
cm.plot(annualsum(out))
cm.plot(annualmax(out))
contourf(out)
cm.plot(annualsum(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

ClimateTools.write(out,"$myvar-cor.nc")

vlat
xrcrds $lv 9.88 49.789

ssup()
using DataFrames
df = wread("$myvar-cor.wa_ctlp",skip=1)
dy=yrsum(df[1:end-1,:])
plot(dy.year,dy.x5)
plot(dy.year,dy.x5, "ro-", label="precsum", linewidth=2)
plot(dy.year,dy.x5, "rs", label="precsum", linewidth=2)
plot(dy.year,dy.x5, "^r:", label="precsum", linewidth=2)

k="/mnt/d/remo/qm/prec/Lohr-Main-Steinbach_pre+rcp85.tsv"
ld=wread(k,skip=5)

dy = yrsum(ld)
plot(dy.year,dy.x5, "^r:", label="precsum", linewidth=2)


cdb()
rglob(regand("ps","nc"))
wrkdir="/mnt/d/remo/qm/ps"
cd(wrkdir)
myvar="ps"

#bash
obslnk="/mnt/d/remo/cordex/eobs/pp_ens_spread_crop.nc"
#-1  pp            sea level pressure [hPa]

simlnk="/mnt/d/remo/cordex/wgs/proj_ps_hist+rcp85_wgs84.nc"
#-1  ps            Surface Air Pressure [mbar]

#sea level pressure [hPa] to Surface Air Pressure [mbar] via cdo ? 
cdo -v expr,'ps=pp/((1-0.0065*850/288.15)^5.25588)' $obslnk ps_obs.nc
vardes ps_obs.nc
griddes ps_obs.nc

#cdo -v chname,fg,sfcWind -setcalendar,365_day -seldate,1970-01-01,2017-12-31 $obslnk obs.nc
cdo -v -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapbil,ps_obs.nc $simlnk simh.nc
cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,ps_obs.nc $simlnk simp.nc
mrz

using NCDatasets
ds = NCDataset("ps_obs.nc")
ds = NCDataset("simh.nc")
ds.attrib
close(ds)

#obs = load(joinpath(wrkdir,"ps_obs.nc") , myvar)
obs = load(joinpath(wrkdir,"pp_ens_spread_crop.nc") , "pp")

simh =load(joinpath(wrkdir,"simh.nc"),myvar)
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

out = qqmap(obs,simh,simp;method="Multiplicative",thresnan=0.001)
cm.plot(annualsum(out))
cm.plot(annualmax(out))
contourf(out)
cm.plot(annualsum(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

ClimateTools.write(out,"$myvar-cor.nc")


##daily averaged relative humidity HU,  ##################
rglob(regand("ps","nc"))
wrkdir="/mnt/d/remo/qm/rh"
cd(wrkdir)
myvar="ps"

#bash
cd $wrkdir
obslnk="/mnt/d/remo/cordex/eobs/hu_ens_mean_crop.nc"
vardes $obslnk 
#-1  hu            mean relative humidity [%]         ##aah wrong!   
simlnk="/mnt/d/remo/cordex/wgs/proj_rh_hist+rcp85_wgs84.nc"
vardes $simlnk
#-1  rh            relative humidity [1]
cdo -v chname,hu,rh -setcalendar,365_day -seldate,1970-01-01,2017-12-31 $obslnk obs.nc
cdo -v -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapbil,obs.nc $simlnk simh.nc
cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs.nc $simlnk simp.nc
mrz
myvar="rh"
obs = load(joinpath(wrkdir,"obs.nc") , myvar)
simh =load(joinpath(wrkdir,"simh.nc"),myvar)
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmax(obs))
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

obs2 = obs*.1          #adjusted to 1/1

cm.plot(annualmean(obs2))
cm.plot(annualmean(simp))

out = qqmap(obs,simh,simp;method="Multiplicative",thresnan=0.001)
cm.plot(annualsum(out))
cm.plot(annualmax(out))
contourf(out)
cm.plot(annualsum(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

ClimateTools.write(out,"$myvar-cor.nc")


# p024x () 
# { 
#     perl -lni -e '@F=split(/\h/,$_);{$F[3]=~s/0/24/g;};print join("\t",@F)' "$1" && printf "changed "$1" to: \n$(head -7 "$1")\n"
# }
# axarg wa p024x


using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)
PyCall.pygui()
#mpl=pyimport("matplotlib")

s = "/mnt/d/remo/cordex/eobs/qq_ens_mean_v27_crop.nc"
s = "/mnt/d/remo/cordex/eobs/proj/qq_ens_spread_utm.nc"
s = "/mnt/d/remo/cordex/eobs/qq_ens_spread_0.1deg_reg_v27.0e.nc"
qg = load(s,"qq")
cm.plot(annualsum(qg))


pt="/mnt/d/remo/cordex/eobs/proj/tg_ens_spread_0.1deg_reg_v26.cropped.nc"
tas = load(pt,"tg")
cm.plot(annualmean(tas))
cm.plot(tt.monthmean(tas))

ts = tt.temporalsubset(tas, (1991,1,1), (1991,12,31))
cm.plot(ts)

N = [8.14986035840001, 48.8498605556]
tt.spatialsubset(ts,N)

# 8.14986035840001,48.8498605556
# gridtype  = lonlat
# gridsize  = 960
# xsize     = 48
# ysize     = 20

pwd()|>cb

myvar = "rsds"
rs = load("rstst.nc",myvar)
@pyimport xarray as xr
ob = xr.open_dataset("rstst.nc")
@pyimport matplotlib.pyplot as plt
display(ob)
mm = ob[myvar].groupby("time.month").mean()
mon = ob.where(ob["time.year"] > 2000, drop=true).groupby("time.month").mean("time")
pygui(true)
mon["rsds"].plot()
# Set the title and axis labels
plt.title("Monthly Mean RSDS")
plt.xlabel("Time")
plt.ylabel("RSDS")
plt.show()

x=ob.where(ob[myvar] > 0,drop=true)
#`...`
# g = x[myvar].groupby("time.dayofyear").mean()
# g.plot()
x.to_netcdf("obs_fix2.nc")

q = tt.load("obs_fix2.nc",myvar)
cm.plot(annualmean(q))
contourf(q)

simh =load(joinpath(wrkdir,"simh.nc"),myvar)
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

out = qqmap(obs,simh,simp;method="Multiplicative",thresnan=0.001)
cm.plot(annualsum(out))
cm.plot(annualmax(out))
contourf(out)
cm.plot(annualsum(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

########
"/mnt/d/remo/qm/corgrids/tas" |>cd
@time setup()
l=glob(r"cor.tsv$")
# dm = mall(l) #this errors!
dfs=[]
for x in l
    df = waread(x)
    push!(dfs,df)
end
first(dfs)
mo = select(first(dfs),2)
for i in dfs
    mo = hcat(mo,select(i,1))
end

first(mo)
op()
wawrite(mo,"tas-cor.wa")
#npp("tas-cor.wa")

println(mo|>propertynames)

@chain mo begin
    groupby(month.(mo.date), sort = true)
    transform(month.(mo.date) => ByRow(monthname) => :month_name)
end
    

pt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/temp_1970.txt"
nm = CSV.read(pt,DataFrame;limit=5)
dfm(pt;ann=false)
dfm(mo;ann=false)

wrkdir="/mnt/d/remo/qm/corgrids/rsds"
cd(wrkdir)
ssup()
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
outfile=tv*"_cor.wa"
using CSV
#@edit waread()
files = map(waread,glob(rgx))
#dfs = innerjoin(files..., on = :date, makeunique=true)
df = files[1] # assuming files is Vector
for i in 2:length(files)
    @show i
    global df = innerjoin(df, files[i], on = :date, makeunique=true)
end

ou = innerjoin(files[24], files[25], on = :date, makeunique=true)

nf = files[20:end]
nf = files[10:23]
ou = innerjoin(nf..., on = :date, makeunique=true)
k="/mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions/saale_airGR_v2.R"
pww(k)
using Conda
Conda.list()

Conda.CONDARC
Conda.ROOTENV
#2.769   GB
@edit fsize()
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
condasize()
# 8GB
Conda.clean()

using Pkg
"/mnt/d/remo/qm/qm/"|>Pkg.activate
pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
include(pt)
Pkg.status()
pwc()
"/mnt/d/remo/genRE"|>cd
#ncrename -O -v pr,pre genRE_main.nc
# @time_imports @pj
m=tt.load("genRE_main.nc","pre")
cm.plot(annualsum(m))
##tt.griddata
    cvar = Symbol("pre")
    # obspt = "/mnt/d/remo/cordex/eobs/v28/sfcWind/sfcWind-obs.nc"
    # obsras = tt.load(obspt,string(cvar))
    obsras = temporalsubset(m,(1997,01,01),(2015,01,01))
    cm.plot(annualsum(obsras))
    cm.contourf(obsras)
    
    simlnk="/mnt/d/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
    sim = tt.load(simlnk,string(cvar))
    @doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    simh = temporalsubset(sim,(1997,01,01),(2015,01,01))
    simp = temporalsubset(sim,(1997,01,01),(2099,01,01))
    simh = tt.griddata(simh,obsras) #2min
    simp = tt.griddata(simp,obsras) #ETA: 6 min


    #tt.write(simh,"simh-jl.nc")
    ###tt.write(simp,"simp-jl.nc") #errors mem-error.
        
    out = qqmap(obsras,simh,simp)
    tt.write(out,"$cvar-cor_raw.nc")


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


run(`conda env remove -p "/home/ubu/.julia/environments/v1.9/.CondaPkg/env"`)
#import Pkg; Pkg.add("EnvironmentMigrators")
using EnvironmentMigrators
EnvironmentMigrators.wizard()

import Pkg;
Pkg.status()
"/mnt/d/remo/qm/qm/"|>cd




using PyCall
@pyimport rioxarray as rx
@pyimport geopandas as gpd
@pyimport xarray as xr
cd("/mnt/d/remo/genRE")
pt="/mnt/d/Wasim/main_basin.geojson"
dataset_path = "genRE_precipitation_hour_1.1.nc"
###dataset_path = "genRE_daily.nc"
geojson_path = "main_basin_32632.geojson"
epsg=32632
ds = xr.open_dataset(dataset_path)
ds = ds.rio.write_crs(epsg) #
ds_daily = ds.drop_vars(["lat", "lon"]).resample(time="D").sum("time")
dl = ds.transpose("time", "y", "x").rio.reproject(dst_crs="EPSG:25832")
for var in dl.data_vars
    dl[var].attrs = delete!(dl[var].attrs, "grid_mapping")
end
# clip it.
gdf = gpd.read_file(geojson_path)
gdf = gdf.to_crs("EPSG:25832")
dl = dl.rio.clip(gdf.geometry) #no memerr.
# Save the dataset to a NetCDF file
dl.to_netcdf("genRE_25832.nc")


"/mnt/d/remo/genRE"
#cdx genRE_43|head     
cdx genRE_main_25832|head     
cdx genRE_main|head     
vardes genRE_main.nc
griddes genRE_main.nc

cdo -settaxis,1996-05-30,12:00:00,1day -setcalendar,365_day genRE_main.nc  gtmp.nc
#>>>>>>>>>>>>>>>moved to qqwin.jl due to memory issues <<<<<<<<<<<<<<<

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
m=tt.load("gtmp.nc","pre")
cm.plot(annualsum(m))
##tt.griddata
    cvar = Symbol("pre")
    
    obsras = temporalsubset(m,(1997,01,01),(2015,01,01))
    cm.plot(annualsum(obsras))
    cm.contourf(obsras)
    
    simlnk="/mnt/d/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
    sim = tt.load(simlnk,string(cvar))
    @doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    simh = temporalsubset(sim,(1997,01,01),(2015,01,01))
    simp = temporalsubset(sim,(1997,01,01),(2100,01,01)) #DateTime("2099-12-31T12:00:00")]
    simh = tt.griddata(simh,obsras) #2min
    simp = tt.griddata(simp,obsras) #ETA: 6 min


    #tt.write(simh,"simh-jl.nc")
    ###tt.write(simp,"simp-jl.nc") #errors mem-error.
        
    out = qqmap(obsras,simh,simp)
    tt.write(out,"$cvar-cor_raw.nc")


    max_obs = annualsum(temporalsubset(obsras,(1990,01,01),(2020,12,31)))
    #cm.plot(max_obs)
    max_modelinterp = annualsum(simp)
    max_modelqqmap = annualsum(out)
    begin
        cm.plot(max_obs, label="genRE")
        cm.plot(max_modelinterp, label="REMO - interpolated")
        cm.plot(max_modelqqmap, label="REMO - bias corrected", 
        titlestr = "Effect of bias correction on annual precipitation values",
        filename = "$cvar-mean.png")
    end
    lat()
    op()
    pwd()
    

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


##parse time from xarrays.
using PyCall
@pyimport xarray as xr
@pyimport datetime
@pyimport cftime
@pyimport pandas as pd
cd("d:/remo/genRE")
rglob("main")
pr = xr.open_dataset("genRE_main.nc")
tvec = map(x->pd.to_datetime(x), pr.time.values )
fn = "D:/remo/qm/corgrids/jlcor/pre-REcor2.nc"
pr = xr.open_dataset(fn)
#tvec = map(x->x.strftime(), pr.t.values )
fn = "D:/remo/qm/corgrids/jlcor/pre-REcor.nc" #<--corrupt.
pr = xr.open_dataset(fn)
tvec = map(x->x.strftime(), pr.time.values )
#parse(DateTime,first(tvec))
fmt = Dates.DateFormat("yyyy-mm-dd HH:MM:SS")   
tdat = [DateTime(d, fmt) for d in tvec]
dat = pr.pre.values;
for i in 1:size(dat)[3]
    println(findmax(dat[:,:,i]))
end 

@vv "contourf"
using Plots
#pr.pre.data[:,:,2]|>Plots.contourf
pr.pre.data[:,1,2]|>Plots.plot