
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
#op()
##rcm crop .. small area
nconly("q")
ob = xr.open_dataset("qq_ens_mean_v28.0e_rcm.nc")
display(ob)
ob = ob.rename(Dict("qq"=>"rsds"))
x = ob.where(ob["rsds"] > 0,drop=true)
x.to_netcdf("rsds_rcm_obs.nc")
x.close()
#mit setcalendar,365_day dubplot(xd) #two dubs!!!!
##now without...
begin
    simlnk="/mnt/d/remo/cordex/wgs/proj_rsds_hist+rcp85_wgs84.nc"
    output_file = "simh.nc"
    #remapfile = "rsds_eobs.nc"
    remapfile = "rsds_rcm_obs.nc"
    ncremap(simlnk,remapfile,output_file,"1970-01-01","2020-12-31")
    simpfile = "simp.nc"
    ncremap(simlnk,remapfile,simpfile,"2021-01-01","2099-12-31")
end
simh = tt.load("simh.nc","rsds")
simp = tt.load("simp.nc","rsds")
obs  = tt.load(remapfile,"rsds")
out = qqmap(obs,simh,simp;)

cm.contourf(out;
    cm="magma",titlestr="qmap rsds")
tt.write(out,"rsds_cor.nc")
begin
    #max_obs = annualmax(temporalsubset(obs,(2000,01,01),(2000,01,01)))
    max_obs = annualmax(obs)
    max_modelinterp = annualmax(temporalsubset(simp,(2000,01,01),(2099,01,01)))
    max_modelqqmap = annualmax(temporalsubset(out,(2000,01,01),(2099,01,01)))
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
    cm.plot(max_obs, label="OBS")
end
r = Raster("rsds_cor.nc";key= "rsds")
xd = ncdf(r)
dubplot(xd) #no dubs!!!!

#now a test with rr
ras = glob("rcm")
#tt.load(ras[4],"rr")
# using NCDatasets
# ds = NCDatasets.read(ras[4])
# tt.renameVar(ds,"rr","pre")
ob = xr.open_dataset(ras[4])
ob = ob.rename(Dict("rr"=>"pre"))
#x = ob.where(ob["rr"] > 0,drop=true)
mkdir("pre")
cd("pre")
ob.to_netcdf("pre_rcm_obs.nc")
ob.close()

#
begin
    simlnk="/mnt/d/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
    output_file = "simh.nc"
    remapfile = "pre_rcm_obs.nc"
    ncremap(simlnk,remapfile,output_file,"1970-01-01","2020-12-31")
    simpfile = "simp.nc"
    ncremap(simlnk,remapfile,simpfile,"2021-01-01","2099-12-31")
end
simh = tt.load("simh.nc","pre")
simp = tt.load("simp.nc","pre")
obs  = tt.load(remapfile,"pre")

qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; 
method="Additive", detrend=true, window::Int=15, 
rankn::Int=50, thresnan::Float64=0.1, 
keep_original::Bool=false, 
interp::Function = Linear(), 
extrap::Function = Flat())

out = qqmap(obs,simh,simp;method="Multiplicative",thresnan=0.001)

cm.contourf(out;
    cm="magma",titlestr="pre rsds")
tt.write(out,"pre_cor.nc")
begin
    #max_obs = annualmax(temporalsubset(obs,(2000,01,01),(2000,01,01)))
    max_obs = annualmax(obs)
    max_modelinterp = annualmax(temporalsubset(simp,(2000,01,01),(2099,01,01)))
    max_modelqqmap = annualmax(temporalsubset(out,(2000,01,01),(2099,01,01)))
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
    cm.plot(max_obs, label="OBS")
end
r = Raster("pre_cor.nc";key= "pre")
xd = ncdf(r)
dubplot(xd) #no dubs!!!!
dfm(xd;fun=monsum,ann=false,mode=:bar) #mit obs vergleichen!
hydromon(xd)
bargroup(xd;fun=monmean)
or = Raster("pre_rcm_obs.nc";key= "pre")
xdo = ncdf(or)
xdo.date .= xdo.date .+ Dates.Hour(12)
map(x->x.pre .= Float64.(x.pre),[xd,xdo])
mx = mall(xd,xdo)
bargroup(mx;fun=monmean)
dpr(mx)
@doc kge2
@doc Main.wa.bargroup
dfm(mx)



########loop for tas rcms ##########################
# cdb()
# pw()
# "/mnt/d/remo/cordex/eobs/v28"
ras = glob("rcm")
myvar = "tas"
mkdir(myvar)
cd(myvar)
oldname="tg";newname=myvar
run(`cdo -v -chname,$oldname,$newname ../tg_ens_mean_v28.0e_rcm.nc tas_obs.nc`)

function ncremap(input_file, remapfile, output_file, startdate, enddate;
    oldname::Union{Nothing, String}=nothing, newname::Union{Nothing, String}=nothing)
    try
    # Construct the CDO command for cropping
    # -setcalendar,365_day 
    # -P <nthreads>  Set the number of OpenMP threads
    cdo_command = `cdo -P 8 -seldate,$startdate,$enddate -remapbil,$remapfile`

    if oldname !== nothing
    cdo_command = `$(cdo_command) -chname,$oldname,$newname`
    end

    cdo_command = `$(cdo_command) $input_file $output_file`

    # Run the CDO command using run
    run(cdo_command)

    println("NetCDF file successfully processed. Output saved to $output_file")

    catch e
    println("Error: $e")
    println("Failed to remap $output_file.")
    end
end

begin
    simlnk="/mnt/d/remo/cordex/wgs/proj_tas_hist+rcp85_wgs84.nc"
    #output_file = myvar*"/simh.nc"
    output_file = "simh.nc"
    remapfile = "tas_obs.nc"
    ncremap(simlnk,remapfile,output_file,
    "1970-01-01","2020-12-31") #;oldname="tg",newname=myvar)
    #simpfile = myvar*"/simp.nc"
    simpfile = "simp.nc"
    ncremap(simlnk,remapfile,simpfile,
    "2021-01-01","2099-12-31") #;oldname="tg",newname=myvar)
end


simh = tt.load("simh.nc",myvar)
simp = tt.load("simp.nc",myvar)
obs  = tt.load(remapfile,myvar)

# qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; 
# method="Additive", detrend=true, window::Int=15, 
# rankn::Int=50, thresnan::Float64=0.1, 
# keep_original::Bool=false, 
# interp::Function = Linear(), 
# extrap::Function = Flat())

out = qqmap(obs,simh,simp;method="Additive",
    thresnan=0.1)

@doc     cm.contourf
cm.contourf(out;center_cs=true,cm="magma",titlestr="pre rsds")
#myvar*"_cor.nc"|>println
tt.write(out,myvar*"_cor.nc")
max_obs = annualmax(obs)
max_modelinterp = annualmax(temporalsubset(simp,(2000,01,01),(2099,01,01)))
max_modelqqmap = annualmax(temporalsubset(out,(2000,01,01),(2099,01,01)))
begin
    max_obs = annualmax(obs)
    
    cm.plot(max_obs, label="OBS")
    cm.plot(annualmax(simh), label="REMO - futre scenario RCP85")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
end
max_modelqqmap = annualmax(temporalsubset(out,(1970,01,01),(2050,01,01)))
cm.plot(max_modelqqmap)

begin
    #max_obs = annualmax(temporalsubset(obs,(2000,01,01),(2000,01,01)))
    cm.plot(max_obs, label="OBS")
    cm.plot(max_modelinterp, label="REMO - futre scenario RCP85")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
end


# fg            Ensemble mean wind speed [m/s]
# cdb()
# pw()
"/mnt/d/remo/cordex/eobs/v28"|>cd
ras = glob("rcm")
myvar = "sfcWind"
mkdir(myvar)
cd(myvar)
oldname="fg";newname=myvar
run(`cdo -v -chname,$oldname,$newname ../fg_ens_mean_v28.0e_rcm.nc $myvar-obs.nc`)

begin
    simlnk="/mnt/d/remo/cordex/wgs/proj_sfcWind_hist+rcp85_wgs84.nc"
    output_file = "simh.nc"
    remapfile = "$myvar-obs.nc"
    ncremap(simlnk,remapfile,output_file,
    "1970-01-01","2020-12-31") #;oldname="tg",newname=myvar)
    #simpfile = myvar*"/simp.nc"
    simpfile = "simp.nc"
    ncremap(simlnk,remapfile,simpfile,
    "2021-01-01","2099-12-31") #;oldname="tg",newname=myvar)
end

begin
    simh = tt.load("simh.nc",myvar)
    simp = tt.load("simp.nc",myvar)
    obs  = tt.load(remapfile,myvar)  
end

# qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; 
# method="Additive", detrend=true, window::Int=15, 
# rankn::Int=50, thresnan::Float64=0.1, 
# keep_original::Bool=false, 
# interp::Function = Linear(), 
# extrap::Function = Flat())

out = qqmap(obs,simh,simp;method="Additive",
    thresnan=0.001)

@doc     cm.contourf
cm.contourf(out;center_cs=true,cm="magma",titlestr="$myvar mean")
#myvar*"_cor.nc"|>println
tt.write(out,myvar*"_cor.nc")
max_obs = annualmax(obs)
max_modelinterp = annualmax(temporalsubset(simp,(2000,01,01),(2099,01,01)))
max_modelqqmap = annualmax(temporalsubset(out,(2000,01,01),(2099,01,01)))
begin   
    cm.plot(max_obs, label="OBS")
    cm.plot(annualmax(simh), label="REMO - futre scenario RCP85")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
end
max_modelqqmap = annualmax(temporalsubset(out,(1970,01,01),(2050,01,01)))
cm.plot(max_modelqqmap)

begin
    #max_obs = annualmax(temporalsubset(obs,(2000,01,01),(2000,01,01)))
    cm.plot(max_obs, label="OBS")
    cm.plot(max_modelinterp, label="REMO - futre scenario RCP85")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
end


"/mnt/d/remo/cordex/eobs/v28"|>cd
ssup()
ras = glob("rcm")
ofile=ras[2]
run(`cdo vardes $ofile`)
#mean relative humidity [%] => rh relative humidity [1]       
myvar = "rh"
mkdir(myvar)
cd(myvar)
oldname="hu";newname=myvar

lk="/mnt/d/remo/cordex/rh_hist+rcp85_19500102-21001231.nc"
#rhw = tt.load(lk,"rh")
run(`cdo vardes $lk`)

cm.plot(annualmax(simh))
# +to percent to 1/1
run(`cdo -v -chname,$oldname,$newname -mulc,0.01 ../$ofile $myvar-obs.nc`)
#cm.plot(annualmax(tt.load("$myvar-obs.nc",myvar)))
begin
    #simlnk="/mnt/d/remo/cordex/wgs/proj_rh_hist+rcp85_wgs84.nc"
    simlnk=lk
    output_file = "simh.nc"
    remapfile = "$myvar-obs.nc"
    ncremapcon(simlnk,remapfile,output_file,
    "1970-01-01","2020-12-31") #;oldname="tg",newname=myvar)
    #simpfile = myvar*"/simp.nc"
    simpfile = "simp.nc"
    ncremapcon(simlnk,remapfile,simpfile,
    "2021-01-01","2099-12-31") #;oldname="tg",newname=myvar)
end

begin
    simh = tt.load("simh.nc",myvar)
    simp = tt.load("simp.nc",myvar)
    #obs  = tt.load(remapfile,myvar)  
end

# qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; 
# method="Additive", detrend=true, window::Int=15, 
# rankn::Int=50, thresnan::Float64=0.1, 
# keep_original::Bool=false, 
# interp::Function = Linear(), 
# extrap::Function = Flat())
#Multiplicative is usually bounded variables such as precipitation and humidity.
out = qqmap(obs,simh,simp;method="Multiplicative",
    thresnan=0.01)

cm.plot(annualmax(simp))
cm.plot(annualmax(obs))
cm.plot(annualmax(out),label="qqmap")

cm.contourf(out;center_cs=true,cm="magma",titlestr="$myvar mean")
cm.contourf(simh;cm="RdYlBu_r",titlestr="$myvar mean")
cm.contourf(simp;cm="cmo.thermal",titlestr="$myvar mean")
#myvar*"_cor.nc"|>println
tt.write(out,myvar*"_cor.nc")
max_obs = annualmax(obs)
max_modelinterp = annualmax(temporalsubset(simp,(2000,01,01),(2099,01,01)))
max_modelqqmap = annualmax(temporalsubset(out,(2000,01,01),(2099,01,01)))
begin   
    cm.plot(max_obs, label="OBS")
    cm.plot(annualmax(simh), label="REMO - futre scenario RCP85")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
end




"/mnt/d/remo/cordex/eobs/v28"|>cd
ras = glob("rcm")
ofile=ras[end]
run(`cdo vardes $ofile`)
myvar = "tg"
shapefile_path="/mnt/d/Wasim/regio/rcm200/v4/ezg.shp"
plot_mean_with_shapefile(ofile, myvar, shapefile_path;titlestr=myvar)

cdof(shapefile_path)
glob(r"ezg")|>println


mkdir(myvar)
cd(myvar)
oldname="hu";newname=myvar

lk="/mnt/d/remo/cordex/rh_hist+rcp85_19500102-21001231.nc"
#rhw = tt.load(lk,"rh")
run(`cdo vardes $lk`)
function loadcdo()
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    @time_imports include(pt)
end
loadcdo()
cm.plot(annualmax(simh))
# +to percent to 1/1
run(`cdo -v -chname,$oldname,$newname -mulc,0.01 ../$ofile $myvar-obs.nc`)
#cm.plot(annualmax(tt.load("$myvar-obs.nc",myvar)))
begin
    #simlnk="/mnt/d/remo/cordex/wgs/proj_rh_hist+rcp85_wgs84.nc"
    simlnk=lk
    output_file = "simh.nc"
    remapfile = "$myvar-obs.nc"
    ncremapdis(simlnk,remapfile,output_file,
    "1970-01-01","2020-12-31") #;oldname="tg",newname=myvar)
    #simpfile = myvar*"/simp.nc"
    simpfile = "simp.nc"
    ncremapdis(simlnk,remapfile,simpfile,
    "2021-01-01","2099-12-31") #;oldname="tg",newname=myvar)
end

begin
    simh = tt.load("simh.nc",myvar)
    simp = tt.load("simp.nc",myvar)
    obs  = tt.load(remapfile,myvar)  
end

# qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; 
# method="Additive", detrend=true, window::Int=15, 
# rankn::Int=50, thresnan::Float64=0.1, 
# keep_original::Bool=false, 
# interp::Function = Linear(), 
# extrap::Function = Flat())
#Multiplicative is usually bounded variables such as precipitation and humidity.
out = qqmap(obs,simh,simp;method="Multiplicative",
    thresnan=0.01)

cm.plot(annualmax(simp))
cm.plot(annualmax(obs))
cm.plot(annualmax(out),label="qqmap")

cm.contourf(out;center_cs=true,cm="magma",titlestr="$myvar mean")
cm.contourf(simh;cm="RdYlBu_r",titlestr="$myvar mean")
cm.contourf(simp;cm="cmo.thermal",titlestr="$myvar mean")
#myvar*"_cor.nc"|>println
tt.write(out,myvar*"_cor.nc") #NaN
max_obs = annualmax(obs)
max_modelinterp = annualmax(temporalsubset(simp,(2000,01,01),(2099,01,01)))
max_modelqqmap = annualmax(temporalsubset(out,(2000,01,01),(2099,01,01)))
begin   
    cm.plot(max_obs, label="OBS")
    cm.plot(annualmax(simh), label="REMO - futre scenario RCP85")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
end



#try in bash
#cdo    genbic: Extrapolation disabled!
cdo griddes rh-obs.nc 
cdo griddes rh-obs.nc > gridfile.txt
lk="/mnt/d/remo/cordex/rh_hist+rcp85_19500102-21001231.nc"
cdo -v genbic,gridfile.txt $lk rwgen.nc
cdo -v -b F32 -remap,gridfile.txt,rwgen.nc $lk sim.nc
xrtsc sim.nc
###julia
myvar
simlnk="sim.nc"
sim=tt.load(simlnk,myvar)
simh=tt.temporalsubset(sim,(1970,01,01),(2020,12,31))
simp=tt.temporalsubset(sim,(2021,01,01),(2099,12,31))
obs=tt.load("rh-obs.nc",myvar)

out = qqmap(obs,simh,simp;
        method="Multiplicative",
        rankn=10,
        thresnan=0.000001)

cm.plot(annualmax(out),label="qqmap")


cm.plot(annualmax(simh))
cm.plot(annualmax(obs))


@doc cm.wbgt
@doc tt.vaporpressure

q=obs   #specific humidity

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
import Plots
Plots.plot(rs)

#options=Dict("COMPRESS"=>"DEFLATE")
Rasters.name(rs) = "orog"
@doc Rasters.name

Rasters.write("dhm.nc",rs;driver="netcdf",
    #options=Dict("NETCDF_VARNAME"=>"orog")
    options=Dict("variable_name"=>"orog"),force=true    )
run(`cdo vardes dhm.nc`)
##err generic cords.
# run(`cdo -chname,unnamed,orog -remapbil,gridfile.txt dhm.nc dhmR.nc`)
# run(`cdo vardes dhm2.nc`)
# or=tt.load("dhm2.nc","orog")        #orography (orog) (m)
# import NCDatasets
# const nx=NCDatasets
# ds = nx.Dataset("dhm.nc")

pt="/mnt/d/remo/cordex/eobs/elev_crop.nc"
run(`cdo vardes $pt`)
#run(`cdo -v -setdate,2000-01-01 -remapbil,gridfile.txt $pt dhm2.nc`)
#run(`cdo -v -settaxis,2000-01-01 -remapbil,gridfile.txt $pt dhm2.nc`)
run(`cdo -v -setdate,2000-01-01, -settime,12:00:00 -setcalendar,365_day -remapbil,gridfile.txt $pt dhm2.nc`)
run(`cdo showdate dhm2.nc`)
@doc tt.load
or=tt.load("dhm2.nc","elevation";data_units="m")        #orography (orog) (m)

@doc tt.vaporpressure
temp=tt.load("../tas/tas_obs.nc","tas")
pp=tt.load("../../pp_ens_spread_crop.nc","pp")
q=tt.load("rh-obs.nc","rh")
vp = tt.vaporpressure(q,pp,or,temp)


using GeoStats
using ClimateTools

target = :rh
n = 30 # max number of neighboring points
solver = Kriging(target => (maxneighbors=n,))
C = regrid(q, obs, solver=solver)

##so regrid == griddata !?

# Interpolate `ClimGrid` A onto the lon-lat grid of `ClimGrid` B,
# where A and B are `ClimGrid`.
#function griddata(A::ClimGrid, B::ClimGrid; method="linear", 
    #min=[], max=[])
#tt.griddata(q,lon=3,lat=4)
@doc tt.griddata
C = tt.griddata(pp,obs)        #5min!
obs|>cm.contourf
lk="/mnt/d/remo/cordex/rh_hist+rcp85_19500102-21001231.nc"
#sim=tt.load(lk,"rh";data_units="%")
sim=tt.load("simh.nc","rh";data_units="%")
sim|>cm.contourf
C = tt.griddata(sim,obs)        #5min!

a,b=ClimateTools.getgrids(obs)
map(extrema,[a,b])
@doc get_timevec
t = get_timevec(obs)
#unique.(Dates.dayofyear.(t))
obs.globalattribs
out.globalattribs|>DataFrame
d = out.globalattribs
DataFrame(d)

obs.globalattribs["history"]
simp.globalattribs["history"]

cm.plot(annualmax(C))
cm.contourf(C)

#ClimateTools.jl/src/biascorrect.jl
#include(joinpath(src_path, "biascorrect_git.jl"))

@doc tt.biascorrect_extremes

simlnk="sim.nc"
sim=tt.load(simlnk,myvar)
simh=tt.temporalsubset(sim,(1970,01,01),(2020,12,31))
simp=tt.temporalsubset(sim,(2021,01,01),(2099,12,31))
#cdo -v chunit,%,1 rh-obs.nc rh-obs1.nc
obs=tt.load("rh-obs1.nc",myvar)
#  Combine the empirial Quantile-Quantile mapping (see qqmap) and Generalized Pareto Distribution
#bias-correction methods
exout=tt.biascorrect_extremes(obs,simh,simp;)
#exout=tt.biascorrect_extremes(sim,simh,simp;)
out = qqmap(obs,simh,simp;method="Multiplicative")
out.globalattribs["history"]
cm.contourf(out)
cm.plot(annualmax(out),label="qqmap")

out.globalattribs|>grep(r"name"i)
out.varattribs
obs.varattribs

#pycall test
using Conda
#Conda.add("")
Conda.pip_interop(true;)
#pip install python-cmethods
Conda.pip("install","python-cmethods")
#import xarray as xr
##https://pypi.org/project/python-cmethods/
#from cmethods import CMethods as cm
@pyimport cmethods

obsh = xr.open_dataset("rh-obs1.nc")
simh = xr.open_dataset("simh.nc")
simp = xr.open_dataset("simp.nc")
myvar = "rh"
x0 = obsh[myvar].to_array()
x0 = obsh[myvar][:,0,0]

x0 = py"$obsh[$myvar][:,0,0]"
x1 = py"$simh[$myvar][:,0,0]"
x2 = py"$simp[$myvar][:,0,0]"

ls_result = cmethods.CMethods.linear_scaling(
    obs = x0,
    simh = x1,
    simp = x2,
    kind = "*" #*> not available for linear_scaling.
)

#ti = get_timevec(sim)
cm.plot(ls_result)

#ValueError("Dimensions {'lat', 'lon'} do not exist. 
#Expected one or more of ('time', 'latitude', 'longitude')")
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
    obs = obsh[myvar],
    simh = simh[myvar],
    simp = simp[myvar],
    n_quaniles = 1000,
    kind = "*"
)
#[myvar]
qdm_result.mean("time").plot()
#ValueError("Dimensions {'lat', 'lon'} do not exist. Expected one or more of ('time', 'latitude', 'longitude')
qdm_result.to_netcdf("rh_cor.nc")

out=tt.load("rh_cor.nc",myvar)
out.varattribs["comment"] = "bias corrected via cmethods"
out.varattribs|>grep("bia")

myvar="rh"
shapefile_path="/mnt/d/Wasim/regio/rcm200/v4/ezg.shp"
println(basename(shapefile_path)," -> ",filesize(shapefile_path) ./ 1024^2, " MB")

plot_mean_with_shapefile("rh_cor.nc", myvar, 
    shapefile_path)
cm.contourf(out;titlestr="obs rh")
cm.plot(annualmax(out),label="qqmap")

simh=tt.temporalsubset(sim,(1970,01,01),(2020,12,31))
simp=tt.temporalsubset(sim,(2021,01,01),(2099,12,31))

max_obs = annualmax(obs)    
max_modelinterp = annualmax(simp)
max_modelqqmap = annualmax(out)
begin   
    cm.plot(max_obs, label="OBS")
    cm.plot(annualmax(simh), label="REMO - futre scenario RCP85")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values of $myvar")
    #filename = "tso.png")
end

cd("/mnt/d/remo/cordex/eobs/v28")
mkdir("biascorr")
v = rglob("_cor")
@doc cp
map(x->cp(x,joinpath("biascorr",basename(x))),v)
cd("biascorr")
#bash
xbash cor vardes
infile="rh_cor.nc"
#ncatted -a bounds,time,d,d $infile outfile.nc
#cdo showtimestamp $infile
##add time_bnds #https://nco.sourceforge.net/nco.html#make_bounds
ncap2 -O -s 'defdim("bnds",2);
time_bnds=make_bounds(time,$bnds,"time_bnds");' $infile rh-t.nc
vardes rh-t.nc #no warning.
pydesc rh-t.nc
griddes rh-t.nc

#make lists.
cp -v /mnt/d/remo/cordex/eobs/mw_eobs_fullwinlist.sh mw_eobs_corlists.sh
a=pre_cor.nc
cdo showtimestamp  $a > ts
sed 's/ /\n/g' ts|sed 's/-\|T\|:/\t/g'|sed '/^$/d'|sed -r 's/\s+\S+$//' > ts2
#ht
#perl -e 'while (my $line = <>) {if ( $line =~ $.==1|eof) {print $line; }}'
#cdo showtimestamp pre_cor.nc |head -2 > z
#perl -e 'while (my $line = <>) {if ( $line =~ $.==1) {print $line; }}' z
#check timestamps.

# for i in *nc;do (cdo showtimestamp $i|tr ' ' '\n'|head -n20|tail +1);done

#    bash mw_eobs_corlists.sh 

##Internal message is: NetCDF: Start+count exceeds dimension bound
##reorder vars
ncks -x -v Regular_longitude_latitude -v lat -v lon  ../tas/tas_cor.nc ta.nc
#cdo -v -f nc reorder,lon,lat,time ta.nc $a
vardes ta.nc
ncpdq -a lon,lat ta.nc tmp.nc

ncpdq --rdr=lon,lat,tas ta.nc tmp.nc
vardes tmp.nc
a=tmp.nc
ncks -A -v lon ta.nc $a
ncks -A -v lat ta.nc $a
ncks -A -v time ta.nc $a
vardes $a

@doc reorder

A = tt.load("tmp.nc","tas")
cm.contourf(A)



ncks -x -v Regular_longitude_latitude,lat,lon  ../tas/tas_cor.nc tas_cor.nc
vardes tas_cor.nc

da = xr.open_dataset("tas_cor.nc")
da.dims
dx = da.transpose("time","longitude","latitude")
dx = dx.drop("Regular_longitude_latitude")
dx.to_netcdf("tas_cor2.nc")
mv tas_cor2.nc tas_cor.nc
vardes tas_cor.nc
#['latitude', 'longitude', 'tas', 'time'] 
ncdump -h tas_cor.nc

tas_cor = xr.open_dataset("../tas/tas_cor.nc")

# Create a new NetCDF dataset
tasre = xr.Dataset(
    coords=Dict(
        "x"=>xr.DataArray(tas_cor.longitude, name="x"),
        "y"=>xr.DataArray(tas_cor.latitude,  name="y"),
        "t"=>xr.DataArray(tas_cor.time, name="t")
        ),
        data_vars=Dict(
            "tas"=>xr.DataArray(tas_cor.tas, name="tas", 
            #dims=("t","x","y"))
            dims=("time","latitude","longitude"))
        ),
        attrs=Dict("history"=>"preprocessed with Julia and PyCall xarray")        
        )


        # attrs=Dict()
        #     , units="m East of reference point"),
        #     , units="m North of reference point"),
        # )

tasre.to_netcdf("tas_cor.nc")
tas_cor.close()
tasre.close()
#mv tas_cor3.nc tas_cor.nc



ncks -x -v Regular_longitude_latitude,lat,lon  ../tas/tas_cor.nc tc.nc
#-d, --dmn, dimension old_dim,new_dim Dimension's old and new names
ncrename -O -d latitude,y -d longitude,x -d time,t tc.nc
ncdump -h tc.nc

ncks -x -v Regular_longitude_latitude,latitude,longitude tc.nc t2.nc

#ncpdq --rdr=lon,lat,tas tc.nc tmp.nc
ncpdq --rdr=t,x,y t2.nc tas_cor.nc
ncdump -h tas_cor.nc
ncatted -O -a tas:missing_value,,o,f,-9999. tas_cor.nc

a=genRE_precipitation_hour_1.1.nc
#ncpdq --rdr=x,y,pr $a tmp.nc #killed
#-x, --xcl, exclude      Extract all variables EXCEPT those specified with -v
#--xcl -v=x,y,pr -d
#ncra -d time,1,24,24 $a tmp.nc #killed
#ncks --mk_rec_dmn time $a interim.nc #killed
ncks --mk_rec_dmn time $a tmp.nc
ncra -y ttl -d time,1,,24 $a tmp.nc


#->works

# ncrename -O -v "longitude","x" tc.nc
# ncrename -O -v "latitude","y" tc.nc
# ncatted -O -a missing_value,,o,f,-9999. tc.nc
# ncatted -O -a _FillValue,,o,f,-9999. tc.nc
# #float tas(t, x, y) ;
# ncdump -h tc.nc


ncks -x -v Regular_longitude_latitude,lat,lon  pre_cor.nc tc.nc
#-d, --dmn, dimension old_dim,new_dim Dimension's old and new names
ncrename -O -d latitude,y -d longitude,x -d time,t tc.nc
ncdump -h tc.nc
ncks -x -v Regular_longitude_latitude tc.nc t2.nc
ncpdq -x -v Regular_longitude_latitude --rdr=t,x,y t2.nc pre_cor.nc
ncdump -h pre_cor.nc
ncatted -O -a pre:missing_value,,o,f,-9999. tas_cor.nc

a=rsds_cor.nc 
ncks -O -x -v Regular_longitude_latitude,lat,lon  $a tc.nc
ncrename -O -d latitude,y -d longitude,x -d time,t tc.nc
ncdump -h tc.nc
ncpdq -x -v Regular_longitude_latitude --rdr=t,x,y tc.nc $a
ncdump -h $a
ncatted -O -a missing_value,,o,f,-9999. $a

a=rh_cor.nc
ncks -O -x -v Regular_longitude_latitude,lat,lon  $a tc.nc
ncrename -O -d lat,y -d lon,x -d time,t tc.nc
ncdump -h tc.nc
ncpdq -x -v lon --rdr=t,x,y tc.nc $a
ncpdq -x -v lat --rdr=t,x,y $a tmp.nc
mv -v tmp.nc $a
ncatted -O -a missing_value,,o,f,-9999. $a
ncdump -h $a

a=sfcWind_cor.nc
ncks -O -x -v Regular_longitude_latitude,lat,lon  $a tc.nc
ncrename -O -d latitude,y -d longitude,x -d time,t tc.nc
ncdump -h tc.nc
ncpdq -x -v Regular_longitude_latitude -v latitude -v longitude --rdr=t,x,y tc.nc $a
ncatted -O -a missing_value,,o,f,-9999. $a
ncdump -h $a


# can you rewrite this bash cmds using PyCall and xarray? 
# a=sfcWind_cor.nc
# ncks -O -x -v Regular_longitude_latitude,lat,lon  $a tc.nc
# ncrename -O -d latitude,y -d longitude,x -d time,t tc.nc
# ncpdq -x -v Regular_longitude_latitude -v latitude -v longitude --rdr=t,x,y tc.nc $a
# ncatted -O -a missing_value,,o,f,-9999. $a
# ncdump -h $a

#restyle eobs
cd("/mnt/d/remo/cordex/eobs/v28")
ls()
ras = glob("rcm")
ras = ras[Not(5)]
using PyCall
#xr = pyimport("xarray")
subprocess = pyimport("subprocess")
# File name
#a = "fg_ens_mean_v28.0e_rcm.nc" #wrong timeaxis.
cdo -v remapbil,rcm_eobs/hu_wa.nc.gridfile fg_ens_mean_0.1deg_reg_v28.0e.nc fg_rcm.nc
# so 1980 is the first year.
a = "fg_rcm.nc"
# Open the netCDF file
ds = xr.open_dataset(a)
# # Drop the specified variables
# ds = ds.drop_vars(["Regular_longitude_latitude", "lat", "lon"])
# # Save to a new netCDF file
# ds.to_netcdf("tc.nc")
# Rename dimensions
ds = ds.rename_dims(Dict("latitude" => "y", "longitude" => "x", "time" => "t"))
#ds.rename_vars
# Reorder dimensions
ds = ds.transpose("t", "x", "y")
# Rename Coordinates <- no. errors. time must hold.
# ds = ds.rename(Dict("latitude" => "y", "longitude" => "x", "time" => "t"))
# Update attribute
ds.attrs["missing_value"] = -9999.
# Shift time by 12 hours
#ds['time'] = ds.indexes['time'].shift(12, 'H')
# Save changes to the original netCDF file
#ds.to_netcdf(a) #'Permission denied')
mkdir("rcm_eobs")
#outfile = replace(a, "_ens_mean_v28.0e_rcm.nc" => "_wa.nc")
outfile = ("fg_wa.nc")
outpt = joinpath("rcm_eobs",outfile)
ds.to_netcdf(outpt)
# Use subprocess to call ncdump -h
subprocess.call(["ncdump", "-h", outpt])
ds["fg"].mean("t").plot()
ds["fg"].mean("x").mean("y").plot()
ds["fg"].max("x").max("y").plot()

#k=tt.load(outpt,"fg")
# # so i have to make a own winlist for fg, starting from 1980.
# cd rcm_eobs
# cdo showtimestamp  $a | sed 's/ /\n/g' |sed 's/-\|T\|:/\t/g'|sed '/^$/d'|sed -r 's/\s+\S+$//' | perl -lne '@F=split(/\t/,$_);{$F[3]=~s/00/12/g;};print join("\t",@F)' > fgt
# var="fg_wa.nc"
# cdo griddes $var > ${var}.gridfile
# onam="${var/.nc/.winlist}"
# touch "$onam"
# echo "GRID LIST of" $var > $onam
# printf 'YY\tMM\tDD\tHH\t9999\n'  >> $onam   #line 2
# printf 'YY\tMM\tDD\tHH\t'  >> $onam 
# awk '/xfirst/ {print substr($0, RSTART +12, RLENGTH + 6);}' ${var}.gridfile >> $onam 
# printf 'YY\tMM\tDD\tHH\t'  >> $onam 
# awk '/yfirst/ {print substr($0, RSTART +12, RLENGTH + 6);}' ${var}.gridfile >> $onam 
# printf 'YY\tMM\tDD\tHH\t'Liste'\n'  >> $onam 
# awk '{OFS=FS}{$NF=FS"D:/remo/cordex/eobs/v28/rcm_eobs/outw.nc<param>"i-1+1;i++;print}' fgt >> $onam
# param=${var//_wa.nc/}	
# perl -i -pe 's/outw.nc/'$var'/g' "$onam"
# perl -i -pe 's/param/'$param'/g' "$onam"
# sed -i '1s,Liste,'"$param"',' "$onam"       #not glob, line 1
# perl -i -pe 's/\h+/\t/g' "$onam"


for nc in ras[2:end]
    ds = xr.open_dataset(nc)
    # Rename dimensions
    ds = ds.rename_dims(Dict("latitude" => "y", "longitude" => "x", "time" => "t"))
    #ds.rename_vars
    # Reorder dimensions
    ds = ds.transpose("t", "x", "y")
    # Rename Coordinates
    #ds = ds.rename(Dict("latitude" => "y", "longitude" => "x", "time" => "t"))
    ds.attrs["missing_value"] = -9999.
    outfile = replace(nc, "_ens_mean_v28.0e_rcm.nc" => "_wa.nc")
    outpt = joinpath("rcm_eobs",outfile)
    ds.to_netcdf(outpt)
    println("$nc done... \n saved to $outpt")
    ds.close()
end



#xbash nc vardes
"/mnt/d/remo/cordex/eobs/v28/rcm_eobs"|>cd
cp("../biascorr/mw_eobs_corlists.sh","mw_eobs_obslists.sh")
mw_eobs_obslists.sh


# @pyimport import xarray as xr
@pyimport import matplotlib.pyplot as plt
@pyimport cartopy.crs as ccrs

da = ds.isel(t=0)
begin
    ax = plt.axes(projection=ccrs.Orthographic(10, 49))
    ax.set_global()
    plt.contourf(da,ax=ax, transform=ccrs.PlateCarree())
    da.coastlines()
end

@pj
"/mnt/d/Wasim/regio/out/rc200/x22/f18-loc"|>cd
fs = @nco "sb0"
pyjl.xrp(fs)
pyjl.xrfacets(fs)

"/mnt/d/Wasim/regio/out/rc200/x30/c2/"|>cd
fs = @nco "sb0"
pyjl.xrp(fs)
pyjl.pyplot_df(r"sb"|>dfr;log=true)
pyjl.pyplot_df(r"qges"|>dfr;log=true)





########### qmap test with full time
pt="/mnt/d/remo/qm/corgrids/pre"

function loadcdo()
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end

loadcdo()
using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)
cd(pt)
begin
    simlnk="/mnt/d/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
    output_file = "simh.nc"
    remapfile = "/mnt/d/remo/cordex/eobs/v28/pre/pre_rcm_obs.nc"
    ncremap(simlnk,remapfile,output_file,"1990-01-01","2020-12-31")
    simpfile = "simp.nc"
    ncremap(simlnk,remapfile,simpfile,"1990-01-01","2099-12-31")
end


simh = tt.load("simh.nc","pre")
simp = tt.load("simp.nc","pre")
obs  = tt.load(remapfile,"pre")

# qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; 
# method="Additive", detrend=true, window::Int=15, 
# rankn::Int=50, thresnan::Float64=0.1, 
# keep_original::Bool=false, 
# interp::Function = Linear(), 
# extrap::Function = Flat())

#31.995314 seconds 
@time out = qqmap(obs,simh,simp;method="Multiplicative",thresnan=0.001)

cm.contourf(out;
    cm="magma",titlestr="pre rsds")

tt.write(out,"pre_cor.nc")

begin
    #max_obs = annualmax(temporalsubset(obs,(2000,01,01),(2000,01,01)))
    max_obs = annualmax(obs)
    max_modelinterp = annualmax(temporalsubset(simp,(1990,01,01),(2099,01,01)))
    max_modelqqmap = annualmax(temporalsubset(out,(1990,01,01),(2099,01,01)))
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")
    #filename = "tso.png")
    cm.plot(max_obs, label="OBS")
end

#passt soweit. muss ich nun im code umschreiben.
#evtl mit regrid versuchen.

#switch to win
pwin()

# rpt = "D:/remo/qm/corgrids/pre/pre_cor.nc"
# r = Raster(rpt;key= "pre")
# xd = ncdf(r)
# dubplot(xd) #no dubs!!!!

pt="/mnt/d/remo/qm/corgrids/tas"
cd(pt)
readdir() 

##tt.regrid test. ==> griddata

cvar = Symbol("tas")
obspt = "/mnt/d/remo/cordex/eobs/v28/tas/tas_obs.nc"
temp=tt.load(obspt,string(cvar))
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

cm.plot(annualmean(temp))

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



A=tt.load("/mnt/d/remo/qm/corgrids/pre-cor.nc","pre")
B=tt.load("/mnt/d/remo/qm/corgrids/pre/pre-cor.nc","pre")
obs = tt.load("/mnt/d/remo/cordex/eobs/v28/pre/pre_rcm_obs.nc","pre")

cm.contourf(A)
cm.contourf(B)

a = annualsum(temporalsubset(A,(2000,01,01),(2022,12,31)))
b = annualsum(temporalsubset(B,(2000,01,01),(2022,12,31)))
c = annualsum(temporalsubset(obs,(2000,01,01),(2022,12,31)))
cm.plot(a, label="full - old")
cm.plot(b, label="jl - interpolated")
cm.plot(c, label="eobs")

plt.rcParams["font.family"] = "cmr10"
plt.rcParams["font.size"] = 22
plt.rcParams|>grep("font")
plt.rcParams|>grep("font.size")


d = plt.rcParams
k,v = string.(keys(d)),string.(values(d))
#DataFrame(hcat(k,v),:auto)
df = DataFrame(nam=k,val=v)
filter(row -> occursin("font",row.nam),df)



nd = filter(kv -> !isempty(first(kv)), d)
d==nd
#|>DataFrame
DataFrame(d,:auto)
#filter(kv -> (keys(kv)|>first), d)

keys(d)|>first
values(d)|>first
(d)|>first


ssup()
pt="/mnt/d/Wasim/Testdata/Control/wasim_soil_props.txt"
nms = DataFrame(CSV.File(pt;header=1,limit=0))
#df = DataFrame(CSV.File(pt;header=false,skipto=2,delim="\t"))
#df = CSV.read(pt, DataFrame;header=false,skipto=2,delim="\t",transpose=true)
df = CSV.read(pt, DataFrame;transpose=true)
dropmissing(df,1)

platform != "windows"
dir="D:/Wasim/regio/control"
wsl_cmd = `wslpath $(dir)`
dir = readchomp(pipeline(wsl_cmd))


findctl("f9")


##radolan 
fn="/mnt/d/Wasim/regio/rcm/radolan/pre_radolan_2016.nc"
p = tt.load(fn,"pre")

p = xr.open_dataset(fn)
p = p.rename(Dict("easting"=>"x","northing"=>"y"))
p = p.rename(Dict("time"=>"t"))
p = p.drop_vars("time_bnds")
p = p.transpose("t","x","y","bnds")
p["pre"].mean("t").plot()

#@pyimport numpy as np
# Add 12 hours
#p["pre"]["t"] = p["pre"]["t"] + np.timedelta64(12, "h")
cdof(fn)
p.to_netcdf("pre_rl_2016.nc")
p.close()
#op()
fn="/mnt/d/Wasim/regio/rcm/radolan/pre_rl_2016.nc"
#p = tt.load(fn,"pre")
macro say_hello(surname)
    return :( println("Hello, ", $(esc(surname))) )
end
@say_hello("John")

pwd()
cd("/mnt/d/Wasim/regio/rcm/radolan")
rglob(".nc")


###corr genre
cd "D:\remo\genRE"
ncdump -h pre-REcor_raw.nc 
ncpdq
a=pre-REcor_raw.nc 
ncpdq --rdr=time,x,y $a pre-REcor.nc
ll
ncdump -h pre-REcor.nc
ncrename -O -d time,t pre-REcor.nc

a=pre-REcor.nc
cdo showtimestamp  $a > ts
cd("D:/remo/genRE/")
a="pre-REcor.nc"
open("ts2", "w") do out
    ts_content = String(read(pipeline(`wsl cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n") #s/ /\n/g
    #ts_content = replace(ts_content, r"\s+" => "\t")
    ts_content = replace(ts_content, r"-|T|:" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "\n" => "", count=1)
    #ts_content = replace(ts_content, r"\s+" => "\t")
    write(out, ts_content)
end


function process_files(;suffix=".winlist")
    files = filter(x -> occursin(r"pre-REcor+.*nc",x) && endswith(x, ".nc"), readdir())
    if any(x -> occursin("cor", x) && endswith(x, suffix), readdir())
        println("Corresponding files already exist. Exiting now!")
        return
    end
    for var in files
        gridfile = string(replace(var,".nc" => ".gridfile"))
        global onam = string(replace(var,".nc" => suffix))

        open(gridfile, "w") do out
            write(out, String(read(pipeline(`wsl cdo griddes $var`))))
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
process_files()
##fix for winlist
#cut -f1-5 tmp > pre-REcor.winlist
#mv -v $a /mnt/d/remo/qm/corgrids/jlcor/pre-REcor.nc



function add_24_after_23(input_file_path::AbstractString, output_file_path::AbstractString)
    lines = readlines(input_file_path)

    for i in 1:length(lines)-1
        if contains(lines[i], "\t23\t")
            parts = split(lines[i], '\t')
            #lines[i] = string(join(join(parts[1:3], '\t'),"24\t", parts[end],'\t'))
            lines[i] = string(join(parts[1:3], '\t')*"24\t"*parts[end])
            #lines[i] = string("24\t", join(parts[2:end], '\t'))
        end
    end

    write(output_file_path, join(lines, '\n'))
end


# Example usage:
input_file_path = "genRE_precipitation_hour_1.1.winlist"
output_file_path = "hrs.winlist"
add_24_after_23(input_file_path, output_file_path)
npp("hrs.winlist")
rm("hrs.winlist")


###the functions were malicious
#this is the fix.
cut -f5 $pt > unix-pt
pt=/mnt/d/remo/qm/corgrids/jlcor/pre-REcor.winlist 
cut -f5 $pt > win-pt
a=ts_edit
b=win-pt 
paste -d'\t' $a $b > tst
hdlat
paste -d'\t' <(cut -f1-4 $a) $b > tst
hdlat
paste -d'\t' <(cut -f1-4 $a) $b > pre-REcor.winlist
c=unix-pt 
paste -d'\t' <(cut -f1-4 $a) $c > pre-REcor.nc_list
hdlat

pt=/mnt/d/remo/qm/corgrids/jlcor/pre-REcor.winlist
cp -v pre-REcor.winlist $pt
pt=/mnt/d/remo/qm/corgrids/jlcor/pre-REcor.nc_list
cp -v pre-REcor.nc_list $pt



cd("D:/remo/genRE/")
a="pre-REcor.nc"