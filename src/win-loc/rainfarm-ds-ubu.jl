cd("/mnt/e/qq/rainfarm")
prj="/home/ubu/.julia/rainfarm"
using Pkg
Pkg.activate(prj)
ENV["PROJ_LIB"]="/home/ubu/.julia/conda/3/share/proj/proj.db"
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


ssup()
@doc rglob
rglob("pr","/mnt/d/remo/cordex/eobs/v28/pre/")

sh="/mnt/d/remo/cordex/eobs/v28/pre/simh.nc"
sp="/mnt/d/remo/cordex/eobs/v28/pre/simp.nc"
ob="/mnt/d/remo/cordex/eobs/v28/pre/pre_rcm_obs.nc"
x_start, x_end, y_start, y_end = 9.348348,10.748348,49.309521,50.989521
pwd()
println("cropping all grids to 14x14...")
osh="simhcrop.nc"
run(`cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $sh $osh`)
osp="simpcrop.nc"
run(`cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $sp $osp`)
oobs="obscrop.nc"
run(`cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $ob $oobs`)
dem="/mnt/e/qq/fabdem.nc"
odem="fabdemcrop.nc" #5040x5688
run(`cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $dem $odem`)

simlnk="/mnt/d/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
osim = "proj_pre_crop.nc"
#remo bounds which gives 14x14 grid:
x_start, x_end, y_start, y_end = 9.278348,10.818348,51.05952,49.409521
run(`cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $simlnk $osim`)

##tt.griddata => PRECI
pwd()
#begin 
    cvar = Symbol("pre")
    obspt = oobs
    obsras = tt.load(obspt,string(cvar))
    #@doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    #"1990-01-01","2099-12-31"  
    cm.plot(annualmean(obsras))
#    simh = tt.load(osh,string(cvar))
#    simp = tt.load(osp,string(cvar))
    sim = tt.load(osim,string(cvar))
    simh = temporalsubset(sim,(1950,01,01),(2015,01,01))
    simp = temporalsubset(sim,(2015,01,01),(2100,01,01)) #changed to 2100
    all = tt.griddata(sim,obsras) 
    tt.write(all,"pre_remo_regridded_1950-2100.nc")

    simh = tt.griddata(simh,obsras) 
    simp = tt.griddata(simp,obsras) 
    contourf(simh)
    contourf(obsras)
    @doc qqmap
    # method::String="Additive", detrend::Bool=true, window::Int64=15, rankn::Int64=50, thresnan::Float64=0.1,
    out = qqmap(obsras,simh,simp,method="Multiplicative")
    tt.write(out,"pre_jlcor.nc")
    contourf(out)
    #Combine the empirial Quantile-Quantile mapping (see qqmap) and Generalized Pareto Distribution bias-correction methods.
    @doc tt.biascorrect_extremes
    outextr = tt.biascorrect_extremes(obsras,simh,simp)
    tt.write(outextr,"pre-crop_extremes.nc")
    contourf(outextr)
#now rainfarm
ls()
rglob("jl","$src_path/tools/")
rf=joinpath("$src_path/tools/","rfarm.jl")
#include(rf)
#in wsl
rf

pwc()
/mnt/e/qq/rainfarm
rf="/mnt/c/Users/Public/Documents/Python_Scripts/julia/tools/rfarm.jl"
rfw="/mnt/c/Users/Public/Documents/Python_Scripts/julia/tools/rfweights.jl"
prj="/home/ubu/.julia/rainfarm"
jrf="julia --threads auto -q --startup-file=no --project=$prj $rf"
$jrf -h
#julia --threads auto -q --startup-file=no --project=$prj $rfw fabdemcrop.nc pre_jlcor.nc
julia --threads auto -q --startup-file=no --project=$prj $rfw fabdemcrop.nc obscrop.nc

#ipython
import xarray as xr
ds = xr.open_dataset("pre-crop_extremes.nc")
ds = ds.drop("Regular_longitude_latitude")
ds.to_netcdf("pre-crop_extremes2.nc")
ds.close()

ipython
import xarray as xr
ds = xr.open_dataset("pre_jlcor.nc")
ds = ds.drop("Regular_longitude_latitude")
ds.to_netcdf("pre_jlcor2.nc")
ds.close()

ncatted -O -a _FillValue,,o,f,-9999 pre_jlcor2.nc

#rp weights.nc 
ref="pre_jlcor2.nc"
ref="pre-crop_extremes2.nc"
$jrf -w weights.nc $ref -o rfjl #breaks
$jrf -s 1.5 -n 2 $ref -o rfjl
cdo -v remapbil,../rainfarm_fine.gridfile pre_jlcor.nc pre_biljl.nc
$jrf -n 2 pre_biljl.nc -o rfjl

ref="obscrop.nc"
$jrf -w weights.nc $ref -o obsrfjl #breaks

#['pre', 'time', 'lon', 'lat'] => ['pre', 'time', 'x', 'y']
z=tt.load("obsrfjl_0001.nc","pre")
cm.contourf(z)

f = towsl("E:/qq/eobs/proj/pre-obs.nc")
z=tt.load(f,"pre")
cm.contourf(z)
z = xr.open_dataset(f)
"pre-crop_extremes.nc"

vgjl("PEST")

julia> pathof(RainFARM)
"/home/ubu/.julia/packages/RainFARM/tWWUJ/src/RainFARM.jl"  
