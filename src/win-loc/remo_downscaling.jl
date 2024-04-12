using Pkg
"D:/remo/qm/qm/"|>Pkg.activate
#Pkg.update()
using OhMyREPL
# pt="C:/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
# include(pt)
using ClimateTools
Pkg.status()
cd("D:/remo/qm/")
pwd()
fn="D:/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
cvar = "pre"
import ClimateTools
m = ClimateTools.load(fn,cvar)
simh = ClimateTools.temporalsubset(m,(1950,01,01),(2015,01,01))
simp = ClimateTools.temporalsubset(m,(1950,01,01),(2100,01,01))
obspt = "D:/remo/cordex/eobs/v28/pre/pre_rcm_obs.nc" 
obsras = ClimateTools.load(obspt,cvar)
# Interpolate ClimGrid A onto the lon - lat grid of ClimGrid B
obslow = ClimateTools.griddata(obsras,simh) #eta 3min
pqm = ClimateTools.qqmap(obslow,simh,simp)  #eta 10min
ClimateTools.write(pqm,"pre-cor-low.nc") 
#vgjl("extremes")
@doc ClimateTools.biascorrect_extremes
#takes forever....
pqex = ClimateTools.biascorrect_extremes(obslow,simh,simp;window=30)
ClimateTools.write(pqex,"pre-cor-low-extremes.nc") 

using ClimatePlots
const ct=ClimatePlots 
amx = annualmax(pqex)
ct.plot(amx)
asm = annualsum(pqex)
ct.plot(asm)


using PackageCompiler
#pt,nm=("C:/Users/chs72fw/.julia/sysimages","climatetools.so")      
pt,nm=("J:/jras","climatetools.so")      
create_sysimage([
    :ClimateTools,
    :DataFrames],sysimage_path=joinpath(pt,nm))
#✔ [07m:53s] PackageCompiler: compiling incremental system image
    #,    :CFTime,:Dates