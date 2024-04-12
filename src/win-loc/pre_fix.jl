#RR check
#/mnt/c/Users/Public/Documents/Python_Scripts/julia/dscale-ubu.jl   
myvar="pre"
#daily precipitation sum RR, 
wrkdir="/mnt/d/remo/qm/prec"
cd(wrkdir)
ept = "/mnt/d/remo/cordex/eobs/"
using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)

ssup()
@gl "pr"
df = wread("Arnstein-Muedesheim_pre+rcp85.tsv";skip=5)
first(df)
plot(df.date,df.x5) #thats a pyplot
dy=yrsum(df)

obs=load("../prec/obsp.nc","pre")
contourf(obs)
pr=load("../prec/pre_delta_up2.nc","pre")
contourf(pr;cm="turbo")

cm.plot(annualsum(pr))
nconly("pre")

pcm=load("../prec/pre_cm.nc","pre")
cm.plot(annualsum(pcm))

#cm.mapclimgrid(pr,"Europe")
#using PyPlot
function mymap()
    central_lon = 10.0
    central_lat = 50.0
    extent = [0, 30, 40, 60]
    #proj_extent = cartopy.crs.PlateCarree()
    projection=cartopy.crs.LambertConformal(central_longitude=central_lon, central_latitude=central_lat)
    resolution = "auto"
    fig, ax = plt.subplots()
    
    ax = plt.axes(projection=projection)
    #ClimatePlots.contourf(pr)
    if ClimateTools.@isdefined extent
        if ClimateTools.@isdefined proj_extent
            ax.set_extent(extent, proj_extent)
        else
            ax.set_extent(extent)
        end
    end
        # Draw the map properties
        ax.gridlines()
        ax.coastlines(resolution=resolution)
    return true, ax, projection

end

mymap(pr)

# @pyimport xarray as xr
# @pyimport matplotlib.pyplot as plt
# ob = xr.open_dataset("obs_mean.nc")
# display(ob)
# x=ob.where(ob[myvar] > 0,drop=true)
# x.to_netcdf("obs_fix2.nc")

#same as cm.plot annualmean
x = xr.open_dataset("obs_mean.nc")
x[myvar].mean("longitude").mean("latitude").groupby("time.year").mean().plot()
x[myvar].isel(longitude=3,latitude=5).groupby("time.year").mean("time").plot()

groupby("time.year").mean("time").plot()


#cdo -v remapsetgrid,$croptg $obslnk obs.nc
simlnk="/mnt/d/remo/cordex/wgs/proj_rsds_hist+rcp85_wgs84.nc"
# cdo -v chname,qq,rsds -setcalendar,365_day -seldate,1970-01-01,2017-12-31 $obslnk obs.nc
# cdo -v -seldate,1970-01-01,2017-12-31 -remapbil,$simlnk obs.nc rstst.nc
cdo -v -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapbil,obs_fix2.nc $simlnk simh.nc
cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs_fix2.nc $simlnk simp.nc
mrz


#q = load(joinpath(wrkdir,"obs_mean.nc") ,myvar)

obs = load(joinpath(wrkdir,"obs_fix2.nc") ,myvar)
# Maps the time-mean average of ClimGrid C. If a filename is provided, the figure is saved in a png format.
contourf(obs;region="EU")

simh =load(joinpath(wrkdir,"simh.nc"),myvar)
contourf(simh; cm="jet")
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simh))

cm.plot(annualmax(obs))
cm.plot(annualmax(simh))

# q = temporalsubset(obs,(2001,6,1),(2001,6,1))
# q = q .* 24
# contourf(q)
# qs = temporalsubset(simh,(2001,6,1),(2001,6,2))
# contourf(qs)
contourf(obs)
contourf(simp)
out = qqmap(obs,simh,simp;)

contourf(simp-out,titlestr="bias correction effect", center_cs=true)

max_obs = annualmax(obs)
max_modelinterp = annualmax(simp)
max_modelqqmap = annualmax(out)
# Plots
plot(max_obs, label="OBS")
plot(max_modelinterp, label="REMO - interpolated")
plot(max_modelqqmap, label="REMO - bias corrected", titlestr = "Effect of bias correction on annual maximum values")


max_obs = annualmax(temporalsubset(obs,(1970,01,01),(2099,01,01)))
max_modelinterp = annualmax(temporalsubset(simp,(1970,01,01),(2099,01,01)))
max_modelqqmap = annualmax(temporalsubset(out,(1970,01,01),(2099,01,01)))

begin
    plot(max_obs, label="OBS")
    plot(max_modelinterp, label="REMO - interpolated")
    plot(max_modelqqmap, label="REMO - bias corrected", 
        titlestr = "Effect of bias correction on annual maximum values",
        filename = "tso.png")
end


##hm ok, naja, e-obs kann iwie nicht stimmen.
## ich versuch es mal mit v28

cm.plot(annualmean(out))
cm.plot(annualmax(out))
contourf(out)
contourf(simp)
cm.plot(annualsum(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

ClimateTools.write(out,"$myvar-cor.nc")

