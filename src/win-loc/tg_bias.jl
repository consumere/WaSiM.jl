#
#from /mnt/c/Users/Public/Documents/Python_Scripts/julia/dscale-ubu.jl   
myvar="tg"
using ClimatePlots
using ClimateTools
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)
using Dates

wrkdir="/mnt/d/remo/qm/tas"
cd(wrkdir)
ept = "/mnt/d/remo/cordex/eobs/"
ssup()
#@gl string.(myvar)
ls()
tg="/mnt/d/remo/cordex/eobs/tg_ens_spread_crop.nc"
tg = load(tg,"tg")
cm.plot(annualmean(tg))


#bash
"/mnt/d/remo/qm/tas"
# obslnk="/mnt/d/remo/cordex/eobs/qq_ens_spread_crop.nc"
obslnk="/mnt/d/remo/cordex/eobs/tg_ens_spread_crop.nc"
#dims of remo data in t.gr
croptg="/mnt/d/remo/cordex/eobs/t.gr"
#cdo -v setgrid,$croptg $obslnk obs.nc
cdo -v chname,tg,tas -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapcon,$croptg $obslnk obs-tg.nc
sf tas
#simlnk="/mnt/d/remo/cordex/wgs/proj_tas_hist+rcp85_wgs84.nc"
#cdo -v -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapbil,obs_fix2.nc $simlnk simh.nc
#cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs_fix2.nc $simlnk simp.nc

simh=load("/mnt/d/remo/qm/tas/simh.nc","tas")
simp=load("/mnt/d/remo/qm/tas/simp.nc","tas")
cm.plot(annualmean(tg),label = "tg")
cm.plot(annualmean(simp),label = "simp")

out = qqmap(tg,simh,simp;method="Additive")
cm.plot(annualmean(tg),label = "tg")
cm.plot(annualmean(out) ,label = "bias")


ClimateTools.write(out,"$myvar-cor.nc")


z = @gl "Bruecken"

read_dlm(z,)
df = DelimitedFiles.readdlm(z)
# Remove first row
#df[1:end.≠1,1:end]
#df[2:end,:]


p = plot(od.geometry);
using CSV,ArchGDAL,GeoFormatTypes

function lpro(x::Union{String,DataFrame})
    if x isa String
        fl = CSV.read(x,DataFrame;limit=4)
    else
        fl = x
    end
    xc = fl[2,5:end]|>collect
    yc = fl[3,5:end]|>collect
    pts = ArchGDAL.IGeometry[]
    for i in 1:length(xc)
        pt = ArchGDAL.createpoint([xc[i],yc[i]])
        pt = ArchGDAL.reproject(pt,
        GeoFormatTypes.EPSG(25832),
        GeoFormatTypes.EPSG(4326))

            # Proj.CRS("EPSG:25832"),
            # Proj.CRS("EPSG:4326"))
        #ProjString("+proj=longlat +datum=WGS84 +no_defs"))
        push!(pts,pt)
    end
    od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)

    return od    
end

od = lpro(z)
#plot(od.geometry)
od.geometry

#bash
xrcrds tg-cor.nc 9.7835 50.308


t = "/mnt/d/remo/cordex/wgs/utm/cropped/proj_tas_hist+rcp85_utm.crop.nc"

tas = load(t,"tas")
cm.contourf(tas)
cm.plot(tas)
#make obs from proj cropped.
lk="/mnt/d/remo/cordex/eobs/proj/tx_ens_mean_utm.nc"
cd "/mnt/d/remo/cordex/wgs/utm/cropped/"
cat gbgrid 
cdo vardes $lk
cdo -v chname,tx,tas -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapcon,gbgrid $lk obs-tas.nc
#Unsupported projection coordinates (Variable: tx)
xm    = 512975.003
ym    = 5545379.031
crop_netcdf
cdo -v chname,tx,tas -setcalendar,365_day -seldate,1970-01-01,2017-12-31 tx-proj.nc obs-tas.nc
cd("/mnt/d/remo/cordex/wgs/utm/cropped")
obs = load("obs-tas.nc","tas")
#su = tt.spatialsubset(obs,tas)
sim = load("proj_tas_histrcp85_utm.crop.nc","tas")
sy=(1970,01,01)
ey=(2016,12,31)
fy=(2017,01,01)
ne=(2099,12,31)
simh = tt.temporalsubset(sim,sy,ey)
simp = tt.temporalsubset(sim,fy,ne)

out = qqmap(obs,simh,simp;method="Additive")
cm.plot(annualmax(simh) ,label = "past")
cm.plot(annualmax(obs),label = "e-obs")
cm.plot(annualmax(simp) ,label = "future")
cm.plot(annualmax(out) ,label = "bias-corrected")

tt.write(out,"tas-cor.nc")

#cm.contourf(out)

cd(wrkdir)
import CSV: read as rd
rd(z,DataFrame)
m = map(x->(x=>Int64),[:YY,:MM,:HH,:DD])
#vcat(m, :5 => Float64)
df = rd(z,DataFrame;types=Dict(m))|>f->dropmissing(f,:YY)
df.date = df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
select!(df,Not(1:4))
DataFrames.metadata!(df, "filename", z, style=:note)

#dateformat=Dates.DateFormat("yy mm dd hh"))
#df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");

function fread(z::Union{Regex,AbstractString})
    # if isa(z,AbstractString)
    #     v = filter(file -> occursin(Regex(z,"i"),file), readdir());
    #     z = v[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt|log"i),v)]|>first        
    if isa(z,Regex)
        v = filter(file -> occursin(z,file), readdir());
        z = v[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt|log"i),v)]|>first
    end 
    println("loading $z ...")   
    m = map(x->(x=>Int64),[:YY,:MM,:HH,:DD])
    #df = CSV.read(z,DataFrame;types=Dict(m))|>f->dropmissing(f,:YY)
    ms = ["-9999","-9999.0","lin", "log", "--"]
    df = CSV.read(z,DataFrame;
        delim="\t", header=1, 
        normalizenames=true, 
        missingstring=ms,
        types=Dict(m),
        stripwhitespace=true)|>f->dropmissing(f,:YY)
    df.date = df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
    select!(df,Not(1:4))
    DataFrames.metadata!(df, "filename", z, style=:note)
end

m="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/schweinhof_1970_0"
df=fread(m)
@pj
pyjl.pyplot_df(df)

cdof(df)
dx = fread(r"Schw")


Dates.year(2010)|>typeof



df = wread(z;skip=5)
first(df)
plot(df.date,df.x5) #thats a pyplot
dy=yrsum(df)
op()
obs=load("../prec/obs.nc","pre")
contourf(obs)
pr=load("../prec/pre_delta_up2.nc","pre")
contourf(pr;cm="turbo")

cm.plot(annualsum(pr))
nconly("pre")

pcm=load("../prec/pre_cm.nc","pre")
cm.plot(annualsum(pcm))

# @pyimport xarray as xr
# @pyimport matplotlib.pyplot as plt
# ob = xr.open_dataset("obs_mean.nc")
# display(ob)
# x=ob.where(ob[myvar] > 0,drop=true)
# x.to_netcdf("obs_fix2.nc")

#same as cm.plot annualmean
@pyimport xarray as xr
nconly("pre")
x = xr.open_dataset("pre_delta_up.nc")
x[myvar].mean("longitude").mean("latitude").groupby("time.year").mean().plot()
x[myvar].isel(longitude=3,latitude=5).groupby("time.year").mean("time").plot()
x[myvar].isel(time=5).plot()


simh =load(joinpath(wrkdir,"simh.nc"),myvar)
contourf(simh; cm="jet")
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simh))

#pr=load("../prec/pre_delta_up2.nc","pre")
out = load("pre_cm.nc","pre")
pt="/mnt/d/remo/qm/corgrids/pre-cor.nc"
out = load(pt,"pre")

contourf(simp-out,titlestr="bias correction effect", center_cs=true)

max_obs = annualmax(obs)
max_modelinterp = annualmax(simp)
max_modelqqmap = annualmax(out)
# Plots
plot(max_obs, label="E-OBS")
plot(max_modelinterp, label="REMO - interpolated")
plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Effect of bias correction on annual maximum values")


max_obs = annualmax(temporalsubset(obs,(1970,01,01),(2015,01,01)))
max_modelinterp = annualmax(temporalsubset(simp,(1970,01,01),(2015,01,01)))
max_modelqqmap = annualmax(temporalsubset(out,(1970,01,01),(2015,01,01)))

begin
    plot(max_obs, label="E-OBS")
    plot(max_modelinterp, label="REMO - interpoliert")
    plot(max_modelqqmap, label="REMO - biaskorrigiert", 
        titlestr = "Auswirkung der Biaskorrektur auf die jährlichen Maximalwerte des Niederschlags",
        filename = "pre-check.png")
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

