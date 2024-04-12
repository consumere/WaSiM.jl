#QQ fix
myvar="rsds"
#now daily mean global radiation QQ.
wrkdir="/mnt/d/remo/qm/rsds"
cd(wrkdir)
ept = "/mnt/d/remo/cordex/eobs/"
using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)

@pyimport xarray as xr
ob = xr.open_dataset("rstst.nc")
@pyimport matplotlib.pyplot as plt
display(ob)
x=ob.where(ob[myvar] > 0,drop=true)
x.to_netcdf("obs_fix2.nc")
q = tt.load("obs_fix2.nc",myvar)
cm.plot(annualmax(q))
contourf(q)

#N = [8.14986035840001, 48.8498605556]
#poly_reg = [[NaN -65 -80 -80 -65 -65];[NaN 42 42 52 52 42]]
lg = q.longrid
lt = q.latgrid
N = [lg;lt]
(tt.spatialsubset)|>methods
su = tt.spatialsubset(q,N)
cm.mapclimgrid(q,"Europe")

#bash
"/mnt/d/remo/qm/rsds"
# obslnk="/mnt/d/remo/cordex/eobs/qq_ens_spread_crop.nc"
obslnk="/mnt/d/remo/cordex/eobs/qq_ens_spread_0.1deg_reg_v27.0e.nc"
#dims of remo data in t.gr
croptg="/mnt/d/remo/cordex/eobs/t.gr"
#cdo -v setgrid,$croptg $obslnk obs.nc
cdo -v chname,qq,rsds -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapcon,$croptg $obslnk obs.nc
cdo -v -mulc,24 -setmissval,-9999 obs.nc obs2.nc
cdo -v -mulc,12 -setmissval,-9999 obs.nc obs_12.nc
mobs="/mnt/d/remo/cordex/eobs/qq_ens_mean_v27_crop.nc"
#mobs="/mnt/d/remo/cordex/eobs/qq_ens_mean_0.1deg_reg_v27.0e"
cdo -v chname,qq,rsds -setmissval,-9999 -setcalendar,365_day $mobs obs_mean.nc

#obs = load(joinpath(wrkdir,"obs2.nc") ,myvar) #this errors
@pyimport xarray as xr
@pyimport matplotlib.pyplot as plt
ob = xr.open_dataset("obs_mean.nc")
display(ob)
x=ob.where(ob[myvar] > 0,drop=true)
#p = x[myvar].groupby("time.year").sum().plot()
x.to_netcdf("obs_fix2.nc")

x
#same as cm.plot annualmean
x["rsds"].mean("longitude").mean("latitude").groupby("time.year").mean().plot()
x["rsds"].isel(longitude=3,latitude=5).groupby("time.year").mean("time").plot()

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

cd("/mnt/d/remo/qm/corgrids")
cm.plot(annualsum(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

ClimateTools.write(out,"$myvar-cor.nc")


##crop it 
xmin,xmax,ymin,ymax= 8.14986035838497, 12.8498603395828 ,48.8498595763884, 50.7498594897639 ; 
xsize     = Int(48)
ysize     = Int(20)

function generate_matrix(xmin, xmax, ymin, ymax, xsize, ysize)
    xrange = range(xmin, stop=xmax, length=xsize)
    yrange = range(ymin, stop=ymax, length=ysize)

    matrix = Matrix{Tuple{Float64, Float64}}(undef, length(xrange), length(yrange))

    for (i, x) in enumerate(xrange)
        for (j, y) in enumerate(yrange)
            matrix[i, j] = x, y
        end
    end

    return matrix
end

# Given input
xmin, xmax, ymin, ymax = 8.14986035838497, 12.8498603395828, 48.8498595763884, 50.7498594897639
xsize = 48
ysize = 20

# Generate the matrix
result_matrix = generate_matrix(xmin, xmax, ymin, ymax, xsize, ysize)

findmax(result_matrix)
findmin(result_matrix)


pt="/mnt/d/remo/cordex/eobs/qq_ens_mean_0.1deg_reg_v28.0e.nc"
qq = load(pt,"qq")


py"""
import subprocess

def crop_netcdf(input_file, output_file, xmin, xmax, ymin, ymax):
    try:
        # Construct the CDO command for cropping
        cdo_command = f"cdo -v sellonlatbox,{xmin},{xmax},{ymin},{ymax} {input_file} {output_file}"

        # Run the CDO command using subprocess
        subprocess.run(cdo_command, shell=True, check=True)

        print(f"NetCDF file successfully cropped. Output saved to {output_file}")

    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        print("Failed to crop NetCDF file.")
"""

py"
input_file = $pt
output_file = 'qq_v28.nc'
xmin, xmax, ymin, ymax = 8.14986035838497, 12.8498603395828, 48.8498595763884, 50.7498594897639
"
py"print(input_file)"
#py"import os;os.chdir('/mnt/d/remo/cordex/eobs/');os.getcwd()"
py"
crop_netcdf(input_file, output_file, xmin, xmax, ymin, ymax)
"

pt="/mnt/d/remo/cordex/eobs/qq_v28.nc"
# qq = load(pt,"qq")#fails
@pyimport xarray as xr
#@pyimport matplotlib.pyplot as plt
ob = xr.open_dataset(pt)
display(ob)
myvar="qq"
x=ob.where(ob[myvar] > 0,drop=true)

x[myvar].mean("time").plot() #same as contourf
x.to_netcdf("obs_fix3.nc")
ds = NCDatasets.load!("obs_fix3.nc")
#NCDatasets.renameDim(ds,myvar,"rsds")
#NCDatasets.close(ds)
q = tt.load("obs_fix3.nc",myvar)
#tt.renameVar(q.data,myvar,"rsds")
cm.plot(annualmax(q))
contourf(q)

cd("/mnt/d/remo/cordex/eobs")
inf="qq_ens_spread_0.1deg_reg_v27.0e.nc"
of =  "qq_ens_spread_0.1deg_reg_v27.0e_crop.nc"  
py"
crop_netcdf(inf,of, xmin, xmax, ymin, ymax)
"

using Cmd

function crop_netcdf(input_file, output_file, xmin, xmax, ymin, ymax)
    try
        # Construct the CDO command for cropping
        cdo_command = `cdo -v sellonlatbox,$xmin,$xmax,$ymin,$ymax $input_file $output_file`

        # Run the CDO command using run
        run(cdo_command)

        println("NetCDF file successfully cropped. Output saved to $output_file")

    catch e
        println("Error: $e")
        println("Failed to crop NetCDF file.")
    end
end

crop_netcdf(inf,of, xmin, xmax, ymin, ymax)

pt="/mnt/d/remo/cordex/eobs/rsds-cor.nc"
ocor=load(pt,"rsds")

max_obs = annualmax(obs)
max_modelinterp = annualmax(simh)
max_modelqqmap = annualmax(ocor)
# Plots
plot(max_obs, label="OBS")
plot(max_modelinterp, label="REMO - interpolated")
plot(max_modelqqmap, label="REMO - bias corrected", titlestr = "Effect of bias correction on annual maximum values")


plot(annualmin(simp),label="REMO - interpolated")
plot(annualmin(ocor),label="REMO - bias corrected")
pt=/mnt/d/remo/qm/corgrids/rsds_rcm85.wa2


#println("taking PS Get-Clipboard and storing wslpath to clipboard ...")

function pew()
    try
        in = clipboard()
        wpath = replace(in, "\\" => "/")
        println("pt=$wpath")
        clipboard("pt=$wpath")
        return wpath
    catch e
        @error "smth errored $e"
    end
end

z = pe()
pt="D:/remo/qm/corgrids/tas_rcm85.wa2"

function pe()
    try
    inp = clipboard()
    wpath = replace(inp, "\\" => "/", "\"" => "")
    cmd = `wslpath -ua $wpath`
    ot = readchomp(pipeline(cmd))
    clipboard("$ot")
    println("$ot in clipboard!")
    return string(ot)
    catch e
        println("Error: $e")
        println("Failed to translate to wslpath.")
    end
end

k = pe()
k|>typeof
df = wread(k;skip=5)
@pj
pyjl.pyplot_df(df)
##all the same methods.
using BenchmarkTools
broadcast(findmax,eachcol(df))
map(findmax,eachcol(df))
findmax.(eachcol(df))
[findmax(df[!, col]) for col in names(df)]

bmarks = @benchmarkable begin
    broadcast(findmax, eachcol(df))
    map(findmax, eachcol(df))
    findmax.(eachcol(df))
    [findmax(df[!, col]) for col in names(df)]
  end
results = run(bmarks)
plot(results.times,label="t"    )
plot(results.gctimes,label="gc"    )

@btime findmax.(eachcol(df))

a = @btime broadcast(findmax, eachcol(df))
b = @btime map(findmax, eachcol(df))
c = @btime findmax.(eachcol(df))
d = @btime [findmax(df[!, col]) for col in names(df)]
#vef = [a,b,c,d]



###measured qq

glob("rsds")
inmet="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/"
v="radiation_ce24.txt"
joinpath(inmet,v)
using CSV
qmeas = waread2(inmet*v)
first(df)
qmeas|>dropmissing
qmeas|>propertynames
qmeas[!,Cols(r"wu"i,:date)]
qkg = qmeas[!,Cols(r"bad"i,:date)]|>dropmissing
cm.plot(annualmax(q))
pyjl.pyplot_df(qkg)
cd("/mnt/d/remo/qm/corgrids")
rdf = waread2("rsds_rcm85.wa")
qsim = rdf[!,Cols(r"bad"i,:date)]|>dropmissing
@wasim
wa.qplot(wa.mall(qkg,qsim))

df = wa.mall(qkg,qsim)
df = df[!,Not(:date)]
s = names(df)
#t2 = string.("q","\n",s[1],"|",s[2])
t2 = string.("obs|sim")
using StatsPlots
StatsPlots.plot( 
StatsPlots.qqplot(df[!,1],df[!,2], qqline = :fit), 
StatsPlots.qqplot(StatsPlots.Cauchy,df[!,1]), 
StatsPlots.qqnorm(df[!,2], qqline = :R),
title = t2)

StatsPlots.qqplot(df[!,1],df[!,2], 
qqline = :fit, markersize=.1, linewidth=2)
StatsPlots.histogram2d(df[!,1],df[!,2])
StatsPlots.histogram(df[!,1],df[!,2])


z = pe()
"/mnt/d/remo/qm/corgrids/rsds-cor.nc"
xo = load(z,"rsds")

annualmax(xo)|>plot
contourf(xo)|>plot
cdof(z)
lat()


myvar="qq"
pt="/mnt/d/remo/cordex/eobs/qq_ens_mean_v27_crop.nc"
q = tt.load(pt,myvar)
cm.plot(annualmax(q))
contourf(q)

#bash tests
ll
a=qq_v28.nc 
du $a
vardes $a
pydesc $a
#lat=latitude
lon=10.11
lat=50.44
#remapnn Nearest neighbor remapping
cdo -v -remapnn,lon=${lon}_lat=$lat $a qtout.nc

cdx qtout|hd
nctowasim
cdotowasim qtout.nc qq |cut -f1-5 > tou.txt
hd tou.txt
tsx tou.txt

dubplot("/mnt/d/remo/cordex/eobs/tou.txt")

#cdo -outputtab,date,lon,lat,value -remapnn,lon=${lon}_lat=$lat $a | xargs -I% nctowasim % qq > tou2.txt
#cdo -outputtab,date,value,lon,lat -remapnn,lon=${lon}_lat=$lat $a > tmp 
cdo -outputtab,date,value -remapnn,lon=${lon}_lat=$lat $a > tmp 
nctowasim tmp qq > tou2.txt
htlat

##julia appr.
sf = "/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
xd = lpro(sf)
wa.stplot(sf)

cd("/mnt/d/remo/cordex/eobs/v28")
xmin,xmax,ymin,ymax= 8.14986035838497, 12.8498603395828 ,48.8498595763884, 50.74985
ncs = nconly("")
for x in ncs
    println("cropping $x ...")
    inf = x
    of = replace(x,"_0.1deg_reg"=>"",".nc"=>"_crop.nc")
    crop_netcdf(inf, of, xmin, xmax, ymin, ymax)
end

sf = "/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/temp_1970.txt"
od = lpro(sf)
wa.stplot(sf)
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
rasterpath= glob(r"tg+.*cr")|>first
#st = Raster(rasterpath;lazy=true) #7GB 
st = Raster(rasterpath;)
dou = collect(extract(st, pnts;atol=0.1))
dx = DataFrame(dou)|>dropmissing
using GeoDataFrames
GeoDataFrames.getgeometrycolumns(dx)
#make station df loop
myvar = propertynames(dx)[2]
df = []
for i in 1:nrow(dx)
	k = tovec(dx,2)[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,myvar ]))
end
df
ds = innerjoin(df..., on= :date, makeunique=true)
nn = map(x->replace(x, 	"\xfc" => "ue",    "\xdf" => "ss",
r",.*" => "",r"-.*" => ""),string.(od.name))

nn = map(p->replace(p, 	"\xfc" => "ue",    "\xdf" => "ss",
r",.*" => "",r"-.*" => "_"),string.(od.name))
#nn = map(y->(y[1:5]*y[end]),string.(nn))
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
outname = string(myvar)*"_ext.wa"
wa.writewa(outname,ds)
op()

##now qq
rasterpath= glob(r"qq+.*cr")|>first
st = Raster(rasterpath;)
dou = collect(extract(st, pnts;atol=0.1))
dx = DataFrame(dou)|>dropmissing
GeoDataFrames.getgeometrycolumns(dx)
#make station df loop
myvar = propertynames(dx)[2]
df = []
for i in 1:nrow(dx)
	k = tovec(dx,2)[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,myvar ]))
end
df
#ds = innerjoin(df..., on= :date, makeunique=true)
ds = innerjoin(unique.(df, :date)..., on = :date, makeunique=true)
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
outname = string(myvar)*"_ext.wa"
wa.writewa(outname,ds)


##now all 
ras = glob(r"crop")
ras = ras[Not(1,3,5)]

for x in ras
    rasterpath = x 
    st = Raster(rasterpath;)
    dou = collect(extract(st, pnts;atol=0.1))
    dx = DataFrame(dou)|>dropmissing
    if nrow(dx) != nrow(od)
        @error "something wrong here..."
        return
    end
    myvar = propertynames(dx)[2]
    df = []
    for i in 1:nrow(dx)
        k = tovec(dx,2)[i]
        push!(df, DataFrame([Rasters.lookup(k,1),
            (k.data)],[:date,myvar ]))
            #Float64.(k.data)],[:date,myvar ]))
    end
    
    ds = innerjoin(df..., on= :date, makeunique=true)
    rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
    outname = string(myvar)*"_ext.wa"
    wa.writewa(outname,ds)
    println("done $outname")
end

cmd = `wsl wslpath -a .`
wsl_path = readchomp(pipeline(wsl_cmd))

run(`head -n 5 $sf > station.wa`)
#run(pipeline(`dothings`, stdout="out.txt", stderr="errs.txt"))
run(pipeline(`head -n 5 $sf`, 
    stdout="station_crds.txt", 
    stderr="errs.txt"))


crd = DataFrame(dx[!,1])
crd.name .= nn
rename!(crd,1=>:lon,2=>:lat)
crd
crdf = permutedims(crd[!,[:name,:lon,:lat]])
writedf("cords.wa",crdf)

dy = dfr(outname)
baryrsum(dy)
#hydromon(dy)
"station_crds.txt"
"cords.wa"
# in julia, i got 2 files: "station_crds.txt" and
# "cords.wa". i need to replace line 3 and 4 from column 5:end from one file to anthoer

# Read the contents of both files
station_crds = readdlm("station_crds.txt")
cords_wa = readdlm("cords.wa",'\t')


# Read the contents of both files
station_crds = readlines("station_crds.txt")
cords_wa = readlines("cords.wa")

# Split each line into fields
station_crds = [split(line, '\t') for line in station_crds]
cords_wa = [split(line, '\t') for line in cords_wa]

# Replace fields in lines 3 and 4 in station_crds with fields from cords_wa
station_crds[3][5:end] = cords_wa[3]
station_crds[4][5:end] = cords_wa[4]
# Join the fields back together
station_crds = [join(fields, '\t') for fields in station_crds]
# Write the modified contents back to "station_crds.txt"
open("station_crds2.txt", "w") do f
    for line in station_crds
        println(f, line)
    end
end

ras = glob(r"crop")
pg = tt.load(ras[4],"rr")
pg = tt.load(ras[5],"tg")
#cm.plot(annualsum(pg))
cm.plot(annualmax(pg),titlestr="temperature")
cm.contourf(pg,titlestr="temperature")
cm.contourf(pg,
region="auto",caxis=(4, 12),
cm="jet",titlestr="mean obs temperature")

#rad = tt.load(ras[3],"qq") #fails, so fix needed
using PyCall
pygui(true)
@pyimport xarray as xr
ob = xr.open_dataset(ras[3])
@pyimport matplotlib.pyplot as plt
display(ob)
ob = ob.rename(Dict("qq"=>"rsds"))
#PyDict(d::Dict{K,V})
x = ob.where(ob["rsds"] > 0,drop=true)
x.to_netcdf("rsds_eobs.nc")
q = tt.load("rsds_eobs.nc","rsds")
cm.plot(annualmax(q))
cm.contourf(q)

ob["rsds"].mean("longitude").mean("latitude").groupby("time.year").mean().plot()
ob["rsds"].mean()
ob["rsds"].mean("longitude").mean("latitude").mean()
ob["rsds"].isel(time=-1).plot()
ob["rsds"].mean("time").plot() #this is xr contourf
#covariance + correlation
dac=xr.cov(ob["rsds"].mean("time"),x["rsds"].mean("time"))
dac=xr.corr(ob["rsds"].mean("time"),x["rsds"].mean("time"))

tt.regrid

#cdo shifttime,-15days infile outfile
# cdo -v -setcalendar,365_day 
# -seldate,1970-01-01,2017-12-31 -remapbil,obs_fix2.nc $simlnk simh.nc

function ncremap_o(input_file, remapfile, output_file, 
    startdate, enddate;oldname::String, newname::String)
    try
        # Construct the CDO command for cropping
        # -setcalendar,365_day 
        #-P <nthreads>  Set number of OpenMP threads
        if oldname
            cdo_command = `cdo -P 4 chname,$oldname,$newname -seldate,$startdate,$enddate -remapbil,$remapfile $input_file $output_file`
        else
            cdo_command = `cdo -P 4 -seldate,$startdate,$enddate -remapbil,$remapfile $input_file $output_file`
        end
        # Run the CDO command using run
        run(cdo_command)
        println("NetCDF file successfully processed. 
            Output saved to $output_file")

    catch e
        println("Error: $e")
        println("Failed to remap $output_file.")
    end
end

simlnk="/mnt/d/remo/cordex/wgs/proj_rsds_hist+rcp85_wgs84.nc"
output_file = "simh.nc"
ncremap(simlnk,output_file,"1970-01-01","2020-12-31")
#cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs_fix2.nc $simlnk simp.nc
output_file = "simp.nc"
ncremap(simlnk,output_file,"2021-01-01","2099-12-31")

simh = tt.load("simh.nc","rsds")
simp = tt.load("simp.nc","rsds")
obs  = tt.load("rsds_eobs.nc","rsds")


out = qqmap(obs,simh,simp;)
dk = ncdf("rsds_eobs.nc")
sim = ncdf("simh.nc")
dubplot(sim)
dubplot(dk)
dk.date .= Dates.Date.(dk.date)
sim.date .= Dates.Date.(sim.date)
sim.rsds .= Float64.(sim.rsds)
dk.rsds .= Float64.(dk.rsds)
#t = mall(sim,dk)
baryrsum(mall(sim,dk))
qplot(mall(sim,dk))
dpr(mall(sim,dk))


max_obs = annualmax(obs)
max_modelinterp = annualmax(simh)
max_modelqqmap = annualmax(simp)

begin
    cm.plot(max_obs, label="OBS")
    cm.plot(max_modelinterp, label="REMO - historical")
    cm.plot(max_modelqqmap, label="REMO - future", 
        titlestr = "a")
end



dfs = broadcast(ncdf,["simh.nc","simp.nc","rsds_eobs.nc"])
df2 = map(x->Dates.Date.(x.date),dfs)
map(x->x.date .= Dates.Date.(x.date),dfs)
map(x->x.rsds .= Float64.(x.rsds),dfs)
dout = innerjoin(unique.(dfs[Not(2)], :date)..., on = :date, makeunique=true)
dpr(dout)
bargroup(dout)



#ext(9.27616867640704, 11.8984838654896, 48.9648300948859, 50.5785625189367)
xmin,xmax,ymin,ymax = 9.27616867640704, 11.8984838654896, 48.9648300948859, 50.5785625189367
xmin,xmax,ymin,ymax = 9, 12, 49, 51
#ceil(xmin)
ncs = nconly("crop")
for x in ncs
    println("cropping $x ...")
    inf = x
    of = replace(x,"_crop.nc"=>"_rcm.nc")
    crop_netcdf(inf, of, xmin, xmax, ymin, ymax)
end



# "/mnt/d/Wasim/regio/out/rc200/x31/loc2/"|>cd
# heat(r"Wol")

"/mnt/c/users/chs72fw/documents/EFRE_GIS/Hydrologie/mhm/mhm-wsl/mhm-v5.12.0/test_domain_2/output"|>cd
ls()
#r = Raster("mRM_Fluxes_States.nc")
r = xr.open_dataset("mRM_Fluxes_States.nc")
r["Qrouted"].mean("time").plot()
r.close()

r = xr.open_dataset("mHM_Fluxes_States.nc")
r["QB"].mean("time").plot()
r["Qsm"].mean("time").plot()
r["preEffect"].mean("time").plot()
r.close()
cdu()
ls()
readdir("input/morph")
r = xr.open_dataset("input/morph/dem.asc")
r["band_data"].plot()

pt="/mnt/c/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/mhm/mhm-wsl/mhm-v5.12.0/test_domain_2/input/gauge/45.txt"

using PyCall
@pyimport xarray as xr
#df = CSV.read(pt,DataFrame;skipto=7,header=1,delim=" ")
@vv "mhm"
cd("/mnt/c/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/mhm/mhm-wsl/mhm-v5.12.0/test_domain_2")
r = xr.open_dataset("input/meteo/ssrd.nc")
r["ssrd"].mean("time").plot()
r["ssrd"].mean("xc").mean("yc").plot()
r.close()

p=r["ssrd"].mean("xc").mean("yc").plot()
@pyimport matplotlib.pyplot as plt
display(only(p))


#fix show backend
@pyimport matplotlib
matplotlib.use("TkAgg")

p=r["ssrd"].mean("xc").mean("yc").plot()
plt.show()

pt=src_path*"/pyjl.jl";
include(pt)
cd("/mnt/d/Wasim/regio/out/rc200/x12/loc6")
pyjl.xrp(r"sb1")
pyjl.xrp(r"sb05")


cd("/mnt/d/Wasim/regio/out/rc200/x12/loc7")
s = nconly("sb0")
r = xr.open_dataset(s[3])
k = r.keys()|>collect |>first
#r[k].mean("x").mean("y")
r[k].transpose.plot()
lyr=0
maskval=0.1
dx = r.isel(t=lyr).to_array()
y=dx.where(dx.values .> maskval)
px = y.transpose().plot(cmap="cividis");
px.colorbar.set_label("mm")        #cbar
plt.grid(color="white", linestyle="-", linewidth=0.5)
plt.title("soil moisture")
plt.xlabel("x")
plt.ylabel("y")


#plt = pyimport("matplotlib.pyplot")


dx.where(dx.values .> maskval).transpose().
    plot(cmap="cividis").colorbar.set_label("mm").figure.title("soil moisture")
    grid(color="white", linestyle="-", linewidth=0.5).xlabel("x").ylabel("y")   