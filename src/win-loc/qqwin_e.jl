raw"E:\qq"|>cd
using Pkg
Pkg.activate(".")
Pkg.gc(; collect_delay=Dates.Day(0))
#Pkg.update()
#Pkg.instantiate() #to install all recorded dependencies.
Pkg.status()
pt="C:/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
include(pt)
pwc()
# function loadcdo()
#     pt="C:/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
#     include(pt)
# end
#loadcdo()

printstyled("Bias Correction of REMO Data by using EOBS grids.\nReprojected to EPSG 25832\n",color=:red, bold=true)

using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall; pygui(true) 
v = rglob(r".nc")
fdi()
#mkdir("jlcor")
cd("jlcor")
grep("obs",v)

#m=tt.load(first(grep("obs",v)),"pre")
#cm.plot(annualsum(m))
##tt.griddata
    cvar = Symbol("pre")
    obspt = "E:/qq/eobs/proj/pre-obs.nc"
    obsras = nothing
    sim = nothing
    
    m = tt.load(obspt,string(cvar))
    obsras = temporalsubset(m,(1960,01,01),(2020,12,31))
    #cm.plot(annualsum(obsras))
    #cm.contourf(obsras)
    #simlnk="D:/remo/cordex/wgs/proj_pre_hist+rcp85_wgs84.nc"
    simlnk="D:/remo/cordex/wgs/utm/proj_pre_hist+rcp85_utm.nc"
    sim = tt.load(simlnk,string(cvar))
    @doc griddata
    #Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
    #where A and B are ClimGrid.
    mytime = [(1960,01,01),(2020,12,31)]
    futime = [(1980,01,01),(2099,12,31)]
    simh = temporalsubset(sim,mytime...)
    simp = temporalsubset(sim,futime...)
    
    simh = tt.griddata(simh,obsras) #
    tt.write(simh,"pre-simh-remo-jl.nc") #
    simp = tt.griddata(simp,obsras) #ETA: 0:10
    tt.write(simp,"pre-simp-remo-jl.nc") #NO mem-error.

    simp.timeattrib

    pwd()
    @doc tt.biascorrect_extremes
    #biascorrect_extremes(obs::ClimGrid, 
    #ref::ClimGrid, fut::ClimGrid; 
    #detrend=false, window::Int=15, 
    #rankn::Int=50, 
    #thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat(), gev_params::DataFrame, frac=0.25, power=1.0)
    @time outextr = tt.biascorrect_extremes(
        obsras,simh,simp;thresnan=0.01) #139.965272 seconds
     #cm.contourf(outextr)
    cm.plot(annualsum(outextr))
    
    
    tt.write(outextr,"pre-jl.nc") #qmextr result
    
    #out = qqmap(obsras,simh,simp) #takes 1h
    out = qqmap(obsras,simh,simp)
    tt.write(out,"pre-qq.nc") #
    lat()
    #pwd()
    cm.contourf(out)
    x0 = out.data[:,:,1]
    #heatmap(x0.data)
    #cm.plot(x0.data)
    
    max_obs = annualsum(obsras)
    max_modelinterp = annualsum(sim)
    max_modelqqmap = annualsum(out)
    max_modelexqqmap = annualsum(outextr)

    begin
        cm.plot(max_modelinterp, label="REMO - interpolated")
        cm.plot(max_modelqqmap, label="REMO - bias corrected")
        cm.plot(max_modelqqmap, label="REMO - exbias corrected", 
        titlestr = "Auswirkung der Bias-Korrektur auf die Jahresniederschlagssummen",
        filename = "$cvar-QQ.png")
    end
    
    begin
        cm.plot(max_obs, label="EOBS")
        cm.plot(max_modelinterp, label="REMO - interpolated")
        #cm.plot(max_modelqqmap, label="REMO - bias corrected")
        cm.plot(max_modelqqmap, label="REMO - exbias corrected", 
        titlestr = "Auswirkung der Bias-Korrektur auf die Jahresniederschlagssummen",
        filename = "$cvar-sum.png")
    end

du()

grep("tas",v)
cvar = Symbol("tas")

s = joinpath("D:/remo/cordex/eobs/v28/proj/","tas-obs.winlist")
readlines(s)[20]
#m=tt.load(first(grep("tas",v)),"tas")

tp="D:/remo/cordex/eobs/v28/proj/tas-obs.nc" #das ist 25832
# tp="D:/remo/cordex/eobs/v28/tas/tas_obs.nc" #das ist 4326
m = tt.load(tp,"tas")
obsras = temporalsubset(m,mytime...)
cm.plot(annualmean(m))
# cm.contourf(m)
#simlnk="D:/remo/cordex/wgs/proj_tas_hist+rcp85_wgs84.nc"
simlnk="D:/remo/cordex/wgs/utm/proj_tas_hist+rcp85_utm.nc" #das ist 25832
sim = tt.load(simlnk,string(cvar))
#cm.contourf(sim)
cm.plot(annualmin(sim))
@doc griddata
#Interpolate ClimGrid A onto the lon-lat grid of ClimGrid B, 
#where A and B are ClimGrid.
simh = temporalsubset(sim,mytime...)
simp = temporalsubset(sim,futime...)
simh = tt.griddata(simh,obsras) #7min
write(simh,"tas-simh-remo-jl.nc") #
simp = tt.griddata(simp,obsras) #7min
write(simp,"tas-simp-remo-jl.nc") #7min
out = qqmap(obsras,simh,simp;method="Additive")
write(out,"tas-qq.nc")          #


max_obs = annualmean(obsras)
max_modelinterp = annualmean(sim)
max_modelqqmap = annualmean(out)

begin
    cm.plot(max_obs, label="EOBS")
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Auswirkung der Bias-Korrektur auf die Jahresmitteltemperatur",
    filename = "$cvar-QQ.png")
end

#now wind #############
#simlnk="D:/remo/cordex/wgs/proj_sfcWind_hist+rcp85_wgs84.nc"
simlnk="D:/remo/cordex/wgs/utm/proj_sfcWind_hist+rcp85_utm.nc"
cvar = Symbol("sfcWind")
sim = tt.load(simlnk,string(cvar))
simh = temporalsubset(sim,mytime...)
simp = temporalsubset(sim,futime...)

m = tt.load("D:/remo/cordex/eobs/v28/proj/sfcWind-obs.nc","sfcWind")
#cm.contourf(m)
obsras = temporalsubset(m,mytime...)

simh = tt.griddata(simh,obsras) 
write(simh,"$cvar-simh-remo-jl.nc") #
simp = tt.griddata(simp,obsras) 
write(simp,"$cvar-simp-remo-jl.nc")
out = qqmap(obsras,simh,simp;method="Additive")
write(out,"$cvar-qq.nc")          #

max_obs = annualmean(obsras)
max_modelinterp = annualmean(sim)
max_modelqqmap = annualmean(out)

begin
    cm.plot(max_obs, label="EOBS")
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Auswirkung der Bias-Korrektur auf die Mittlere Windgeschwindigkeit",
    filename = "$cvar-QQ.png")
end


lat()
#now rh
#simlnk="D:/remo/cordex/wgs/proj_rh_hist+rcp85_wgs84.nc"
sim=nothing
simlnk="D:/remo/cordex/wgs/utm/proj_rh_hist+rcp85_utm.nc"
cvar = Symbol("rh")
sim = tt.load(simlnk,string(cvar))
cm.plot(annualmean(sim))

simh = temporalsubset(sim,mytime...)
simp = temporalsubset(sim,futime...)
#m = tt.load("D:/remo/cordex/eobs/v28/proj/rh-obs.nc","rh")
m = tt.load("E:/qq/eobs/proj/rh-obs.nc","rh")
obsras = temporalsubset(m,mytime...)
cm.plot(annualmean(obsras))
simh = tt.griddata(simh,obsras)
#cm.plot(annualmean(simh))
#pwd()
write(simh,"$cvar-simh-remo-jl.nc") #
simp = tt.griddata(simp,obsras)
write(simp,"$cvar-simp-remo-jl.nc")
#out = qqmap(obsras,simh,simp;method="Multiplicative")
cm.plot(annualmean(out)) #empty
out = nothing
#out = qqmap(obsras,simh,simp;method="Additive")
#out = tt.biascorrect_extremes(obsras,simh,simp;thresnan=0.001) #empty
out = nothing
out = tt.biascorrect_extremes(obsras,simh,simp)
out = qqmap(obsras,simh,simp;method="Additive",detrend=false,thresnan=0.01)
cm.plot(out) #klappt
cm.plot(annualmean(out))
write(out,"$cvar-qq.nc")          #

cm.plot(annualmean(obsras))
cm.plot(annualmean(simh))
cm.plot(annualmean(simp))


max_obs = annualmean(obsras)
max_modelinterp = annualmean(sim)
max_modelqqmap = annualmean(out)
begin
    cm.plot(max_obs, label="EOBS")
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Auswirkung der Bias-Korrektur der relativen Luftfeuchtigkeit",
    filename = "$cvar-QQ.png")
end

pfix()

#radiation 
#simlnk="D:/remo/cordex/wgs/proj_rsds_hist+rcp85_wgs84.nc"

sim = nothing
simlnk="D:/remo/cordex/wgs/utm/proj_rsds_hist+rcp85_utm.nc"
cvar = Symbol("rsds")
sim = tt.load(simlnk,string(cvar))
simh = temporalsubset(sim,mytime...)
simp = temporalsubset(sim,futime...)

obsras = nothing
m = tt.load("D:/remo/cordex/eobs/v28/proj/rsds-obs.nc","rsds")
obsras = temporalsubset(m,mytime...)

simh = tt.griddata(simh,obsras)
write(simh,"$cvar-simh-remo-jl.nc") #
simp = tt.griddata(simp,obsras)
write(simp,"$cvar-simp-remo-jl.nc")
out = qqmap(obsras,simh,simp;method="Multiplicative")
write(out,"$cvar-qq.nc")          #

max_obs = annualmean(obsras)
max_modelinterp = annualmean(sim)
max_modelqqmap = annualmean(out)
begin
    cm.plot(max_obs, label="EOBS")
    cm.plot(max_modelinterp, label="REMO - interpolated")
    cm.plot(max_modelqqmap, label="REMO - bias corrected", 
    titlestr = "Auswirkung der Bias-Korrektur der Globalstrahlung",
    filename = "$cvar-QQ.png")
end


##########check####################
# using PyCall
# @pyimport rioxarray as rx
# @pyimport geopandas as gpd
# @pyimport xarray as xr

@pyjl

"E:/qq/jlcor/"|>cd
pwd()
pwc()
op()
@pyimport xarray as xr
@doc pyjl.doy
#v = glob("qq")
@edit pyjl.fdoy(v[3])
v = nconly("qq")
pyjl.fdoy(v[1])
pyjl.fdoy(v[2])
pyjl.fdoy(v[3])
pyjl.fdoy(v[4])
pyjl.fdoy(v[5])

#rm(v[2],force=true)

ds = xr.open_dataset("tas-qq.nc")
ds = ds.drop(["lat","lon"])
#ds.Regular_longitude_latitude
#ds = ds.drop("Regular_longitude_latitude")
#ds.variables["grid_mapping"]
@pyimport rioxarray as rx
# if haskey(ds.attrs, "grid_mapping")
#     delete!(ds.attrs, "grid_mapping")
# end
#ds = ds.rio.write_crs("EPSG:25832")
# ds.Regular_longitude_latitude
# ds.Regular_longitude_latitude.attrs
# ds = ds.rename(Dict("Regular_longitude_latitude"=>"spatial_ref"))
# ds.spatial_ref.attrs

# if haskey(ds.attrs, "grid_mapping")
#     ds = ds.drop_vars("grid_mapping", errors="ignore")
# end

ds = xr.open_dataset("tas-qq.nc")
ds = ds.drop(["lat","lon"])
ds = ds.rename(Dict("time" => "t"))
ds = ds.transpose("x","y","t")
ds = ds.rio.write_crs("EPSG:25832")
if haskey(ds["tas"].attrs, "grid_mapping")
    ds["tas"].attrs = delete!(ds["tas"].attrs, "grid_mapping")
end
ds["tas"].attrs
ds.attrs
ds = ds.rename(Dict("Regular_longitude_latitude"=>"spatial_ref"))
ds.spatial_ref.attrs
ds.to_netcdf("tas-qq2.nc")

"E:/qq/jlcor/"|>cd
#R"ls()"
v=glob(r"-qq2[.]")
broadcast(rm,v)

using PyCall
xr = pyimport("xarray")
rio = pyimport("rioxarray")

#ds = xr.open_dataset("tas-qq2.nc")
# for var in ds.variables
#     ds[var].encoding["_FillValue"] = -9999
# end


#ds["tas"].encoding = setindex!(ds["tas"].encoding, -9999, "_FillValue")

#double pr(time, x, y) ; !!!!

function nctowasim(filename::String;)
    ds = xr.open_dataset(filename)
    variable = last(collect(ds.keys()))
    if variable=="Regular_longitude_latitude"
        @error "var is Regular_longitude_latitude"
        return
    end
    ds = ds.drop(["lat","lon"])
    ds = ds.rio.write_crs("EPSG:25832")

    # fill_value = -9999
    # ds[variable].encoding =  setindex!(ds[variable].encoding, fill_value, "_FillValue")
    # ds[variable].attrs = setindex!(ds[variable].attrs, fill_value, "missing_value")
    
    if "Regular_longitude_latitude" in ds.variables
        ds = ds.drop_vars("Regular_longitude_latitude")
    end
    if haskey(ds[variable].attrs, "grid_mapping")
        ds[variable].attrs = delete!(ds[variable].attrs, "grid_mapping")
    end
    #ds = ds.rename(Dict("Regular_longitude_latitude"=>"spatial_ref"))
    #ds = ds.rename(Dict("time" => "t"))
    #ds = ds.transpose("x","y","t")
    ds = ds.transpose("time","x","y")
    #ds = ds.transpose("x","y","time")
    outfn = replace(filename, ".nc" => "2.nc")
    ds.to_netcdf(outfn)
    ds.close()
    @info "stored to $outfn !"
end

v=nconly(r"-qq[.]")
#v = v[Not(end)] #remove last element(tas) cause i already did this
for x in v
    nctowasim(x)
end
#now C:\Users\Public\Documents\Python_Scripts\julia\make_lists-jlcor.jl

function xrlist2(;cwd::AbstractString=pwd(),xcrd="x",ycrd="y",tvec="t",gridres::Float64=8317.01, xmatch="qq2",suffix=".winlist")
    files = filter(x -> contains(x, xmatch) && endswith(x, ".nc"), readdir())
    if any(x -> occursin(xmatch, x) && endswith(x, suffix), readdir())
        println("Corresponding files already exist. Exiting now!")
        return
    end

    for var in files
        global onam = string(replace(var,".nc" => suffix))
        ds = xr.open_dataset(var)
        println("Processing $var ...")
        # # Access the coordinate information
        #xf = ds["lon"].values|>first
        xf = ds[xcrd].values|>first
        #yf = ds["lat"].values|>first
        #xf = round(xf, digits=3)
        yf = ds[ycrd].values|>first
        #yf = round(yf, digits=3)
        gridfile_content = ["GRID\tLIST\tof\t$var", #var is var in files
        "YY\tMM\tDD\tHH\t$xf",
        "YY\tMM\tDD\tHH\t$yf",
        "YY\tMM\tDD\tHH\t$gridres",
                            "YY\tMM\tDD\tHH\tListe\n"]
        open(onam, "w") do out
            write(out, join(gridfile_content, "\n"))
        end
        #timevec = ds["time"].values
        timevec = ds[tvec].values
        dts = nothing
        try
            dts = [Dates.DateTime.(t.year, t.month, t.day, t.hour, t.minute) for t in timevec]
        catch
            @warn "timevec is not ctftime"
            dts = timevec.astype("M8[m]").tolist()
        end      
        param = ds.keys()|>collect|>last
        dk = DataFrame(YY=year.(dts),MM=month.(dts),DD=day.(dts),HH=hour.(dts),
            #link="stateini/jlcor/$var",
            #link="D:/remo/qm/corgrids/jlcor/$var",
            link = joinpath(cwd,var),
            param = "<$param>",
            cnt = (1:length(dts)).-1)
        dout = hcat(dk[:, 1:4], [join(row, "") for row in eachrow(dk[:, 5:end])])
        #append to onam file
        CSV.write(onam,dout, 
            transform = (col, val) -> something(val, missing),
            delim="\t",append=true)  
            println("$param written to $onam !")
            println("$onam done!")
    end
    println("Lists for WaSiM Input are written to the current directory. All done!")
end

xrlist2(;tvec="time",gridres=8317.01012437603)

#68731/8317.01012437603

#rmlat()
npplat()
v=glob(r"-qq2.winlist")
broadcast(rm,v)


#rcm_y5.yml läuft.