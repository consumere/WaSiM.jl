#pycall
using PyCall
hombr()
cd("m6")
using Conda
Conda.add("statsmodels")
sm = pyimport("statsmodels.api")

spector_data = sm.datasets.spector.load()
spector_data["exog"] = sm.add_constant(spector_data["exog"], prepend=false)
model = sm.OLS(spector_data["endog"], spector_data["exog"])
res = model.fit()
res.summary()


lk="/mnt/d/Wasim/Goldbach/stateini/evaporation_modis_500/2022/ctlcum"
lk="d:/Wasim/Goldbach/stateini/evaporation_modis_500/2022/ctlcum"
using DataFrames, StatsPlots, GLM, CSV, Dates
nd = readdf(lk)
df = rename(nd,1=>"modis",2=>"wasim")
model = lm(@formula(modis ~ wasim), df)
#gr()
data = df[!,Not(:date)]
rename!(data,1=>"x",2=>"y")
r = lm( @formula(y ~ x), data)
pred = predict(r, data, interval = :confidence, level = 0.95)
pd = @df data Plots.scatter(:x, :y, leg = false)
# sort data on x
pred_s = pred[sortperm(data[!,:x]), : ]
x_s = sort(data[!, :x])
Plots.annotate!([(7,7,text("hey", 14, :left, :top, :green))])
Plots.annotate!([(10,3,"some R²...")])
Plots.plot!(pd, x_s, pred_s.prediction, linewidth = 2,
        ribbon = (pred_s.prediction .- pred_s.lower, pred_s.upper .- pred_s.prediction))

#Plots.annotate!([(7,3,"(7,3)"),(3,7,text("hey", 14, :left, :top, :green))])
#r.model|>DataFrame

df = readdf("tout")

rename!(df,2=>"y",3=>"x1",4=>"x2")

using GLM
model = glm(@formula(y ~ x1 + x2), df, Poisson(), LogLink())
fitted = predict(model) # a vector of fitted values
residuals = predict(model, residuals=true) # a vector of residuals
df_fitted = hcat(df, fitted=fitted)
df_residuals = hcat(df, residuals=residuals)

backends()
#import Pkg; Pkg.add("PyPlot")
xr = pyimport("xarray")
ENV["PYTHON"]="/usr/bin/python3.8"
pd = pyimport("pandas")
#xr = pyimport("xarray")
pyplot()

rename!(df,1=>"x1",2=>"x2")
model = sm.OLS(df[!,"x1"], df[!,"x2"])
res = model.fit()
res.summary()

#import hydroeval as he
# Conda.add("hydroeval") #gibts nich tin conda
# he = pyimport("hydroeval.api")
Conda.add("xarray")
Conda.add("netcdf4")
xr = pyimport("xarray")
vv=ct("sb")
xr.open_dataset(vv[1])
#xr.open_dataset("sb05rcm_1000")
import Pkg
Pkg.add("PyPlot")
#using Pyplot
#pyplot()
using PyCall
ENV["PYTHON"]=""
Pkg.build("PyCall")

#thats the mpl version
PyPlot.version  

xr = pyimport("xarray")
"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m6"|>cd
ds = xr.open_dataset("tempfab.2012.nc")
##interactive plotting
pygui(true)
ds["t"].plot()|>only|>display

display(ds["tempfab"].plot())

using Conda
Conda.add("rasterio")
# pygui(false) #dann kommt nix.
pygui(true)
ds = xr.open_rasterio("tempfab.2012.nc")
display(ds[0].plot())
mm = ds[0].where(ds[0]>0)
display(mm.plot())
#ds.plot()


ct("stack")
ds = xr.open_dataset("/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m6/hgeofab__stack.2012.nc")
display(ds.sel(t=2,drop=true).plot())

ds.sel(t=2,drop=true)["hgeofabstack"].plot()|>display
vrs = list(ds.data_vars)|>only

ds[vrs].sel(t=3).where(ds[vrs]>0).plot()|>display
ds[vrs].where(ds[vrs]>0).sel(t=2).plot()|>display

##interpo geht so ... lieber mit gr backend
rr=readras(r"hgeo")
#Plots.plot(rr[Band(3)])
vgjl("Dim")
Plots.plot(rr[t=3])

gr()

rr[Dim{:t}(Rasters.Between(4,end))] |>contourf

vgpy("viridis")
ds[vrs].where(ds[vrs]>0).sel(t=2).transpose().plot(cmap="jet")|>display

##mhm stuff
cd("/mnt/c/users/chs72fw/documents/EFRE_GIS/Hydrologie/mhm/mhm-wsl/mhm/pre-proc/test_cut_mhm_input/333.0/morph")
dem=Raster("dem.asc")
Plots.plot(dem)

lk="/mnt/c/users/chs72fw/documents/EFRE_GIS/Hydrologie/mhm/mhm-wsl/mhm/pre-proc/test_cut_mhm_input/398.0/gauge/00398.txt"


function readmhm(x::AbstractString)
    "--- main reader ---"
    ms=["-9999","lin","log"]
    #x=lk
    df::DataFrame = CSV.read(x,DataFrame,
    header=false,
    missingstring = ms,
    types =  Dict(6=>Float64,
    1=>Int64,
    2=>Int64,
    3=>Int64,
    4=>Int64,
    5=>Int64), 
    skipto = 6,
    delim=" ",
    ignorerepeated=true,
    silencewarnings=false,
    normalizenames=true)
    DataFrames.metadata!(df, "filename", x, style=:note);
    #DataFrames.metadata!(df, "basename", basename(x), style=:note);
    #split(basename(x),".")[1]
    #    drop=(i, nm) -> i == 4) |> dropmissing
    df = df[!,Cols(1:3,end)]      #subset to dayres only
    newnames = ["YY","MM","DD",split(basename(x),".")[1]]
    rename!(df, newnames)
    df.YY=map(x ->Int(x),df.YY);
    df.MM=map(x ->Int(x),df.MM);
    df.DD=map(x ->Int(x),df.DD);
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]    
end

df = readmhm(lk)
dfp(df)

#missingstring=ms,
df = CSV.read(x,DataFrame,
    header=false,
    types =  Dict(6=>Float64,
    1=>Int64,
    2=>Int64,
    3=>Int64,
    4=>Int64,
    5=>Int64), 
    skipto = 6,
    delim=" ",
    ignorerepeated=true,
    silencewarnings=false,
    normalizenames=true)


res = Raster("/mnt/c/users/chs72fw/documents/EFRE_GIS/Hydrologie/mhm/mhm-wsl/mhm/test_basin/output_b1/mHM_Fluxes_States.nc";name="QD")
res[Ti=12]|>contourf
#Plots.plot()
rn="/mnt/c/users/chs72fw/documents/EFRE_GIS/Hydrologie/mhm/mhm-wsl/mhm/test_basin/output_b1/mHM_Fluxes_States.nc"
#Raster(rn;name="recharge")|>Plots.plot #2d
rc = Raster(rn;name="recharge")
#rc[:,1,:]
#tim=map(parent, dims(rc, (Ti)))
tim = dims(rc, (Ti))|>collect
Plots.plot(rc)
#nd=map(parent,rc)
vgjl("parent")
#dd = rc.data|>collect|>DataFrame
#nd=DataFrame(date=tim,dat=rc.data)

ds = xr.open_dataset(rn)
list(ds.data_vars)|>show
vrs = "recharge"
ds|>show

#ds[vrs].isel(time=3).where(ds[vrs]>0).plot()|>display
ds[vrs].plot()|>display

hombr()
mg = ct("qbas")
lplot(mg[2])
lplot(mg[1])

"/mnt/d/temp/saale/out_smf200/cdx_v3/qbassmf200.cordex.v3.1973"|>lplot
"/mnt/d/temp/saale/out_smf200/cdx_v3/qbassmf200.cordex.v3.1973"|>dfp
dfp(r"/mnt/d/temp/saale/out_smf200/cdx_v3/qges")

"/mnt/d/temp/saale/out_smf200/cdx_v4/"|>cd
dfp(r"qgessmf2")
"/mnt/d/temp/saale/out_smf200/cdx_v4/prec_station"|>cd
dfp(r"qges")

"/mnt/d/temp/saale/control/"|>cd
vg("prec_station","ctl")
ll()

"/mnt/d/temp/franken_s580_ugm/spin_09_02"|>cd
dfp(r"qg")
dfp(r"qbas")
lplot(r"qbas")
cdb()

vio(r"qges")
vgjl("shap")
pt="/mnt/d/temp/franken_s580_ugm/All_HydrologicResponseUnits.shp"
using Shapefile
shp = Shapefile.Handle(pt)
rr = readras("/mnt/d/temp/franken_s580_ugm/All_HydrologicResponseUnits.nc")
cut = trim(mask(rr; with=shp.shapes[1:11]))
Plots.plot(cut)