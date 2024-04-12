"D:/remo/qm/qm/"|>cd

#before activating the env, one could try...
# wapt=raw"C:\Users\chs72fw\.julia\dev\WaSiM\src\wa.jl"
# @time include(wapt)

using Pkg
Pkg.activate(".")
Pkg.status()
Pkg.update()
using Base.Threads
nthreads()
#ENV["PROJ_LIB"] = "c:/OSGeo4W64/bin"
ENV["PROJ_LIB"]
using Conda
Conda.ROOTENV

using ClimatePlots
using ClimateTools
#https://juliaclimate.github.io/ClimateTools.jl/stable/datasets/

wrkdir="D:/remo/qm/prec"
wrkdir="/mnt/d/remo/qm/prec"
cd(wrkdir)
ssup()
ls()
latx()
methods(qqmap)

ob = load("obsh.nc","pre")
osp = load("obsp.nc","pre")
rf = load("simh.nc","pre")
fut = load("simp.nc","pre")

const cm = ClimatePlots


cm.plot(ob,filename="tst.png")
using PyCall
PyCall.pygui()
mpl=pyimport("matplotlib")
pygui(true)
cm.plot(ob)
cm.plot(annualsum(ob))
cm.plot(annualsum(rf))


@vv "savefig"
@vv "backend"
@vv "module"

out = qqmap(ob,rf,fut;) #method=
#method::String = "Additive" (default) or "Multiplicative". Additive is used for most climate variables. Multiplicative is     
#usually bounded variables such as precipitation and humidity.
-out = qqmap(ob,rf,fut;method="Multiplicative",interp=Quadratic())
#->nested task error: only Linear, Constant, and NoInterp supported, got Gridded(Quadratic(Line(OnGrid())))
out = qqmap(ob,rf,fut;method="Multiplicative",interp=Linear())
ClimateTools.write(out,"pre-mult.nc")


ClimatePlots.contourf(out,region="EU",filename="tst.png");
cm.plot(out,filename="tst.png");
z=lat()
wa.irfan(z)

old=ClimateTools.annualsum(fut) 
ind=ClimateTools.annualsum(out) 
cm.plot(ind,filename="pas.png");
cm.plot(old,filename="pas-nr.png");
lat()|>wa.irfan


#https://juliaclimate.github.io/ClimateTools.jl/stable/biascorrection/
#Quantile-quantile mapping (Themeßl et al. 2012, Grenier et al. 2015) is 
#provided with ClimateTools.jl through the function qqmap.
qqmap(obs::ClimGrid, ref::ClimGrid, 
fut::ClimGrid; method::String="Additive", 
detrend::Bool=true, window::Int=15, 
rankn::Int=50, thresnan::Float64=0.1, 
keep_original::Bool=false, interp = Linear(), extrap = Flat())


df = wa.ncdf("simp.nc")
dfm(df)
ClimateTools.write(out,"pre-add.nc")
df2 = wa.ncdf("pre-add.nc")
wa.dfm(df2;fun=yrsum)
df.date .= DateTime.(string.(df.date))
df2.date .= DateTime.(string.(df2.date))
dfp(df)
dfp!(df2)

ma = mall(df,df2)
describe(ma)
wa.qplot(ma)

r=Raster("pre-mult.nc")
plot(r[Ti=50])

d3 = ncdf(r)
ma = mall([df,df2,d3])
rename!(ma,1=>"no_corr",2=>"date",3=>"add",4=>"mult")
qplot(select(ma,[1,4]))
qplot(select(ma,[3,4]))
wa.dfm(ma) #<-das wäre eine gute Abbildung für die Präsentation



vgjl("MPLBACKEND")
ENV["MPLBACKEND"] = "TkAgg"
#ENV["MPLBACKEND"] = "Qt5Agg"

px = cm.plot(ind);
ClimatePlots.show(px)

using PyCall
matplotlib=pyimport("matplotlib")
matplotlib.use("TkAgg")
default(show=true) 
ox = ClimatePlots.contourf(ind);
display(ox)
show(ox)
pygui(true)   #very important
cm.plot(out)
@vv "pygui"

###18,09,2023
using ClimatePlots
using ClimateTools
"D:/remo/qm/prec/"|>cd
r=Raster("pre85.nc")
plot(r[Ti=1])
df=wa.ncdf(r)
df.date[1]|>typeof

split(string(df.date[1]),"(")

parse(DateTime,string(df.date[1]))

dfp(df)




using PyCall
plt = pyimport("matplotlib.pyplot")
PyCall.pygui() #Return the current GUI toolkit as a symbol. 
pygui(true)
xr = pyimport("xarray")

simh=xr.open_dataset("simh.nc")
simp=xr.open_dataset("simp.nc")
obsh=xr.open_dataset("obsh.nc")

begin
    fig = plt.figure(figsize=(10,5), dpi=216)
    simh["pre"].groupby("time.dayofyear").mean(...).plot(label="\$P_{sim,h}\$")
    simp["pre"].groupby("time.dayofyear").mean(...).plot(label="\$P_{sim,p}\$")
    obsh["pre"].groupby("time.dayofyear").mean(...).plot(label="\$P_{obs,h}\$")
    plt.title("Historical modeled and observed Precipitation and predicted values")
    plt.xlim(0, 365)
    plt.gca().grid(alpha=0.3)
    plt.legend()
end


x, y1, y2 = rand(10),rand(10),rand(10)

np = pyimport("numpy")
# plot
fig, ax = plt.subplots()
ax.fill_between(x, y1, y2, alpha=.5, linewidth=0)
ax.plot(x, (y1 + y2)/2, linewidth=2)
ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
       ylim=(0, 8), yticks=np.arange(1, 8))

plt.show()



##########load from a git repo
pt=raw"C:\Users\Public\Documents\Python_Scripts\julia\win\ClimateTools_functions.jl"
include(pt)

using GeoStats
modelinterp = regrid(fut,ob)

nx, ny = 20, 10
vX = [x for x in range(0,10,length=nx), j in 1:ny]
vY = sin.(vX) .+ [0.5j for i in 1:nx, j in 1:ny]

g = GeoStats.StructuredGrid(vX, vY)

plot(g)



#### daily mean temperature TG######
cdb()
fdi()
wrkdir="D:/remo/qm/tmp"
cd(wrkdir)

using ClimatePlots
using ClimateTools
const tt = ClimateTools
const cm = ClimatePlots

using PyCall
pygui(true)
PyCall.pygui()
mpl=pyimport("matplotlib")

#first regrid...
tas = load("D:/remo/cordex/wgs/proj_tas_hist+rcp85_wgs84.nc","tas")
op()
nconly(".")
tg="D:/remo/cordex/eobs/tg_enc_nck_crop.nc"
tg = load(tg,"tg")

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


#doesnt work
pt=raw"C:\Users\Public\Documents\Python_Scripts\julia\win\ClimateTools_functions.jl"
include(pt)
using GeoStats
regrid(ob,rf)

##so i need to save the data first and regrid it via cdo...
#NCDatasets.renameDim is a Function.
#renameDim(ds::NCDataset, oldname::Union{AbstractString, Symbol}, newname::Union{AbstractString, Symbol})
#renameDim(ob,:tg,:tas)
pwd()
ClimateTools.write(ob,"obsh.nc") #also errors :(
tt.renameDim(ob,:tg=>:tas)


out = qqmap(ob,rf,fut;method="Additive") 


#method::String = "Additive" (default) or "Multiplicative". Additive is used for most climate variables. Multiplicative is     
#usually bounded variables such as precipitation and humidity.
-out = qqmap(ob,rf,fut;method="Multiplicative",interp=Quadratic())
#->nested task error: only Linear, Constant, and NoInterp supported, got Gridded(Quadratic(Line(OnGrid())))
out = qqmap(ob,rf,fut;method="Multiplicative",interp=Linear())
ClimateTools.write(out,"pre-mult.nc")


ClimatePlots.contourf(out,region="EU",filename="tst.png");
cm.plot(out,filename="tst.png");
z=lat()
wa.irfan(z)

old=ClimateTools.annualsum(fut) 
ind=ClimateTools.annualsum(out) 
cm.plot(ind,filename="pas.png");
cm.plot(old,filename="pas-nr.png");
lat()|>wa.irfan



"D:/remo/qm/prec/"|>cd
glob("cor")
r=Raster("pre-cor.nc")
dx=wa.ncdf(r)

plot(r)

df=wa.pyread(glob("cor")|>last)
dfm(df)
dfm(df;fun=yrsum,mode=:line)

#:markersize ∈ plotattr(:Series)
#map(col->tryparse(Float64, col),)


"D:/remo/qm/"|>cd
tree()
"D:/remo/qm/qm"|>cd
using Pkg
Pkg.activate(".")

"D:/remo/qm/rsds"|>cd
ssup()
nconly(".")
using ClimateTools
ob = load("obs.nc","rsds")         #<-THIS errors. wrong datatype
simh=load("simh.nc","rsds")
simp=load("simp.nc","rsds")


using Conda
Conda.CONDA_EXE
Conda.add("cartopy")
#Conda.rm("cartopy")
#Conda.pip("install","cartopy")
using ClimatePlots
using ClimateTools
using PyCall
pygui(true)
const cm = ClimatePlots
#pyimport_conda("cartopy", "scipy")
#this is a faulty dataset
#lk=raw"D:\remo\cordex\eobs\qq_ens_spread_0.1deg_reg_v26.0e.nc"
#daily mean global radiation QQ.
lk=raw"D:\remo\cordex\eobs\qq_ens_mean_0.1deg_reg_v27.0e.nc"
#lk = "D:/remo/cordex/eobs/proj/qq_ens_spread_utm.nc"    

#cdo -v -P 8 remapbil,$gr $lv qq_ens_mean_v27_crop.nc
#remapcon: YAC first order conservative weights from lonlat (48x20) to gaussian (128x64) grid
#cdo -P 6 remapcon,n32  qq_ens_mean_v27_crop.nc qq_con.nc #nope
lk = "D:/remo/cordex/eobs/qq_ens_mean_v27_crop.nc"    
#q2 = wa.ncdf(lk)
q2 = load(lk,"qq")
#"qq_ens_xr_crop.nc"|>rm
cd(dirname(lk))
pwc()
#lk = "D:/remo/cordex/eobs/qq_con.nc"    
lk = "D:/remo/cordex/eobs/qq_ens_xr_crop.nc"
q2 = load(lk,"qq";)
cm.plot(q2|>annualmean)

xr=pyimport("xarray")
ds=xr.open_dataset(lk)
ds['qq'] = ds['qq'].astype('float64')
ds.to_netcdf("qq_ens_xr_crop.nc") 
ds.close() 

pwc()
####
obslnk="qq_ens_xr_crop.nc"
simlnk="/mnt/d/remo/cordex/wgs/proj_rsds_hist+rcp85_wgs84.nc"
cdo -v chname,qq,rsds -setcalendar,365_day -seldate,1970-01-01,2017-12-31 $obslnk obs.nc
cdo -v -setcalendar,365_day -seldate,1970-01-01,2017-12-31 -remapbil,obs.nc $simlnk simh.nc
cdo -setcalendar,365_day -seldate,1970-01-01,2100-01-01 -remapbil,obs.nc $simlnk simp.nc
mrz
myvar="rsds"
wrkdir=pwd()
obs = load(joinpath(wrkdir,"obs.nc") ,myvar) 
println([maximum(obs), minimum(obs)])
simh =load(joinpath(wrkdir,"simh.nc"),myvar)
simp =load(joinpath(wrkdir,"simp.nc"),myvar)
cm.plot(annualmean(obs))
cm.plot(annualmean(simp))

#C:\Users\chs72fw\.julia\packages\ClimateTools\cUGFX\src\biascorrect.jl
#out = qqmap(obs,simh,simp;method="Additive",thresnan=2.99)
out = qqmap(obs,simh,simp;method="Multiplicative",thresnan=0.1)

cm.plot(annualmax(out))
cm.plot(annualmax(obs))
cm.plot(annualmax(simp))
#contourf(annualsum(out))
cm.plot(annualmax(out),filename="$myvar-corr.png")
contourf(annualsum(out),filename="$myvar-sum.png")

ClimateTools.write(out,"$myvar-cor.nc")


x="D:/remo/cordex/eobs/rsds-cor.nc"
ds=xr.open_dataset(x)  
ds["rsds"].isel(time=234).plot() 

##for wasimgrids
n="D:/temp/saale/out_30m/v1-new/Soil_Temperature_Stack.nc"
o=xr.load_dataset(n)
plt=pyimport("matplotlib.pyplot")
o.isel(t=2).transpose().to_array().plot(robust=true,cmap=plt.cm.RdYlBu_r) 

#o.where(o.groupby("t").mean("t"))

m = o.keys() |>collect |>last
o[m].transpose().plot(col="t",col_wrap=4,robust=true,cmap=plt.cm.RdYlBu_r)

da = o.where(o.t > 17, drop=true).transpose()    #subset layers.
da[m].plot(col="t",col_wrap=4,robust=true,
    cmap=plt.cm.RdYlBu_r,cbar_kwargs=Dict(
    "orientation"=> "horizontal",
    "shrink"=> 0.8,
    "aspect"=> 40,
    "pad"=> 0.1))


ds.where(ds["time.year"]==2014,drop=true).mean("time").qq.plot()
ds.where(ds.groupby("time.month").mean("time")).qq.plot()

a="D:/temp/saale/out_smf200/v7/thetsmf180_stack.2017.nc"
ad=xr.open_dataset(a)
m = ad.keys()|>collect|>last    
ad[m].where(ad[m]>0).transpose().plot(col="t",col_wrap=4,robust=true,cmap=plt.cm.RdYlBu_r)

function xrfacets(a::AbstractString)
    """
    uses xarray to plot a 4xn grid of facets.
    """
    xr  = pyimport("xarray")
    plt = pyimport("matplotlib.pyplot")
    ad = xr.open_dataset(a)
    m = ad.keys()|>collect|>last    
    p1 = ad[m].where(ad[m]>0).transpose().plot(col="t",col_wrap=4,
        robust=true,
        cmap=plt.cm.RdYlBu_r);
    return p1.fig #now we can see the plot inside vscode.
end


mon = ds.where(ds.groupby("time.month").mean("time"))   
#mon = mon.drop_dims("time")
mean_mon = mon.mean(dim="month")
fg = mean_mon.plot(
    col="month",
    col_wrap=4,
    robust=true,
    cmap=plt.cm.RdYlBu_r,
    cbar_kwargs=Dict(
        "orientation"=> "horizontal",
        "shrink"=> 0.8,
        "aspect"=> 40,
        "pad"=> 0.1))
##Only 1d and 2d plots are supported for facets in xarray.


mon.qq.isel(longitude=1).plot(
    x="month",    robust=true,    cmap=plt.cm.RdYlBu_r,
    cbar_kwargs=Dict(
        "orientation"=> "horizontal",
        "shrink"=> 0.8,
        "aspect"=> 40,
        "pad"=> 0.1))


using RCall
tr=rimport("terra")
ra = tr.rast("$myvar-cor.nc")
@rput ra
R"
library(terra)
plot(ra[[5]])
"

R"""
# Generate data for the sine and cosine waves
x <- seq(0, 2 * pi, length.out = 100)
sine_wave <- sin(x)
cosine_wave <- cos(x)
# Create the plot
plot(x, sine_wave, type = "l", col = "blue", xlab = "X-axis", ylab = "Y-axis", 
     main = "Sine and Cosine Waves", ylim = c(-1, 1))
lines(x, cosine_wave, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Sine Wave", "Cosine Wave"), 
        col = c("blue", "red"),
        lty = 1, cex = 1.1)

# Add a grid
grid()

"""


cm.bar()
as = annualsum(ob)
ar = getindex(as,1)
typeof(ar)
data = ar.data

ti = ar.axes|>last|>collect
v = []
for i in 1:length(ti)
    d = data[:,:,i]|>collect|>vec
    push!(v,d)
end
#DataFrame(data[:,:,1],:auto)
df = DataFrame(v,:auto)|>permutedims
df = hcat(df, parent(ti), makeunique=true)
rename!(df,ncol(df)=>:year)

cm.bar(df.year,df.x3)
rmean = DataFrames.combine(df, names(df)[1:end-1] .=> mean)
rmean = DataFrames.transform(df, names(df)[1:end-1] .=> mean)
#mean(df,dims=2)
# Assuming df is your DataFrame
df[!,:mean] = mean.(eachrow(df[!,Not(:year)]))
cm.bar(df.year,df.mean)



data_vectors = collect(row -> [x for x in row], data[:,:,1])
data_vectors = collect(row -> [x for x in row], data)
df = DataFrame(data_vectors)

run(`/mnt/c/Windows/explorer.exe .`)
run(`cmd.exe /c start .`)

