
#jl --threads auto -q
using Base.Threads
nthreads()
ENV["PROJ_LIB"] = "c:/OSGeo4W64/bin"
using Conda
Conda.ROOTENV

using ClimatePlots
using ClimateTools
#https://juliaclimate.github.io/ClimateTools.jl/stable/datasets/

wrkdir="D:/remo/qm/"
cd(wrkdir)
ssup()
glob|>methods
readdir(".")
latx()

glob("sim")

ob = load("obsh.nc","pre")
rf = load("simh.nc","pre")
fut = load("simp.nc","pre")

@vv "geht"
@vv "savefig"
@vv "backend"
@vv "module"

out = qqmap(ob,rf,fut;) #method=
#using Plots
using PyPlot
pyplot()

ClimatePlots.contourf(out,region="EU",filename="tst.png");
ClimatePlots.contourf(out,region="Gr",filename="tst.png");
latx()
import Images
im=Images.load("tst.png")
ind=ClimateTools.annualsum(out) 
x=load("tas_delta.nc","tas")

fn="tasmax_A2.Commit_1.PCM1.atmd.2000-01-01_cat_2099-12-31.nc"
x=load(fn,"tasmax")
ind=ClimateTools.annualmax(x) 
ClimatePlots.plot(ind,filename="tst.png")
latx()
im=Images.load("tst.png")

#ClimatePlots.contourf(rf-fut,"tst2.png")
x=annualmax(fut)
ClimatePlots.plot(x,filename="tst2.png")
im=Images.load("tst2.png")
#C = ClimateTools.merge(ob,rf)
ClimatePlots.plot(annualmax(ob),filename="tst2.png")

# #display(p1)
# show(p1)
# outname="tst.png"
# Pyplots.savefig(p1,outname,width=1600,height=800)

# using PyPlot
# # create a figure and add a colorbar
# fig, ax = subplots()
# cbar = colorbar()
# # save the figure
# savefig("colorbar.png")


#https://juliaclimate.github.io/ClimateTools.jl/stable/biascorrection/
#Quantile-quantile mapping (Themeßl et al. 2012, Grenier et al. 2015) is 
#provided with ClimateTools.jl through the function qqmap.
qqmap(obs::ClimGrid, ref::ClimGrid, 
fut::ClimGrid; method::String="Additive", 
detrend::Bool=true, window::Int=15, 
rankn::Int=50, thresnan::Float64=0.1, 
keep_original::Bool=false, interp = Linear(), extrap = Flat())




ENV["MPLBACKEND"] = "TkAgg"
ENV["MPLBACKEND"] = "Qt5Agg"
ClimatePlots.contourf(am)
ClimatePlots.show()
#using Plots
mapclimgrid(C,"World")

ClimatePlots.plot(C)
#using PyPlot
#show()
#C = load(fn,"hu";data_units="%",start_date=parsetuple("2014"),end_date=parsetuple("2015"))
using Dates
t = Date(2014, 1, 31)

t=(2014, 1, 31)

