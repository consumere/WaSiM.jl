
#jl --threads auto -q
using Base.Threads
nthreads()
#Pkg.add("ClimatePlots")
using ClimatePlots
#import Pkg; Pkg.Registry.update(); Pkg.activate(); Pkg.instantiate();
#ENV["PYTHON"]="C:\\Users\\chs72fw\\Miniconda3\\python.exe"
#ENV["PYTHON"]=""
# using Conda
# Conda.create(raw"C:\Users\chs72fw\.julia\conda")
# Conda.add("cmocean")
#Conda.add("numpy")
using ClimateTools
#https://juliaclimate.github.io/ClimateTools.jl/stable/datasets/
#C = load(filename::String, vari::String; poly::Array, data_units::String, start_date::Tuple, end_date::Tuple, dimension::Bool=true)

#  Temporal subsetting can be done by providing start_date and end-date Tuples of length 1 (year), length 3 (year, month, day) or 6 (hour, minute, second).
#C = load(fn,"hu";data_units="%",start_date=Tuple("2014"),end_date=Tuple("2015"))
parsetuple(s::AbstractString) = Tuple(parse.(Float64, split(s, ',')))
parsetuple("234,456,654")


"""
    C = regrid(A::ClimGrid, B::ClimGrid; solver=Kriging(), min=[], max=[])

Interpolate `ClimGrid` A onto the lon-lat grid of `ClimGrid` B, using the GeoStats.jl methods.

Min and max optional keyword are used to constraint the results of the interpolation. For example, interpolating bounded fields can lead to unrealilstic values, such as negative precipitation. In that case, one would use min=0.0 to convert negative precipitation to 0.0.

"""
function rgrid(A::ClimGrid, B::ClimGrid; solver=GeoStats.Kriging(), min=[], max=[])

    # ---------------------------------------
    # Get lat-lon information from ClimGrid B
    londest, latdest = ClimateTools.getgrids(B)

    # Get lat-lon information from ClimGrid A
    lonorig, latorig = ClimateTools.getgrids(A)
    points = hcat(lonorig[:], latorig[:])

    # -----------------------------------------
    # Get initial data and time from ClimGrid A
    dataorig = A[1].data
    timeorig = get_timevec(A) #[1][Axis{:time}][:] # the function will need to loop over time

    # ---------------------
    # Allocate output Array
    OUT = zeros(Float64, (size(B.data, 1), size(B.data, 2), length(timeorig)))

    # ------------------------
    # Interpolation
    interp!(OUT, timeorig, dataorig, lonorig, latorig, londest, latdest, solver=solver, msk=B.msk)

    if !isempty(min)
        OUT[OUT .<= min] .= min
    end

    if !isempty(max)
        OUT[OUT .>= max] .= max
    end

    # -----------------------
    # Construct AxisArrays and ClimGrid struct from array OUT
    latsymbol = Symbol(B.dimension_dict["lat"])
    lonsymbol = Symbol(B.dimension_dict["lon"])
    dataOut = AxisArray(OUT, Axis{lonsymbol}(B[1][Axis{lonsymbol}][:]), Axis{latsymbol}(B[1][Axis{latsymbol}][:]), Axis{:time}(timeorig))

    C = ClimateTools.ClimGrid(dataOut, longrid=B.longrid, latgrid=B.latgrid, msk=B.msk, grid_mapping=B.grid_mapping, dimension_dict=B.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=A.experiment, run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=B.latunits, lonunits=B.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)

end



fn="d:/remo/cordex/REMO11_2011_2015/prec_wgs_day_20110101-20151231.nc"
C = load(fn,"pr";)
using Dates
t = Date(2014, 1, 31)
#https://juliaclimate.github.io/ClimateTools.jl/stable/biascorrection/
#Quantile-quantile mapping (Themeßl et al. 2012, Grenier et al. 2015) is provided with ClimateTools.jl through the function qqmap.
qqmap(obs::ClimGrid, ref::ClimGrid, 
fut::ClimGrid; method::String="Additive", 
detrend::Bool=true, window::Int=15, 
rankn::Int=50, thresnan::Float64=0.1, 
keep_original::Bool=false, interp = Linear(), extrap = Flat())

fobs="d:/remo/cordex/eobs/proj/rr_ens_spread_utm.nc"
obsg = load(fobs,"rr";)
qqmap(obsg, obsg, C)
#@edit ClimateTools.qqmap()
lnk=raw"C:\Users\Public\Documents\Python_Scripts\julia\win\climgrid2023funcs.jl"
include(lnk)
using GeoStats
# modelinterp = regrid(C, obsg)
# modelinterp = rgrid(C, obsg)

fobs="d:/remo/cordex/eobs/proj/rr_ens_spread_utm.nc"
obsg = load(fobs,"rr";)
fn="d:/remo/cordex/eobs/rr2016.nc"
C = load(fn,"rr";)

modelinterp = griddata(C, obsg,min=0.0)        #DAS!
#C.region = "europe"

ClimatePlots.mapclimgrid(C)

outg = annualmax(modelinterp)
outg = annualmax(obsg)
outg.data|>unique

ClimatePlots.plot(outg,filename="test.png")
ClimatePlots.contourf(C,filename="test.png")
ClimatePlots.contourf(modelinterp,filename="test.png")
using Images 
x=Images.load("test.png")           #geht auch so...


@vv "pyplot()"
#using PyPlot
using Plots
pyplot()
pt = ClimatePlots.plot(outg)
display(pt)

ClimatePlots.plot(outg,filename="test.png") #das geht.
ClimatePlots.plot(modelinterp,filename="test.png") #das geht.

#using IJulia
#IJulia.magic("%matplotlib inline")

#C = resample(C::ClimGrid, season::String) # hardcoded seasons -> "DJF", "MAM", "JJA" and "SON"
A = resample(C, "SON") # hardcoded seasons -> "DJF", "MAM", "JJA" and "SON"
ind = annualmax(C)

mapclimgrid(obsg; region = "World")

modelinterp = regrid(model, obs)

tst = prcp1(C)
#tst[:,:,4]|>mapclimgrid
tst|>mapclimgrid
mapclimgrid(C; region = "Europe")


# using Pkg
# # specify the package name, repository URL, and version
# pkg = "ClimateTools"
# url = "https://github.com/JuliaClimate/ClimateTools.jl.git"
# #version = "v0.22.0"
# version = "0.22"
# # add the package at the specified version
# Pkg.add(PackageSpec(name=pkg, url=url, version=version))
# add ClimateTools@0.22

x="D:/temp/saale/in_mf/smf30/smf.maxpond"
z=load(x)
