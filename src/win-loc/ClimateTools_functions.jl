using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)
#####################
# My own functions
if Sys.isapple()
    platform = "osx"
    const homejl = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
    const mybash = "/Users/apfel/.bash_aliases"
    src_path = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
elseif Sys.iswindows()
    platform = "windows"
    src_path = "C:\\Users\\Public\\Documents\\Python_Scripts\\julia"
    macro wasim() pt="C:\\Users\\chs72fw\\.julia\\dev\\WaSiM\\src\\wa.jl";include(pt);end
else
    platform = "unix"
    winpt = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
    pcld = "~/pCloud Drive/Stuff/Python_Scripts/julia"
    src_path = isdir(winpt) ? winpt : pcld
    println("sourcepath is $src_path")
    if isdir(winpt)
        macro wasim() pt="/mnt/c/Users/chs72fw/.julia/dev/WaSiM/src/wa.jl";include(pt);end
    end
end  

"""
load ClimateTools functions
#pt=/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl
"""
function loadcdo()
    #include(joinpath(src_path, "/win/ClimateTools_functions.jl"))
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end

#loadcdo()

"""
ncremapbil 8 Threads.
"""
function ncremap(input_file, remapfile, output_file, startdate, enddate;
    oldname::Union{Nothing, String}=nothing, newname::Union{Nothing, String}=nothing)
    try
    # Construct the CDO command for cropping
    # -setcalendar,365_day 
    # -P <nthreads>  Set the number of OpenMP threads
    cdo_command = `cdo -P 8 -seldate,$startdate,$enddate -remapbil,$remapfile`

    if oldname !== nothing
    cdo_command = `$(cdo_command) -chname,$oldname,$newname`
    end

    cdo_command = `$(cdo_command) $input_file $output_file`

    # Run the CDO command using run
    run(cdo_command)

    println("NetCDF file successfully processed. Output saved to $output_file")

    catch e
    println("Error: $e")
    println("Failed to remap $output_file.")
    end
end

"""
remapnn 8 Threads.
"""
function ncremapnn(input_file, remapfile, output_file, startdate, enddate;
    oldname::Union{Nothing, String}=nothing, newname::Union{Nothing, String}=nothing)
    try
    # Construct the CDO command for cropping
    # -setcalendar,365_day 
    # -P <nthreads>  Set the number of OpenMP threads
    cdo_command = `cdo -P 8 -seldate,$startdate,$enddate -remapnn,$remapfile`

    if oldname !== nothing
    cdo_command = `$(cdo_command) -chname,$oldname,$newname`
    end

    cdo_command = `$(cdo_command) $input_file $output_file`

    # Run the CDO command using run
    run(cdo_command)

    println("NetCDF file successfully processed. Output saved to $output_file")

    catch e
    println("Error: $e")
    println("Failed to remap $output_file.")
    end
end

"""
remapcon 8 Threads.
"""
function ncremapcon(input_file, remapfile, output_file, startdate, enddate;
    oldname::Union{Nothing, String}=nothing, newname::Union{Nothing, String}=nothing)
    try
    # Construct the CDO command for cropping
    # -setcalendar,365_day 
    # -P <nthreads>  Set the number of OpenMP threads
    cdo_command = `cdo -P 8 -seldate,$startdate,$enddate -remapcon,$remapfile`

    if oldname !== nothing
    cdo_command = `$(cdo_command) -chname,$oldname,$newname`
    end

    cdo_command = `$(cdo_command) $input_file $output_file`

    # Run the CDO command using run
    run(cdo_command)

    println("NetCDF file successfully processed. Output saved to $output_file")

    catch e
    println("Error: $e")
    println("Failed to remap $output_file.")
    end
end


"""
remapdis 8 Threads.
Distance-weighted average remapping
Performs a distance-weighted average remapping of the nearest neighbors value on all
input ﬁelds. The default number of nearest neighbors is 4.
"""
function ncremapdis(input_file, remapfile, output_file, startdate, enddate;
    oldname::Union{Nothing, String}=nothing, newname::Union{Nothing, String}=nothing)
    try
    # Construct the CDO command for cropping
    # -setcalendar,365_day 
    # -P <nthreads>  Set the number of OpenMP threads
    cdo_command = `cdo -P 8 -seldate,$startdate,$enddate -remapdis,$remapfile`

    if oldname !== nothing
    cdo_command = `$(cdo_command) -chname,$oldname,$newname`
    end

    cdo_command = `$(cdo_command) $input_file $output_file`

    # Run the CDO command using run
    run(cdo_command)

    println("NetCDF file successfully processed. Output saved to $output_file")

    catch e
    println("Error: $e")
    println("Failed to remap $output_file.")
    end
end

using PyCall
@pyimport xarray as xr
@pyimport geopandas as gpd
@pyimport matplotlib.pyplot as plt

# @pyimport pyproj
function plot_mean_with_shapefile(pt::AbstractString, 
    myvar::AbstractString, shapefile_path::AbstractString,
    mycrs::AbstractString="EPSG:4326"; kwargs...)

    # Load the shapefile
    shp = gpd.read_file(shapefile_path)
    # Reproject the shapefile to EPSG:4326
    #project = pyproj.Transformer.from_crs("EPSG:25832", "EPSG:4326", always_xy=true)
    shp = shp.to_crs(crs=mycrs)
    
    ob = xr.open_dataset(pt)
    ds = ob.where(ob[myvar] > 0, drop=true)
    plt.rc("font", family="serif", serif=["cmr10"])
    # Plot the mean of the variable
    ds[myvar].mean("time").plot()
    # Plot the outline of the shapefile
    shp.boundary.plot(ax=plt.gca(), 
        color="red", linewidth=1.75, linestyle="dotted")

    plt.show()
end

# pt="simh_bil.nc"
# pt=@gl "obs"
# myvar="rh"
# #shapefile_path="/mnt/d/Wasim/regio/rcm200/ezg.shp"
# shapefile_path="/mnt/d/Wasim/regio/rcm200/v4/ezg.shp"
# filesize(shapefile_path)
# plot_mean_with_shapefile(pt, myvar, shapefile_path;titlestr="obs rh")

# "D:/remo/cordex/eobs/v28/proj/"|>cd
# pt="sfcWind-obs.nc"
# myvar="sfcWind"
# shapefile_path="D:/Wasim/regio/rcm200/v14/catchment_v14.shp"
# Main.pyjl.rglob(".shp","D:/Wasim/regio/rcm200/v0/")
# shapefile_path="D:/Wasim/regio/rcm200/v0/basins-v0.shp"
# plot_mean_with_shapefile(pt, myvar, 
#     shapefile_path, 
#     mycrs="EPSG:25832";
#     titlestr="obs sfcWind")

# shp = gpd.read_file(shapefile_path)
# shp = shp.to_crs(crs="EPSG:25832")

# fn="D:/remo/cordex/eobs/v28/proj/pre-obs.nc"
# ds = xr.open_dataset(fn)
# myvar = "pre"
# ds[myvar].mean("time").plot()
# shp.boundary.plot(ax=plt.gca(), color="red", 
#     linewidth=1.75, linestyle="dotted")

# ds[myvar].mean(["x","y"]).groupby("time.year").sum().plot()
# ds.close()
# ###das ist auch sehr komisch.
# "D:/Wasim/Tanalys/DEM/brend_fab/out/c9/c9-eobs4/"
# "D:/Wasim/Tanalys/DEM/brend_fab/out/c10/s2/"
# ## es waren xc und yc float in der gridlist....
# ##nc scaling_factor checken!!!!
# fn="D:/remo/cordex/eobs/v28/biascorr/sfcWind_cor.nc"
# ds = xr.open_dataset(fn)
# myvar = "sfcWind"
# #Dimensions without coordinates: t, x, y
# ds[myvar].mean("t").plot()
# #ds.t as Array{Dates.DateTime,1} ?
# @pyimport datetime
# @pyimport numpy as np
# pyslice = pybuiltin("slice")
# pyint = pybuiltin("int")
# #pyint(12354.654)
# dt = [datetime.datetime(1990, 1, 1) + datetime.timedelta(days=pyint(i)) for i in ds.t]
# ds.close()
# fn="D:/remo/qm/corgrids/jlcor/pre-cor_extremes_raw.nc"
# ds = xr.open_dataset(fn)
# myvar = "pre"

# shp = shp.to_crs(crs="EPSG:4326")
# collect(ds.keys())
# collect(ds.dims)
# ds[myvar].mean("time").plot()
# shp.boundary.plot(ax=plt.gca(), color="red", 
#     linewidth=1.75, linestyle="dotted")
# ds.close()

using PyCall
function pysetup()
    nm = "/mnt/c/Users/Public/Documents/Python_Scripts/myfuncs-jl.py"
    py"""
    with open($nm) as f:
        exec(f.read())
    """
end

println("formatting: alt + shift + f")
println("new vscode window: alt + shift + n")
println("win/ClimateTools_functions.jl loaded")
println("used Threads: ", Threads.nthreads())

#pysetup()

# cd("/mnt/d/Wasim/regio/out/rc200/x22/f18-spin")
# ssup()
# k=@gl "qg"

# py"""
# df=waread3($k)
# vio(df)
# """

# py"""heat($k)"""
# v = glob(r"qoutjl$")

# py"""heat($v[2])"""
# #Conda.pip("install","hydroeval")
# py"""theplot($v[2])"""
# py"""theplot($v[3])"""
# v[3]|>cb
# py"""ftp('Wolf','Q')"""
# cd("/mnt/d/Wasim/regio/out/rc200/x22/eobs/loc/")
# pysetup()
# py"""ftp('Wolf','M')"""
# py"""ftp('Wolf','Y')"""