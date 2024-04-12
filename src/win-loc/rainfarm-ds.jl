#rainfarm calcs

#######

prj="/home/ubu/.julia/rainfarm" 
julia --project=$prj

cd "/mnt/d/Wasim/Tanalys/DEM/brend_fab/in4"
julia --startup-file=no --project=$prj rfarm.jl  
orofile="/mnt/d/Relief_DGMs/FABDEM/dem01.nc"
#"/mnt/d/Relief_DGMs/FABDEM"
#gdalwarp -te 9.38 50.07 10.69 50.58 -tr 0.001 0.001 ufra30_fabdem.tif dem01.tif
#gdalwarp -te 9.38 50.07 10.69 50.58 -tr 0.001 0.001 -of NetCDF ufra30_fabdem.tif dem01.nc
julia --startup-file=no --project=$prj tools/rfweights.jl $orofile l.nc

f="/mnt/d/remo/cordex/eobs/v28/rr_ens_mean_v28.0e_crop.nc"
cdo selyear,2000 rr_ens_mean_v28.0e_rcm.nc ou.nc
f="/mnt/d/remo/cordex/eobs/v28/ou.nc"
julia --startup-file=no --project=$prj tools/rfarm.jl $f
#julia --startup-file=no --project=$prj tools/rfarm.jl $f --region "1/1/1/1"



reffile="pre2000.nc"
#ww = rfweights(orofile, reffile, nf; weightsfn="", varname="", fsmooth=false)
ww = rfweights(orofile, reffile, 2;varname="pre", fsmooth=true)

ENV["PROJ_LIB"]="/home/ubu/.julia/conda/3/share/proj/proj.db"
using ClimateTools,ClimatePlots

function vector_to_square_matrix(vec::Vector{T}) where T
    n = isqrt(length(vec))
    if n * n != length(vec)
        m = ceil(Int, sqrt(length(vec)))
        padding = m^2 - length(vec)
        vec = vcat(vec, fill(0.0, padding))
        n = m
    end
    return reshape(vec, n, n)
end
l = "/mnt/d/Relief_DGMs/FABDEM/dem01.nc"
l= "/mnt/d/remo/cordex/wgs/pr_remap_rr.nc"
l="/mnt/d/Wasim/Tanalys/DEM/brend_fab/in4/l.nc"

v = read_netcdf2d(l)[1] 
vx = replace(v,NaN=>missing)|>skipmissing|>collect|>vector_to_square_matrix 
import ClimatePlots as cpl 
cpl.hist(vx[:,:,1])  
c = rainfarm(vx,1.7,2,1)
cpl.hist(c[:,:,1])  
import ClimateTools as tt
a = tt.load(l,"pre")

function calculate_grid_boundaries(xfirst::Float64, xinc::Float64, yfirst::Float64, yinc::Float64, xsize::Int, ysize::Int, grid_size::Int)
    # Calculate ending longitude and latitude
    xlast = xfirst + xinc * xsize
    ylast = yfirst + yinc * ysize

    # Calculate total range and midpoint
    x_range = xlast - xfirst
    y_range = ylast - yfirst
    x_midpoint = xfirst + x_range / 2
    y_midpoint = yfirst + y_range / 2

    # Calculate starting and ending boundaries for grid_size x grid_size grid
    x_start = x_midpoint - xinc * grid_size / 2
    x_end = x_midpoint + xinc * grid_size / 2
    y_start = y_midpoint - yinc * grid_size / 2
    y_end = y_midpoint + yinc * grid_size / 2

    return x_start, x_end, y_start, y_end
end

xfirst = 9.04986035478427
xinc = 0.0999999995999529
yfirst = 49.0498597218401
yinc = 0.0999999958273361
xsize = 30
ysize = 17
grid_size = 15

x_start, x_end, y_start, y_end = calculate_grid_boundaries(xfirst, xinc, yfirst, yinc, xsize, ysize, grid_size)

println("Starting Longitude: ", x_start)
println("Ending Longitude: ", x_end)
println("Starting Latitude: ", y_start)
println("Ending Latitude: ", y_end)
j = [x_start,x_end,y_start,y_end]
join(round.(j;digits=4),",")|>cb
# cdo sellonlatbox,9.7999,11.2999,49.1499,50.6499 ou.nc ot.nc
# griddes ot.nc

# tools="/mnt/d/Wasim/Tanalys/DEM/brend_fab/in4/tools"
tools=$jlpt/tools 
prj="/home/ubu/.julia/rainfarm" 
f="ot.nc"
julia --startup-file=no --project=$prj $tools/rfarm.jl $f

f="D:/remo/cordex/eobs/v28/rainfarm_0001.nc"
shapefile_path="D:/Wasim/regio/rcm200/v0/basins-v0.shp"
pyjl.plot_mean_with_shapefile(f, shapefile_path;myvar="rr")

xfirst    = 8.09
xinc      = 0.11
yfirst    = 50.93 
xsize     = 45                                                                                                 
ysize     = 18
x_start, x_end, y_start, y_end = calculate_grid_boundaries(xfirst, xinc, yfirst, xinc, xsize, ysize, grid_size)
j = [x_start,x_end,y_start,y_end]
j = [xfirst,x_end,yfirst,y_end]
join(round.(j;digits=4),",")|>cb

(11.39 - 8.09) / 0.11
50.93 - (30*0.11) #weil negativ inc.
cdo -v sellonlatbox,8.09,11.39,47.63,50.93 $fn p.nc
xl = xfirst + (17*.11)
yl = yfirst - (19*.11)
cdo -v sellonlatbox,8.09,9.96,48.84,50.93 $fn p.nc #18x18
f="D:/remo/cordex/wgs/rainfarm_0001.nc"
pyjl.plot_mean_with_shapefile(f, shapefile_path;myvar="pre")
ds=xr.open_dataset(f)
plt.imshow(ds["pre"].mean(dim="time"))
ds["pre"].groupby("time.year").sum().plot()
ds["pre"].groupby("time.dayofyear").mean("lon").mean("lat").plot(label="Reference")
shapefile_path="D:/Wasim/regio/rcm200/y11/catchment-y11-ezk.shp"
pyjl.plot_mean_with_shapefile(f,shapefile_path;myvar="pre")
@pyimport geopandas as gpd
gdf = gpd.read_file(shapefile_path)
gdf = gdf.to_crs("EPSG:4326")
gdf.to_file(joinpath(dirname(shapefile_path),"catchment-y11-4326.geojson"), driver="GeoJSON")  
lpath = joinpath(dirname(shapefile_path),"catchment-y11-4326.geojson")
d = GeoDataFrames.read(lpath)
ed = reduce(GeoDataFrames.union,d.geometry)
cds = GeoInterface.coordinates(ed)
extrema.(cds)
flat_coords = vcat(cds...)
min_x = minimum(x -> x[1], flat_coords)
max_x = maximum(x -> x[1], flat_coords)
min_y = minimum(x -> x[2], flat_coords)
max_y = maximum(x -> x[2], flat_coords)

cd("D:/remo/cordex/wgs")
inf = "proj_pre_hist+rcp85_wgs84.nc"
outf= "pre_rcp85-saale.nc"
run(`wsl.exe -e cdo -v sellonlatbox,$min_x,$max_x,$min_y,$max_y $inf $outf`)
latx() #16mb

f="pre_rcp85-saale.nc"
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")

vgjl("centroid")
import LibGEOS: centroid
centerpoints = [centroid(d.geometry[i]) for i in 1:length(d.geometry)]
cen = centroid(ed)
#LibGEOS.getX(cen)ERROR: MethodError: no method matching getX(::Point)
GeoInterface.coordinates(cen)

function calculate_grid_boundaries_center(center_lon::Float64, center_lat::Float64, xinc::Float64, yinc::Float64, grid_size::Int)
    # Calculate starting and ending boundaries for grid_size x grid_size grid
    x_start = center_lon - xinc * (grid_size - 1) / 2
    x_end = center_lon + xinc * (grid_size - 1) / 2
    y_start = center_lat - yinc * (grid_size - 1) / 2
    y_end = center_lat + yinc * (grid_size - 1) / 2
    return x_start, x_end, y_start, y_end
end

xinc = 0.11
yinc = -0.1099997
grid_size = 14
a,b = GeoInterface.coordinates(cen)
x_start, x_end, y_start, y_end = calculate_grid_boundaries_center(a,b,xinc,yinc,grid_size)
#j = [x_start,x_end,y_start,y_end]
#join(round.(j;digits=6),",")|>cb
x_start = x_start - (0.11 * 1)
y_end = y_end + (0.11 * 4)
inf = "proj_pre_hist+rcp85_wgs84.nc"
outf= "pre_rcp85-saale.nc"
run(`wsl.exe -e cdo -v sellonlatbox,$x_start,$x_end,$y_start,$x_end $inf $outf`)
latx()
run(`wsl.exe -e cdo griddes $outf`)
f="pre_rcp85-saale.nc"
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")

function cb_c(center_lon::Float64, center_lat::Float64, xinc::Float64, yinc::Float64, grid_size::Int)
    # Calculate starting and ending boundaries for grid_size x grid_size grid
    x_start = center_lon - xinc * (grid_size - 1) / 2
    x_end = center_lon + xinc * (grid_size - 1) / 2
    y_start = center_lat - yinc * (grid_size - 1) / 2
    y_end = center_lat + yinc * (grid_size - 1) / 2
    return x_start, x_end, y_start, y_end
end
x_start, x_end, y_start, y_end = cb_c(a,b,xinc,yinc,grid_size)
#y_end = y_end + (0.11 * 4)
y_end = y_end + (0.11 * 1)

inf = "pre_rcp85-saale.nc"
outf = "pre_rcp85-saale2.nc"
run(`wsl.exe -e cdo -v -P 8 sellonlatbox,$x_start,$x_end,$y_start,$y_end $inf $outf`)
f="pre_rcp85-saale3.nc"
run(`wsl.exe -e cdo griddes $f`)
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")

#now wsl rainfarm
pwc()
f="pre_rcp85-saale.nc"
julia --startup-file=no --project=$prj $tools/rfarm.jl $f

f="/mnt/e/qq/prec14.nc"
f = towin(f)
cdof(f)
@pyjl
cd()
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")
rm(f)

import GeoDataFrames,GeoInterface
gdf = gpd.read_file(shapefile_path)
gdf = gdf.to_crs("EPSG:4326")
#gdf.to_file(joinpath(dirname(shapefile_path),"catchment-y11-4326.geojson"), driver="GeoJSON")  
lpath = joinpath(dirname(shapefile_path),"catchment-y11-4326.geojson")
d = GeoDataFrames.read(lpath)
ed = reduce(GeoDataFrames.union,d.geometry)
cds = GeoInterface.coordinates(ed)
extrema.(cds)
flat_coords = vcat(cds...)
min_x = minimum(x -> x[1], flat_coords)
max_x = maximum(x -> x[1], flat_coords)
min_y = minimum(x -> x[2], flat_coords)
max_y = maximum(x -> x[2], flat_coords)

cd("D:/remo/cordex/wgs")
inf = "proj_pre_hist+rcp85_wgs84.nc"
#outf= "pre_rcp85-saale.nc"


import LibGEOS: centroid
cen = centroid(ed)
a,b=GeoInterface.coordinates(cen)
#10.04834792372862      #midpoint saale
#50.28952065588712
run(`wsl.exe -e cdo griddes $inf`)
xfirst    = 8.09
xinc      = 0.11
yfirst    = 50.93
yinc      = -0.1099999
grid_size = 15
outf = "/mnt/e/qq/prec_sq3.nc"
#min_x = min_x - 0.11
#y_start = y_start + (0.11 * 4 )

##DAS WAR SCHWER
x_start, x_end, y_start, y_end = cb_c(a,b,xinc,yinc,grid_size)
y_end = y_end - (0.11 * 1)
run(`wsl.exe -e cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $inf $outf`)
j = [x_start,x_end,y_start,y_end]
join(round.(j;digits=6),",")|>cb #9.278348,10.818348,51.05952,49.409521

#latx() #16mb
f = towin(outf)
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")
f = "E:/qq/rainfarm_0001.nc"
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")


/mnt/d/Relief_DGMs/FABDEM/franken_selection/
#gdal_translate -of NetCDF mosaic.vrt fabdem.nc
#1186.80MB  
inf="/mnt/d/Relief_DGMs/FABDEM/franken_selection/fabdem.nc"
outf="/mnt/e/qq/fabdem.nc"
run(`wsl.exe -e cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $inf $outf`)
run(`wsl.exe -e cdo griddes $outf`)
run(`wsl.exe -e cdo vardes $outf`)

#in wsl 
"/mnt/e/qq"
tools=$jlpt/tools 
prj=/mnt/d/remo/cordex/wgs/utm/cropped/rainfarm
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl
cdo griddes 
oro="/mnt/e/qq/fabdem.nc"
ref="/mnt/e/qq/prec_sq3.nc"
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfweights.jl $oro $ref
# cdof(towin(ref))
# lat()
#generated "weights.nc"
# errors julia --threads auto -q  --startup-file=no --project=$prj $tools/rfsmooth.jl --gaussian $ref psmooth.nc
#40 sec
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl --weights "weights.nc" $ref -o "rainf_weights.nc"
mv rainf_weights.nc_0001.nc rainf_weights.nc
f="rainf_weights.nc"
latx()
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")

julia --threads auto -q  --startup-file=no --project=$prj $tools/rfweights.jl -n 8 $oro $ref
llf
# n = "Subdivisions for downscaling"
#             arg_type = Int
#             default = 2
p="julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl"
$p --weights "weights.nc" -n 8 -o "rainf_weights2" $ref
#140.579336 second , 2639.50MB
f=lat()
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")


raw"/home/ubu/.julia/logs/repl_history.jl"
using ClimateTools
a=load("../rainf_weights2_0001.nc","pre")
import ClimatePlots as cps
using ClimatePlots
plot(annualmax(a))
contourf(a)

##cut pre-cor and add to rainfarm with smooth
"/mnt/e/qq"
#inf="/mnt/e/qq/pre/pre-cor_extremes.nc" #key weight err.
#outf="/mnt/e/qq/rainfarm/pre-cor2.nc"
inf="/mnt/e/qq/pre/pre-cor2.nc"
outf="/mnt/e/qq/rainfarm/pre-cor3.nc"
pwd()
run(`wsl.exe -e cdo griddes $inf`)

xfirst    = 9.04986035478427
xinc      = 0.0999999995999529
yfirst    = 49.0498597218401
yinc      = 0.0999999958273361
xsize     = 30
ysize     = 17

grid_size = 15
x_start, x_end, y_start, y_end = cb_c(a,b,xinc,yinc,grid_size)
#y_end = y_end - (yinc * 3)
#x_start, x_end, y_start, y_end = 9.278348,10.818348,49.409521,51.05952
#y_end = y_end + (yinc * 2)
#y_start = y_start + (yinc * 2)
y_start = 49.409521 - yinc
run(`wsl.exe -e cdo -v -P 4 sellonlatbox,$x_start,$x_end,$y_start,$y_end $inf $outf`)
run(`wsl.exe -e cdo griddes $outf`)
f = towin(outf)
f = towin(inf)
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")
j = [x_start,x_end,y_start,y_end]
join(round.(j;digits=6),",")|>cb
#9.348348,10.748348,49.309521,50.989521 #für corgrids

cd "/mnt/e/qq/rainfarm"
##Creating weights from file /mnt/e/qq/fabdem.nc
oro="/mnt/e/qq/fabdem.nc"
#ref="/mnt/e/qq/rainfarm/pre-cor2.nc"
ref="/mnt/e/qq/rainfarm/pre-cor3.nc"
jp="julia --threads auto -q --startup-file=no --project=$prj $tools/rfweights.jl"
$jp $oro $ref
#LoadError: KeyError: key "lat" not found
ncvarlst $ref
# cp -v ../weights.nc . #geht nicht
# rm weights.nc
jrf="julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl"
#$jrf $ref -o "rainf_cor3"
ref="/mnt/e/qq/rainfarm/pre-cor2.nc"
$jrf $ref -o "rainf_cor2"
#falsch weil transponiert.
f = "/mnt/e/qq/rainfarm/rainf_cor2_0001.nc"
f = towin(f)
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")

#make new biascorr using ClimateTools.
f="E:/qq/rainfarm/pre_biljl.nc"
pyjl.plot_mean_with_shapefile(f,lpath;myvar="pre")


####03.04.24


@pyjl
using .pyjl
f=raw"E:\qq\rainf_weights2_0001.nc"

import GeoDataFrames as gdf
ge = gdf.read("D:/Wasim/regio/rcm/ezg_4326.json")
rcmoutline = reduce(gdf.union,ge.geometry)
crds = GeoInterface.coordinates(rcmoutline)
flat_coords = vcat(crds...)
min_x = minimum(x -> x[1], flat_coords)
max_x = maximum(x -> x[1], flat_coords)
min_y = minimum(x -> x[2], flat_coords)
max_y = maximum(x -> x[2], flat_coords)

#j = [min_x, max_x, min_y, max_y]
#roundup to 4 digits
#j = round.(j, RoundUp; digits=4)
mins = [min_x, min_y]
maxs = [max_x, max_y]
mins = round.(mins, RoundDown; digits=4)
maxs = round.(maxs, RoundUp; digits=4)
cd(raw"E:\qq\rf2")
pwc()
/mnt/e/qq/
#cdo -v -f nc -sellonlatbox,9.278348,10.818348,49.409521,51.05952 -random,r1440x720 randomwgs.nc
#cdo -v -f nc -sellonlatbox,9.278348,10.818348,49.409521,51.05952 -random,r2880x1440 randomwgs2.nc
yxinc      = 0.125
#xmin - (2*yxinc)
xmin,xmax,ymin,ymax = mins[1],maxs[1],mins[2],maxs[2]
Xlowerleft = 9.05
xmin = Xlowerleft
xmax = xmin + 15 * yxinc
#ymin = mins[2] - (yxinc)
ymin = 49.0
ymax = ymin + 15 * yxinc
randomout="randomwgs4.nc"
#15x15
run(`wsl.exe -e cdo -v -f nc -sellonlatbox,$xmin,$xmax,$ymin,$ymax -random,r2880x1440 $randomout`)
#run(`wsl.exe -e cdo -v -f nc -sellonlatbox,$xmin,$xmax,$ymin,$ymax -random,r2880x1440 $randomout`)
#9.0724 / 0.125
#cdo -v -f nc -sellonlatbox,9.25,10.818348,49.41,51.06 -random,r2880x1440 randomwgs4.nc
#mkcd rf2
#13x13
#cdo -v -f nc -sellonlatbox,9.25,10.818348,49.41,51.06 -random,r2880x1440 randomwgs3.nc
#cdo griddes randomwgs2.nc > grid84.txt

cdo -v griddes randomwgs4.nc > grid84.txt
#xref="/mnt/e/qq/prc-sq13.nc"
#raw REMO
xref = "D:/remo/cordex/pre_hist+rcp85_19500102-21001231.nc"
xref = towsl(xref)
#xref="/mnt/e/qq/prc-sq13.nc"
#cdo -P 4 remapbil,grid84.txt $xref prc-sq13-wgs84.nc
ofl = "pre-15x15-wgs84.nc"

pwd()
#Bilinear weights from curvilinear (34x25) to lonlat (15x15) grid
run(`wsl.exe -e cdo -P 4 remapbil,grid84.txt $xref $ofl`)
#oro="/mnt/e/qq/fabdem.nc"
#ref="/mnt/e/qq/prc-sq13-wgs84.nc"
#ref="/mnt/e/qq/rf2/prc-sq13-wgs84.nc"
ref="/mnt/e/qq/rf2/pre-15x15-wgs84.nc"
tools=$jlpt/tools 
prj="/home/ubu/.julia/rainfarm" 
inf="/mnt/d/Relief_DGMs/FABDEM/franken_selection/fabdem.nc"
outf="/mnt/e/qq/rf2/fabdem2.nc"
#run(`wsl.exe -e cdo -v -P 4 -sellonlatbox,9.278348,10.818348,49.409521,51.05952 $inf $outf`)
#run(`wsl.exe -e cdo -v -P 4 -settaxis,2000-01-01 -sellonlatbox,9.25,10.818348,49.41,51.06 $inf $outf`)
# join([xmin,xmax,ymin,ymax],",")
# run(`wsl.exe -e cdo -v -P 4 -settaxis,2000-01-01 -sellonlatbox,8.8,10.675,48.8732,50.7482 $inf $outf`)

join([xmin,xmax,ymin,ymax],",")
run(`wsl.exe -e cdo -v -P 4 -settaxis,2000-01-01 -sellonlatbox,9.05,10.925,49.0,50.875 $inf $outf`)
#cdo(1) sellonlatbox: nlon=6750 nlat=6750
run(`wsl.exe -e cdo griddes $outf`)
#run(`wsl.exe -e cdo vardes $outf`)
oro="/mnt/e/qq/rf2/fabdem2.nc"
ref="pre-15x15-wgs84.nc"
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfweights.jl -n 4 $oro $ref
#fact4
lat()
run(`wsl.exe -e cdo griddes $ref`)
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl -h
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl -n 4 --weights "weights.nc" $ref -o "pre_remo_rainf"
#Computed spatial spectral slope: 1.777813158575936 #fact8
#Computed spatial spectral slope: 1.7010206183585082 #fact4
#154.398580 seconds (52.91 M allocations: 313.224 GiB, 3.88% gc time, 11.72% compilation time)

#Terzago haben 80 ensemble members

#ich mach mal 8
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl --nens 8 -n 4 --weights "weights.nc" $ref -o "pre_80_ens"

lv="pre_remo_rainf_0001.nc"
#
run(`wsl.exe -e cdo -v selyear,2013 $lv t.nc`)
sp = "D:/Wasim/regio/rcm/ezg_4326.json"
#pyjl.plot_mean_with_shapefile("t.nc",sp)
f="E:/qq/rf2/pre_80_ens_0006.nc"
pyjl.plot_mean_with_shapefile(f,sp)

#rm1 = mean(Raster(f,key=:pre),dims=Ti) #Variable 'time_bnds' not found

import Pkg; Pkg.add("Unitful")
using Unitful
time_in_seconds = 154.398580u"s"  # Declare a quantity in seconds
time_in_hours = uconvert(u"hr", time_in_seconds * 80)  # Convert to hours
time_in_hours = uconvert(u"hr", time_in_seconds * 8)
#3.43 hrs

s=246u"J/kg"
uconvert(u"m^2/s^2", s)
uconvert(u"km^2/hr^2", s)
a=9.81u"m/s^2"
b=0.125u"kg"
#F = m*a
F = b*a
uconvert(u"N", F) #Newton = kg m/s^2
uconvert(u"kN", F) #kiloNewton = 1000 N

julia --threads auto -q  --startup-file=no --project=$prj $tools/rfweights.jl -h
echo $ref
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfweights.jl -n 8 $oro $ref -o "weights_8.nc"
p="julia --threads auto -q  --startup-file=no --project=$prj $tools/rfarm.jl"
$p --weights "weights_8.nc" -n 8 -o "pre_rf_8" $ref
#313.224 GiB #Gibibyte binary 
#3029.97MB

#first regrid it to WGS84 via cdo remapbil
#then regrid it to 25832 via rio
#then rename the dimensions
#then save it
#then close it

#r = Raster(f,key=:pre)
xr = pyimport("xarray")
rio = pyimport("rioxarray")
#f="E:/qq/rf2/pre_v8_0001.nc"
#f="E:/qq/rf2/prc-sq13-wgs84.nc"
f=lat()
ds = xr.open_dataset(f)
#transform to 25832
ds = ds.rio.write_crs("EPSG:4326")
#ds = ds.transpose("t","x","y")
ds = ds.rio.reproject("EPSG:25832")
dk = ds.isel(time=100)
dk["pre"].plot()
#dm = ds["pre"].mean(dim="time")
#dm.plot()
#ds = ds.rename(Dict("lon"=>"x","lat"=>"y","time"=>"t"))
ds = ds.rename(Dict("time"=>"t"))
ds = ds.transpose("t","x","y")
#fo=replace(f,"rainf_weights2_0001" => "pre_rf")
fo=replace(f,"pre_rf_8_0001" => "pre_rf8_25832")
# differences in dimensions
println(ds.dims) #ist auchlt wasimdoku wurst
println(ds["pre"].dims) #das ist wichtig
ds.to_netcdf(fo)
ds=ds.close()
#now make lists and run it....
f
pyjl.plot_mean_with_shapefile(f,sp)

@doc pyjl.plot_mean_with_shapefile(fo,sp,"EPSG:25832")
pyjl.plot_mean_with_shapefile(fo,sp;mycrs="EPSG:25832") #dim is t
@doc pyjl.xrlist
@doc pyjl.xrlist2

run(`wsl.exe -e cdo griddes $fo`)

gres=1462.06578447128
@edit pyjl.xrlist2(;gridres=gres,xmatch="pre_rf8_25832",suffix=".winlist")
npplat()

fn=raw"E:\qq\rf2\pre_rf8_25832.winlist"
df = tsread(fn)
ncol(df)
#replace 12 with 24 in col 4
names(df)
#df[!,4] .= replace(df[!,4],12=>24)
replace!(df[!,4],"12"=>"24")
hd(df)
#write it back
#fn = replace(fn,".winlist24" => ".winlist_tst")
fn=raw"E:\qq\rf2\pre_rf8_25832.winlist24"
writedf(fn,df)
npplat()


fn="E:/qq/rf-cor/tas-qq2.nc"
pyjl.plot_mean_with_shapefile(fn,sp;myvar="tas")

cd("E:/qq/rf-cor")
cp("E:/qq/rf2/pre_rf8_25832.winlist24","rf-dly.list_24")
ls()
#filename::String, old_path::String, new_path::String)
pyjl.replace_path_in_file("rf-dly.list_24","E:\\qq\\rf2\\","stateini/")
##
npplat()

# old RE: $inpath_meteogrids//pr-dly.list_24
# new RF: $inpath_meteogrids//rf-dly.list_24

###test for tas clipping
inf="tas-qq2.nc"
# ext(540041.61, 625441.61, 5539909.291, 5606709.291) #y14
outf="tas-qq2-y14.nc"
run(`wsl.exe -e cdo griddes $inf`)
run(`wsl.exe -e cdo -v -P 4 -sellonlatbox,540041.61,625441.61,5539909.291,5606709.291 $inf $outf`)
xinc=8317.01
540041 + ( 26 * xinc )
# idx1 INTEGER Index of ﬁrst longitude (1 - nlon)
# idx2 INTEGER Index of last longitude (1 - nlon)
# idy1 INTEGER Index of ﬁrst latitude (1 - nlat)
# idy2 INTEGER Index of last latitude (1 - nlat)
#gives me 23x23 grid
run(`wsl.exe -e cdo -v -P 4 -selindexbox,1,23,1,23 $inf $outf`)
run(`wsl.exe -e cdo griddes $outf`)
# cdo sinfo $lv
#inwsl:
oro="/mnt/e/qq/rf2/fabdem2.nc"
ref="tas-qq2-y14.nc"
julia --threads auto -q  --startup-file=no --project=$prj $tools/rfweights.jl -n 4 $oro $ref
#cdo sinfo tas_cor_raw.nc|head -10

#1st step#StaionData+GenRE->PEST->eval
#2nd step#remo->eobscor->rainfarm->project->wasim

fn = @gl "unix"
cp(fn,"rf-dly.eobs-unixlist")
pyjl.replace_path_in_file("rf-dly.eobs-unixlist","stateini/","stateini/eobs/")

#weiter mit Temperature 
jlrfr rftemp.jl -h 
oro="/mnt/e/qq/rf2/fabdem2.nc"
ref="/mnt/e/qq/tas/tas_cor2.nc"

griddes $ref
# -l, --lapse LAPSE     Lapse rate [K/km]
# -r, --radius RADIUS   Smoothing radius (in grid units)
jlrfr rftemp.jl -l 1.2 -r 0.5 $oro $ref "tas_rf8_25832"



