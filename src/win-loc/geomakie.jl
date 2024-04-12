using GeoMakie, CairoMakie
import GeoMakie as gm
#using GeoMakie.GeoJSON

#test
# fig = Figure()
# ga = GeoAxis(fig[1, 1]; dest = "+proj=wintri")  # Specify the CRS for plotting
# gm.lines!(ga, GeoMakie.coastlines())  # Plot coastlines from Natural Earth as a reference
# gm.scatter!(ga, -120:15:120, -60:7.5:60; color = -60:7.5:60, strokecolor = (:black, 0.2))
# fig
#GeoDataFrames.read(raw"D:\Wasim\main_basin.geojson")
#https://geo.makie.org/stable/
# First, make a surface plot
lons = 9:12
lats = 48:51
field = [exp(cosd(l)) + 3(y/90) for l in lons, y in lats]

fig = Figure()
ax = GeoAxis(fig[1,1])
sf = gm.surface!(ax, lons, lats, field; shading = NoShading)
cb1 = gm.Colorbar(fig[1,2], sf; label = "field", height = Relative(0.65))
#countries_file = GeoMakie.assetpath("vector", "countries.geo.json")
#fr=raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\saale\elevation_data\merit_adapted_basins.geojson"
#countries_file = raw"D:\Wasim\main_basin.geojson"
#countries_file = fr
#countries_file = "D:/Relief_DGMs/FABDEM/wasim/lowres/saale.geojson"
countries_file = "D:/Wasim/regio/rcm/ezg_4326.json"
countries = GeoMakie.GeoJSON.read(read(countries_file, String))
n = length(countries)
hm = gm.poly!(ax, countries; color= 1:n, colormap = :dense,
    strokecolor = :black, strokewidth = 0.5,
)
gm.translate!(hm, 0, 0, 100) # move above surface plot

fig

fn = raw"D:\remo\qm\corgrids\prec\pre"
fn = "D:/remo/proj_main/APRR_2000-2017_25832_daysum.nc"
fn = "D:/remo/remo_for_wasim/2000/APRR_2000-2000_ufr_daymean.nc"
#@pyimport xarray as xr
#ds = xr.open_dataset(fn)
#ds.close()
r = Raster(fn)
dims(r)
#rmo = mean(r,dims=Dim{:t})
rmo = sum(r,dims=Ti)

#fig = cmk.mkrheat(rmo);
fig = Rasters.rplot(rmo;);
gm.poly!(countries)
fig


lons = lookup(rmo,:X)
lats = lookup(rmo,:Y)
field = rmo[:,:,1]
fig = Figure()
ax = GeoAxis(fig[1,1])
sf = gm.surface!(ax, parent(lons), parent(lats), field; shading = NoShading)

field = rmo[:,:,1]
fig = Figure()
#ax = Axis(fig[1,1])
ax = GeoAxis(fig[1, 1]; dest = "+proj=merc", limits=((9, 12), (48, 52)))
sf = gm.heatmap!(field)
fig
cb1 = gm.Colorbar(fig[1,2], sf; label = "Precipitation",
 height = Relative(0.65))
fn = "D:/Wasim/regio/rcm/ezg_4326.json"
pol = GeoMakie.GeoJSON.read(read(fn, String))
n = length(pol)
hm = gm.poly!(ax, pol; color= 1:n, colormap = :dense,
    strokecolor = :black, strokewidth = 0.5,
)
gm.translate!(hm, 0, 0, 100) # move above surface plot
fig


fn = "D:/Wasim/regio/rcm/ezg_4326.json"
pol = GeoMakie.GeoJSON.read(read(fn, String))
fig = Figure()
#ga = GeoAxis(fig[1, 1]; dest = "+proj=merc", limits=((9, 12), (48, 52)))
#ga = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84", limits=((9, 12), (48, 52)))
poly!(ga, pol; strokewidth = 0.7, color=:gold, rasterize = 5)
fig


r = Raster(fn)
dims(r)
#rmo = mean(r,dims=Dim{:t})
rmo = sum(r,dims=Ti)
fn = "D:/Wasim/regio/rcm/ezg_4326.json"
pol = GeoMakie.GeoJSON.read(read(fn, String))

lons = lookup(rmo,:X)|>parent
lats = lookup(rmo,:Y)|>parent
field = rmo[:,:,1]|>collect
fig = Figure()
#ax = GeoAxis(fig[1,1]; dest = "+proj=merc", limits=((9, 12), (48, 52)))
ax = GeoAxis(fig[1,1]; limits=((8, 13), (48, 52)))
sf = gm.surface!(ax, lons, lats, field; shading = NoShading)
cb1 = gm.Colorbar(fig[1,2], sf; label = "Precipitation",
 height = Relative(0.65))
n = length(pol)
hm = gm.poly!(ax, pol; #color= 1:n, #colormap = :dense, #Reverse(:plasma)
color = :transparent,    
strokecolor = :black, strokewidth = 1.5,
)
gm.translate!(hm, 0, 0, 100) # move above surface plot
fig


#make a func out of it:

"""
rmo = sum(x,dims=Ti)
fn = "D:/Wasim/regio/rcm/ezg_4326.json"
"""
function gmplot(rmo::Raster, fn::Union{String,GeoMakie.GeoJSON.FeatureCollection};lims::Tuple{Tuple{Int,Int},Tuple{Int,Int}}=((8, 13), (48, 52)),c=Reverse(:plasma))
    if fn isa GeoMakie.GeoJSON.FeatureCollection
        pol = fn
    else
        pol = GeoMakie.GeoJSON.read(read(fn, String))
    end
    #rdims = dims(rmo)
    lons = lookup(rmo,:X)|>collect #parent
    lats = lookup(rmo,:Y)|>collect #parent
    
    field = rmo[:,:,1]|>collect
    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; 
        limits=lims)
    sf = GeoMakie.surface!(ax, lons, lats, field; 
        #colormap = Reverse( :plasma, 0.5 ),     
        #colormap = ["grey", "red", "orange", "yellow", "blue"],
        colormap = c,
        shading = NoShading)
    cb1 = GeoMakie.Colorbar(fig[1,2], sf; 
        label = "Precipitation",
        height = Relative(0.65))
    n = length(pol)
    hm = GeoMakie.poly!(ax, pol; color = :transparent,    
    strokecolor = :black, strokewidth = 1.5,
    )
    GeoMakie.translate!(hm, 0, 0, 100) # move above surface plot
    return fig
end

fn="D:/remo/qm/corgrids/pre/pre-cor.nc"
pr = Raster(fn,key=:pre)
pr.dims|>last|>size|>first|>typeof
tdim=pr.dims|>last
yrs = pr.dims|>last|>size|>first 
yrs = yrs ÷ 365
cy = sum(pr,dims=tdim) ./ yrs
dims(cy)
gmplot(cy, "D:/Wasim/regio/rcm/ezg_4326.json")

cmt = gdf.read("D:/Wasim/regio/rcm200/v13/cmtv13.shp")
#@doc project
#cmt = project(cmt;dst="+proj=longlat +datum=WGS84")
@rimport terra as tr
ve = tr.vect("D:/Wasim/regio/rcm200/v13/cmtv13.shp")
ve = tr.project(ve, "+proj=longlat +datum=WGS84 +no_defs")
tr.writeVector(ve,"D:/Wasim/regio/rcm200/v13/cmtv13_4326.json")

gmplot(cy, "D:/Wasim/regio/rcm200/v13/cmtv13_4326.json";
    lims=((8, 11), (49, 51)))

mycm=["grey","orange", "yellow", "blue"]

k=GeoMakie.GeoJSON.read(read(fn, String))
gmplot(cy, k;
    lims=((8, 12), (49, 51)),c=mycm)


########################wasimrasters#########
@doc readras
r="D:/Wasim/regio/out/rc200/n3/thetrcm__stack.2010.nc"
#r = Raster(r;crs=Rasters.EPSG(25832),mappedcrs=Rasters.EPSG(25832),missingval=-9999)
r = readras(r)
ra = r[t=3]
ra = replace_missing(ra,0)
ra = rst.project(ra;dst="+proj=longlat +datum=WGS84")
cmk.mkrheat(ra)
gmplot(ra, k;
    lims=((8, 12), (49, 51)),c=mycm)
#GeoMakie.heatmap(ra)

    #lons = lookup(ra,:X)|>collect
    #lats = lookup(ra,:Y)|>collect
    #field = ra[:,:]|>collect
    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1];)
    sf = GeoMakie.heatmap!(ax,ra;       
    colormap = ["grey", "red", "orange", "yellow", "blue"],
    shading = NoShading)
    cb1 = GeoMakie.Colorbar(fig[1,2], sf; 
        label = "SB",
        height = Relative(0.65))
    n = length(pol)
    hm = GeoMakie.poly!(ax, pol; color = :transparent,    
    strokecolor = :black, strokewidth = 1.5,
    )
    GeoMakie.translate!(hm, 0, 0, 100) # move above surface plot
    fig

"""
wasim 
"""
function gmwas(rmo::Union{String,Raster}, fn::Union{String,GeoMakie.GeoJSON.FeatureCollection};
        lyr::Int = 1,
        lims::Tuple{Tuple{Int,Int},Tuple{Int,Int}}=((8, 13), (48, 52)),c=Reverse(:plasma))
    if fn isa GeoMakie.GeoJSON.FeatureCollection
        pol = fn
    else
        pol = GeoMakie.GeoJSON.read(read(fn, String))
    end
    if rmo isa Raster
        ra = rmo[t=lyr]
        ra = replace_missing(ra,0)
        ra = rst.project(ra;dst="+proj=longlat +datum=WGS84")
    else
        r = Raster(rmo)
        ra = r[t=lyr]
        ra = replace_missing(ra,0)
        ra = rst.project(ra;dst="+proj=longlat +datum=WGS84")
    end
    
    fig = GeoMakie.Figure()
    ax = GeoMakie.GeoAxis(fig[1,1]; 
        limits=lims)
    sf = GeoMakie.heatmap!(ax,ra;
        #colormap = Reverse( :plasma, 0.5 ),     
        #colormap = ["grey", "red", "orange", "yellow", "blue"],
        colormap = c,
        shading = NoShading)
    cb1 = GeoMakie.Colorbar(fig[1,2], sf; 
        label = "",
        height = Relative(0.65))
    n = length(pol)
    hm = GeoMakie.poly!(ax, pol; color = :transparent,    
    strokecolor = :black, strokewidth = 1.5,
    )
    GeoMakie.translate!(hm, 0, 0, 100) # move above surface plot
    return fig
end

s="D:/Wasim/regio/out/rc200/n3/thetrcm__stack.2010.nc"
k="D:/Relief_DGMs/FABDEM/wasim/lowres/saale.geojson"
gmwas(s, k;lyr=2,
    lims=((8, 12), (49, 51)),c=mycm)


fn="D:/remo/qm/corgrids/pre/pre-cor_extremes.nc"
pr = Raster(fn,key=:pre)
dims(pr)
yrs = pr.dims|>last|>size|>first
yrs = yrs ÷ 365
cy = sum(pr,dims=tdim) ./ yrs
dims(cy)
gmplot(cy, "D:/Wasim/regio/rcm/ezg_4326.json",
    lims=((8, 12), (49, 51)),c=["white","orange","blue","black"])

fn = "D:/remo/genRE/genRE_precipitation_hour_1.1.nc"
using PyCall
@pyimport xarray as xr
ds = xr.open_dataset(fn)
ds.keys()
da = ds.where(ds["time.year"]==2016,drop=true)
dx = ds.where(
    (ds["time.year"] >= 2015) &
    (ds["time.year"] <= 2016), drop=true)

dx = dx.drop_vars(["lat", "lon"])
ds = dx.resample(time="D").sum()
dsn = ds.rename(Dict("pr"=>"pre"))
dsn.attrs["units"] = "mm"
#dsn.to_netcdf("D:/remo/genRE/pre_15-16.nc")
dsn.to_zarr("D:/remo/genRE/pre_15-16.zarr") #50mb
da.close()
ds.close()

dx = xr.open_dataset(fn).drop_vars(["lat", "lon","crs"])
dx = dx.where(
    (ds["time.year"] >= 2000) &
    (ds["time.year"] <= 2016), drop=true)
dx = dx.resample(time="D").sum()
dsn = dx.rename(Dict("pr"=>"pre"))
dsn.attrs["units"] = "mm"
dsn.to_zarr("D:/remo/genRE/pre_00-16.zarr") #50mb
#dsn.to_netcdf("D:/remo/genRE/pre_dn.nc")
df = dsn.resample(time="Y").sum()
#dsn.resample(time="Y").sum().to_netcdf("D:/remo/genRE/X.nc")



pr = Raster(fn,key=:pr;lazy=true,crs=Rasters.EPSG(25832),mappedcrs=Rasters.EPSG(25832))
dims(pr)
#yrs = pr.dims|>last|>size|>first
prmean = mean(pr,dims=Ti) #takes too long...



s="D:/Wasim/regio/out/lowres/adj/v0/snowrcm_1600.sum.nc"
s="D:/Wasim/regio/out/lowres/adj/v0/rad_rcm_1500.mit.nc"
k="D:/Relief_DGMs/FABDEM/wasim/lowres/saale.geojson"
gmwas(s, k;lyr=1,lims=((8, 12), (49, 51)))


import Downloads, GeoJSON
geoger = GeoJSON.read(read(Downloads.download("https://raw.githubusercontent.com/isellsoap/deutschlandGeoJSON/main/2_bundeslaender/4_niedrig.geo.json"), String))
first(geoger.features).properties
by = first(filter(v -> occursin("Bayern", v.properties[:name]), geoger.features))

lakes = GeoJSON.read(read(Downloads.download("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_lakes_europe.geojson"), String))
#first(lakes.features).properties|>println
lakes.geometry|>first

fig = gm.Figure()
ga = gm.GeoAxis(fig[1, 1]; dest = "+proj=merc", limits=((9, 15), (47, 51)),
title="Simulation topology") 
#xgridvisible=false, ygridvisible=false, 
# xticklabelsvisible=false, yticklabelsvisible=false, 
# xticksvisible=false, yticksvisible=false)
gm.poly!(ga, geoger; strokewidth = 1, color=:peru, rasterize=5)
gm.poly!(ga, lakes; strokewidth = 1, color=:blue, rasterize = 5,  xautolimits=false, yautolimits=false)
fig



rglob(r"json","D:/Relief/DLM250/dlm250_ufrasubset")
rglob("shp","D:/Relief/DLM250/dlm250_ufrasubset/")
rglob("nc","kadf")
#make json from dlm
@rimport terra as tr
ve = tr.vect("D:/Relief/DLM250/dlm250_ufrasubset/DLM250_Fliessgew_Ufra.shp")
ve = tr.project(ve, "+proj=longlat +datum=WGS84 +no_defs")
tr.writeVector(ve,"D:/Relief/DLM250/dlm250_ufrasubset/DLM250_Fliessgew_Ufra_4326.json")
using GeoMakie, Makie
gew = GeoMakie.GeoJSON.read("D:/Relief/DLM250/dlm250_ufrasubset/DLM250_Fliessgew_Ufra_4326.json")
# Convert the MultiLineString objects to arrays of coordinates
lines = [GeoInterface.coordinates(line) for line in gew];
#The let variant creates an extra hard scope, but it doesn’t change anything if it’s within a function scope.
#It is commonly used for temporary bindings or to create a local context.
mymap = let   
    # Create a new figure
    fig = Figure()
    # Create a new GeoAxis
    ga = GeoAxis(fig[1, 1]; dest = "+proj=merc", limits=((9, 12), (49, 51)),
    title="Simulation topology") 
    poly!(ga, geoger; strokewidth = 1, color=:peru, rasterize=5)
    # Plot the lines
    for line in lines
        # Each line is an array of arrays of coordinates, so we need to flatten it
        flattened_line = reduce(vcat, line)
        # Separate the x and y coordinates
        x_coords = first.(flattened_line)
        y_coords = last.(flattened_line)
        # Plot the line
        lines!(ga, x_coords, y_coords;color=:blue)
    end
    fig    
end
mymap


#pegel
@vv "lpro"
peg = wa.lpro("D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt")
pts = [GeoInterface.coordinates(p) for p in peg.geometry];
#x,y = first.(GeoInterface.coordinates.(peg.geometry)), last.(GeoInterface.coordinates.(peg.geometry));
ug = GeoJSON.read("D:/Relief_DGMs/FABDEM/wasim/lowres/saale.geojson")

#https://geo.makie.org/stable/examples/most_projections/
fig = Figure()
ga = GeoAxis(fig[1, 1]; dest = "+proj=comill", limits=((9, 12), (49, 51)),
#ga = GeoAxis(fig[1, 1]; dest = "+proj=merc", limits=((9, 12), (49, 51)),
title="Simulation topology") 
# for line in lines
#     # Each line is an array of arrays of coordinates, so we need to flatten it
#     flattened_line = reduce(vcat, line)
#     # Separate the x and y coordinates
#     x_coords = first.(flattened_line)
#     y_coords = last.(flattened_line)
#     # Plot the line
#     lines!(ga, x_coords, y_coords;color=:blue, alpha=0.5)
# end
gm.poly!(ga, ug; strokewidth = 1.5, color=:transparent)
gm.scatter!(ga, first.(pts), last.(pts)) #elegant...
text!(first.(pts), last.(pts), text = string.(peg.name), align = (:right, :top), color = :black)
fig


k="D:/Wasim/regio/rcm/pestpp/res3/main-v2/out/gwstrcm.c6.2016"
cmk.tsp(readdf(k))

df = readdf(k)
cmk.cloudplot2(df)

z=raw"D:\Wasim\regio\rcm\pestpp\res3\main-v2\out\v4\snowrcm_1600.sum.nc"
z=raw"D:\Wasim\regio\rcm\pestpp\res3\main-v2\out\v4\sb05rcm_1600.mit.nc"
z=raw"D:\Wasim\regio\rcm\pestpp\res3\main-v2\out\v4\gwbalance1rcm.sum.nc"
mycm = ["grey","orange", "yellow", "blue"]
mycm = ["white","grey", "blue"]
gmwas(z, "D:/Wasim/regio/rcm/ezg_4326.json";lyr=1,
    lims=((8, 12), (49, 51)),c=mycm)

cdof(z)
map(dfr,glob(r"qoutjl$")) |>  λ -> map(skipyr, λ) |>  λ ->map(byear, λ) |>  λ -> plot_grouped_metrics( λ;col=:ve)
jdf = GeoJSON.read("D:/Wasim/regio/rcm200/v13/cmtv13_4326.json")
z="D:/Wasim/regio/out/rc200/x34/rad_rcm_1500.mit.nc"
z="D:/Wasim/regio/out/rc200/x34/sb05rcm_1400.mit.nc"
gmwas(z, jdf;lyr=1,lims=((9, 11), (50, 51)),c=mycm)
rplot(z)

@vv "pre_rcm85.wa2"
k="D:/remo/qm/corgrids/pre_rcm85.wa2"
#stps,df = wa.dfr(k;station=true,pts_and_data=true)
stps,df = wa.dfr(k;pts_and_data=true)
addplot(stps)

#fl = CSV.read(k,DataFrame;limit=4)
#xc = fl[2,5:end]|>collect

r = readras(z)
Plots.plot(r;
c=cgrad(:matter),xlabel="",
ylabel="",title = "",xaxis=false,yaxis=false,legend=false,grid=false)
addplot(stps)

wk=selt(df,r"kuppe")
ki=selt(df,r"Kissi")
#@edit hyeval(mall(wk,ki);freq="Y")
wa.hyeval(mall(wk,ki);freq="Y",grid=false,xaxis=false)
wa.hyeval(mall(wk,ki);freq="Q",grid=false,xaxis=false)

ri = nconly("ifl")|>last|>agread
nconly("ifl")|>last|>cmk.mkrheat
nconly("ifl")|>last|>agheat

fn="D:/remo/qm/corgrids/jlcor/tas-cor2.nc"
@pj
@edit 
pyjl.fdoy(fn;txy=true)


fn="D:/Wasim/regio/rcm200/v13/cmtv13_4326.json"
import GeoDataFrames as gdf
p = gdf.read(fn)
r=Raster(raw"E:\jl_LC\EarthEnv\LandCover\without_DISCover\Consensus_reduced_class_9.tif")
r = crop(r,to=p.geometry)
plot(r)
cd(raw"E:\jl_LC\EarthEnv\LandCover\without_DISCover")
mkdir("lc_crop")
cd("lc_crop")
lst = rglob("tif","../")
for i in lst
    r = Raster(i)
    r = crop(r,to=p.geometry)
    otname=replace(basename(i),"Consensus_reduced"=>"lc")
    write(otname,r)
    println("$otname written")
end
ser = RasterStack(glob(r"lc"))
write("all.nc",ser)
r = Raster("all.nc")
plot(r)
plot!(p.geometry,fillcolor=false)
map(rm,lst)
cd("..")
du()
csize()