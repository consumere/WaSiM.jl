## Modules 
using NCDatasets     # open and manipulate NetCDFs
using Plots          # generate simple plots
using Dates          # to work with dates and time indices

# Open the file and show the metadata 
nc = NCDataset("./path_data/GLO_BIO_Dec2022.nc")
println(nc) 
 
# Read dimensions
lon = nc["longitude"][:]
lat = nc["latitude"][:]
time = nc["time"][:]
depth = nc["depth"][:]

# Load all data of a variable (here pH)
ph = nc["ph"][:]

# Load the attributes
ph.attrib

# Define day and depth level of interest
date = DateTime(2022, 12, 15, 12)
depth_level = 50

# Extract the indices 
time_indices = findall(date .== time)[1]   
depth_indices = findall(depth_level .<= depth)[1]

# Generate the plot 
heatmap(lon, lat, ph[:,:, depth_indices, time_indices]', 
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string("",ph.attrib["long_name"],, " at ", depth_level," meters, on ", date))


########################
#crop many nc files
"D:/wslconda/cmip6/"|>cd
xl = nconly("lnd")

@vv "poly"
using GeoDataFrames
shp="D:/Relief_DGMs/FABDEM/franken_fabdem_domain.shp"
gdf = GeoDataFrames.read(shp)
p = gdf[!,1][1]
@vv "reproject"
using ArchGDAL
const AG = ArchGDAL
p = AG.reproject(p,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
tmp=Raster(first(xl))

mask_trim2(raster, poly) = Rasters.trim(Rasters.mask(raster; with=poly); pad=0)
ot = mask_trim2(tmp,p)
@vv "Where"
@vv "At("
ot[Z(At(0.05))]|>plot
AG.write(ot,"D:/wslconda/cmip6/lnd/trim.nc")
using Makie
using CairoMakie
ot[Z(At(0.05))]|>Rasters.rplot
ot[Ti(At(DateTime("2015-01-16T12:00:00")))]|>Rasters.rplot


#Base.write(filename::AbstractString, A::AbstractRaster; kw...)
Base.write("tst.nc", ot;)
op()

#now in loop
using GeoDataFrames
shp="D:/Relief_DGMs/FABDEM/franken_fabdem_domain.shp"
gdf = GeoDataFrames.read(shp)
p = gdf[!,1][1]
using Rasters, ArchGDAL
const AG = ArchGDAL
p = AG.reproject(p,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
mask_trim2(raster, poly) = Rasters.trim(Rasters.mask(raster; with=poly); pad=0)
# ot = mask_trim2(tmp,p)
"D:/wslconda/cmip6/"|>cd
xl = nconly("lnd")
#ot.name
xl[5]
import NCDatasets
for i in xl
        tmp=Raster(i)
        ot = mask_trim2(tmp,p)
        Base.write(string(i), ot;force=true)
        println(string(i)," over-written")
end

@edit pfix()





lk="D:/Wasim/regio/rcm200/v4/rcm.lin"
lk="D:/Wasim/regio/rcm200/v6/rcm.lin"
lk="D:/Wasim/regio/rcm200/v9/rcm.lin"
lk="D:/Wasim/regio/rcm200/v2/rcm.lin"

lk="D:/Wasim/regio/rcm200/rcm.lin"
@edit rplot(lk)
rplot(lk)

lk="D:/Wasim/regio/rcm200/v11/rcm.lin"
rx=Raster(lk)
Rasters.rplot(rx)

plotly()
gr()
zp(rplot)|>println

xr = read(Raster(lk;crs=EPSG(25832),missingval=-9999))
plot(xr)
wa.rplot(lk;naval=0)
wa.rplot(lk;naval=-9999)

