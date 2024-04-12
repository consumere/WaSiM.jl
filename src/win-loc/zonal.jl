#zonal stats

#from https://github.com/rafaqz/Rasters.jl/issues/312
using Rasters
import ArchGDAL
using Shapefile
using Statistics
using DataFrames
using PyCall

#shapefilepath = "D:/Wasim/Tanalys/DEM/brend_fab/in3/catchment-in3.shp"
shapefilepath = "/mnt/d/Wasim/Tanalys/DEM/brend_fab/in3/catchment-in3.shp"
#rasterpath = "D:/Wasim/Tanalys/DEM/brend_fab/in3/fab.slp"
rasterpath = "/mnt/d/Wasim/Tanalys/DEM/brend_fab/in3/fab.dhk"

# ri = pyimport("rasterio")
# d = ri.open(rasterpath)


py"""
import rasterio
import numpy as np
from rasterio.mask import mask
import geopandas as gpd
def pyzonal(shapefilepath, rasterpath, agg=np.ma.sum):
    shp = gpd.read_file(shapefilepath)  
    raster = rasterio.open(rasterpath)
    def geom_mask(geom, dataset=raster, **mask_kw):
        mask_kw.setdefault('crop', True)
        mask_kw.setdefault('all_touched', False)
        mask_kw.setdefault('filled', False)
        masked, mask_transform = mask(dataset=dataset, shapes=(geom,),
                                    **mask_kw)
        return masked
    return shp.geometry.apply(geom_mask).apply(agg)
"""

function jlzonal(shapefilepath, rasterpath;agg=mean)
    shp = Shapefile.Table(shapefilepath)
    raster = Raster(rasterpath)
    Float64.(zonal(agg, raster; of=shp, boundary=:center))
    #Int.(zonal(sum, raster; of=shp, boundary=:center))
end
#boundary=:touches))

@time python_zonal = convert(Array{Int64}, py"""pyzonal($shapefilepath, $rasterpath)""")
#  3.100880 seconds (156 allocations: 4.016 KiB)

@time julia_zonal = jlzonal(shapefilepath, rasterpath)
#  0.820940 seconds (14.85 k allocations: 175.825 MiB, 0.88% gc time)

difference = python_zonal .- julia_zonal;

DataFrame(; julia_zonal, python_zonal, difference)




######
using ImageCore
using Statistics
using Plots
using Rasters
using DimensionalData
const DD = DimensionalData

function eachband(r::Raster)
    bands = dims(r, Band)
    return (view(r, Band(b)) for b in bands)
end

function normalize!(raster, low=0.1, high=0.9)
    for band in eachband(raster)
        l = quantile(skipmissing(band), low)
        h = quantile(skipmissing(band), high)
        band .-= l
        band ./= h - l + eps(float(eltype(raster)))
        band .= clamp.(band, zero(eltype(raster)), one(eltype(raster)))
    end
    return raster
end

function plot_raster(r::Raster; bands=[1,2,3], low=0.02, high=0.98)
    img = float32.(copy(r[Band([bands...])]))
    normalize!(img, low, high)
    img = permutedims(img, (Band, X, Y))
    img = DimensionalData.reorder(img, DD.ForwardOrdered)
    x = DD.index(reorder(dims(img, X), DD.ForwardOrdered))
    y = DD.index(reorder(dims(img, Y), DD.ForwardOrdered))
    plottable_img = colorview(RGB, parent(img))
    Plots.plot(x,y,plottable_img,
               title = string(name(r)),
               xlabel = Rasters.label(dims(r, X)),
               ylabel = Rasters.label(dims(r, Y)),
               )
end

plot_raster(r)





#########


import Shapefile
shp="D:/Wasim/regio/rcm200/v12/catchment.shp"
rpt="D:/Wasim/regio/rcm200/v12/rcm.use"

cd("C:/Users/chs72fw/Documents/EFRE_GIS/Landcover/20181016_Deutschland_LC_clip_for_Ullmann/")
rglob("tif")
rpt=".\\v2\\LULC_DE_2014_nbg_200m.tif"
#r = agread(rpt)
r = Raster(rpt)

rn = rst.project(r;dst="EPSG:25832",
    src="EPSG:3034",imethod="near")

shp = Shapefile.Table(shp)
raster = rn
Float64.(zonal(median, raster; of=shp, boundary=:center))
zq = Int64.(zonal(minimum, raster; of=shp, boundary=:center))
zmx = Int64.(zonal(maximum, raster; of=shp, boundary=:center))


#c = Shapefile.Handle(shp)
#cc = wa.reproject(c.shapes;crs=3034)
using GeoDataFrames
c = GeoDataFrames.read(shp)


@edit wa.reproject(c.geometry;crs=EPSG(3034))
#@edit jlzonal(shp,rpt)
zz = rst.jlzonal(shp,rpt;agg=sum)

