@time setup()

#@time_imports 
using GeoDataFrames
using GeoInterface
const gdf = GeoDataFrames
n="D:/Wasim/regio/rcm200/v12/catchment.shp"
g = gdf.read(n)
plot(g.geometry, fillcolor=false)
pol = ArchGDAL.IGeometry[]
for i in 1:size(g,1)
    polygon = g.geometry[i]  # assuming the first geometry is a polygon
    src = EPSG(25832)
    dst = EPSG(4326)
    reprojected_polygon = ArchGDAL.reproject(
        polygon, src, dst) #;order=:compliant
    push!(pol,reprojected_polygon)
end

ArchGDAL.bounds(pol[1])
ArchGDAL.gety(pol[1],0)

function reverse_coords(polygon)
    # Extract the coordinates from the polygon
    crds = GeoInterface.coordinates.(pol)
    # Reverse each pair of coordinates
    reversed_crds = [[[(lon, lat) for (lat, lon) in ring] for ring in polygon] for polygon in crds]
    # Create a new polygon with the reversed coordinates
    reversed_polygon = ArchGDAL.createpolygon.(reversed_crds)
    return reversed_polygon
end

reversed_polygon = reverse_coords(pol)

crds = GeoInterface.coordinates.(pol)
# Reverse each pair of coordinates
reversed_crds = [[[(lon, lat) for (lat, lon) in ring] for ring in polygon] for polygon in crds]
pls = AG.createpolygon.(reversed_crds)
plot(pls)

g.geometry .= pls
plot(g.geometry, fillcolor=false)


AG.createpolygon.(reversed_crds[end-4:end])|>plot
map(size,reversed_crds)
map(typeof,pol)
ArchGDAL.getgeom(pol[5], 1)
H = map(x->ArchGDAL.getgeom(x, 1),pol)
#ArchGDAL.empty!.(H)
typeof.(H)
filter(x->typeof(x)==ArchGDAL.IGeometry{ArchGDAL.wkbPolygon},H)

ArchGDAL.IGeometry{ArchGDAL.wkbUnknown}

@doc ArchGDAL.flattento2d!
ArchGDAL.flattento2d!(H[1])
reverse_coords()


os = [AG.createpolygon.(fa[1]) for fa in reversed_crds]

fn="D:/Relief/Modis_Data/Franken_2016-2020/MCD15A2H.006_500m_aid0001.nc"
laif = Raster(fn; name = "Lai_500m") #, crs=CRS("+init=epsg:4326"))
# descr(laif)
# plot(laif[Ti=5])
# using Shapefile
# poly = Shapefile.Handle("D:/Wasim/regio/rcm200/v12/catchment.shp")
# #poly = reproject(poly, EPSG(4326))
# @vv "polygon"
# #lai = wa.mask_trim(laif, poly.shapes[1:end])
# lai = wa.mask_trim(laif, pol)
# replace_missing!(lai, 0)
# plot(lai[Ti=5])
# plot(lai[X(10),Y(10)])

da = wa.ncdf(laif)
dfp(da)
dfm(da;fun=yrmean,mode=:bar)

x="D:/Wasim/Tanalys/DEM/Input_V2/meteo/met0/wind_1970.txt"
p = wa.stp(x;repro=true)
p.geometry
dou = collect(extract(laif, p.geometry;atol=.1)|>skipmissing)
dx = DataFrame(dou)
GeoDataFrames.getgeometrycolumns(dx)
# k = dx.Lai_500m[3]
# plot(k, title=dx.geometry[1])
# typeof(k)
# #make station df
# df = DataFrame([Rasters.lookup(k,1),
# 	Float64.(k.data)],[:date,:Lai_500m])
# dfp(df)
#make station df loop ################
df = []

for i in 1:nrow(dx)
	k = dx.Lai_500m[i]
    replace_missing!(k,0)
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,:Lai_500m]))
end
df

ds = innerjoin(df..., on= :date, makeunique=true)
#ds = map(x->replace.(x, 25.0 => 0),ds)
#ds = transform(ds, propertynames(ds)[2:ncol(ds)] .=> (x -> replace.(x, 25.0 => 0.0)) .=> propertynames(ds)[2:ncol(ds)])
#ds = transform(ds, names(ds)[2:ncol(ds)] .=> 
#    (x -> eltype(x) <: Number ? replace.(x, 25.0 => 0.0) : x) .=> names(ds)[2:ncol(ds)])
# for c in propertynames(ds)[2:ncol(ds)]
#     ds[!, c] .= replace.(ds[!, c], 25.0 => 0.0)
# end

collect(eachcol(df))
replace!.(eachcol(ds[!,Not(1)]), 25.0 => 0.0)  ##!!!!
dfp(ds)
findmax.(eachcol(ds[!,Not(1)])) 
mean.(eachcol(ds[!,Not(1)]))
nn = map(y->(y[1:5]*y[end-1]) ,string.(p.name))
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
dfp(ds)
wa.wawrite(ds,"lai_ext.wa")
op()

da = qrtr(ds;fun=mean)
using StatsPlots
# Melt the DataFrame to a long format suitable for boxplots
#df_long = stack(da, Not(1:3))
df_long = stack(da[:,4:9])

# Create the boxplots
# you can set the y-axis to the same scale for all 
# subplots by using the `link` argument in the `boxplot` function
@df df_long boxplot(:variable, :value, 
    group = :variable, 
    layout = (1, length(unique(df_long.variable))), 
    link = :y, 
    legend = false)


#crop laif to catchment
#lai = wa.mask_trim(laif, poly.shapes[1:end])
lai = wa.mask_trim(laif, g.geometry[1:end])
#replace_missing!(lai, 0)
plot(lai[Ti=15])
lai[X=Near(10),Y=Near(10)]|>plot
#mean(lai)
dfe = ncdf(lai)|>dropmissing
dfp(dfe)

plot(lai[Ti=231])
plot!(g.geometry, fillcolor=false)

#xm = mask(laif[Ti=5];with=g.geometry[6])
xm = mask_trim(laif[Ti=5],g.geometry[6])
plot(xm)
plot!(g.geometry[6], fillcolor=false)



##again
using GeoDataFrames
using GeoInterface
const gdf = GeoDataFrames
n="D:/Wasim/regio/rcm200/v12/catchment.shp"
g = gdf.read(n)
plot(g.geometry, fillcolor=false)
pol = ArchGDAL.IGeometry[]
for i in 1:size(g,1)
    polygon = g.geometry[i]  # assuming the first geometry is a polygon
    src = EPSG(25832)
    dst = EPSG(4326)
    reprojected_polygon = ArchGDAL.reproject(
        polygon, src, dst) #;order=:compliant
    push!(pol,reprojected_polygon)
end
pol
crds = GeoInterface.coordinates.(pol)
reversed_crds = [[[(lon, lat) for (lat, lon) in ring] for ring in polygon] for polygon in crds]
# Create a new polygon with the reversed coordinates
#reversed_polygon = [ArchGDAL.createpolygon(c) for c in reversed_crds]
ArchGDAL.createpolygon(reversed_crds[1])



#reversed_polygon = typeof(pol)[]
# reversed_polygon = []
# for c in enumerate(reversed_crds)
#     push!(reversed_polygon,ArchGDAL.createpolygon(c))
# end
k = wa.reverse_coords(pol)
g.geometry .= 

polygon = pol

GeoInterface.coordinates(polygon)

function rcdf(x::DataFrame)
    polys = []
    for g in x.geometry
        a = g|>first|>GeoInterface.coordinates|>only
        #rc = [(lon, lat) for (lat, lon) in a]
        rc = map(reverse,a)
        pg = ArchGDAL.createpolygon(rc)
        push!(polys,pg)
    end
    return polys
end

k = rcdf(g)


rc = ArchGDAL.createpolygon.(map(reverse,GeoInterface.coordinates.(g.geometry)))


using ArchGDAL
using GeoInterface

function reverse_coords(geometries)
    reversed_geometries = Vector{ArchGDAL.IGeometry{ArchGDAL.wkbPolygon25D}}(undef, length(geometries))

    for (i, geometry) in enumerate(geometries)
        # Extract the coordinates from the geometry
        crds = GeoInterface.coordinates(geometry)

        # Reverse each pair of coordinates
        reversed_crds = [(lon, lat) for (lat, lon) in crds]

        # Flatten the reversed_crds array and convert it to a string
        fc = join([join("$lon $lat" for (lon, lat) in ring) for ring in reversed_crds], ", ")
        # Create the WKT string
        s = "POLYGON (($fc))"
        # Create a new geometry from the WKT string
        geometry = ArchGDAL.fromWKT(s)        
        # Create a new geometry from the WKT representation
        reversed_geometries[i] = geometry
    end

    return reversed_geometries
end

# Use the function on your geometries
rc = reverse_coords(g.geometry)

plot(rc)


s = "POLYGON ((" * join(["$lon $lat" for (lon, lat) in reversed_crds], ", ") * "))"
ArchGDAL.fromWKT(s)
ArchGDAL.createpolygon(rc[1])


pol = gdf.read("D:/Wasim/regio/rcm200/v9/ezp2.shp")
pol.geometry .= reverse_coords(pol.geometry)
fl = "D:/Wasim/Tanalys/DEM/brend_fab/out/m8/evapfab.2016.nc"
r = readras(fl)
plot(r)
plot!(pol.geometry, fillcolor=false)



function rc2(polygon)
    # Extract the coordinates from the polygon
    crds = GeoInterface.coordinates.(polygon)

    # Reverse each pair of coordinates
    reversed_crds = [[(lon, lat) for (lat, lon) in ring] for ring in crds]

    # Create a new polygon with the reversed coordinates
    reversed_polygon = Vector{ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}}(undef, length(reversed_crds))
    for i in 1:length(reversed_crds)
        reversed_polygon[i] = ArchGDAL.createpolygon(reversed_crds[i])
    end

    return reversed_polygon
end

rc2(pol.geometry)

wa.@pyjl
wa.pyjl.xrfacets(fl)
wa.pyjl.xrp(fl)
xrp(fl)

using RCall
@rimport terra as tr
fn=("D:/Wasim/regio/rcm200/v9/ezp2.shp")
r = tr.vect(fn)
R"""file.exists("C:/Program Files/R/R-4.2.2/library/stats/libs/x64/stats.dll")"""

