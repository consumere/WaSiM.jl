# using ArchGDAL
# using Plots

# # read in geojson file
# geojson_file="/mnt/d/ClimateExplorer/pressure/stations.geojson"
# geojson_layer = ArchGDAL.read(geojson_file)
# lyr = ArchGDAL.getlayer(geojson_layer)

# # plot map
# plot()
using ArchGDAL
using Plots
geojson_file="/mnt/d/ClimateExplorer/pressure/stations.geojson"
dataset = ArchGDAL.read(geojson_file) # read the geojson file
layer = ArchGDAL.getlayer(dataset,0) # get the first layer

ArchGDAL.iterate(layer) do feature # iterate over features
    # do something with feature
end

plot() # create an empty plot
ArchGDAL.iterate(layer) do feature # iterate over features
    geom = ArchGDAL.getgeom(feature) # get the geometry
    coords = ArchGDAL.getcoords(geom) # get the coordinates
    plot!(coords,fill = true) # plot the polygon and fill it
end
display(plot()) 


layer = ArchGDAL.getlayer(dataset,0) # get the first layer
ArchGDAL.iterate(layer) do feature # iterate over features
    # do something with feature
end

ds = ArchGDAL.iterate(layer) |>first

#typeof(ds)
ArchGDAL.iterate(layer) do feature # iterate over features
    geom = ArchGDAL.getgeom(feature) # get the geometry
    coords = ArchGDAL.getcoords(geom) # get the coordinates
    plot!(coords,fill = true) # plot the polygon and fill it
end

ArchGDAL.iterate(layer) do feature # iterate over features
    ArchGDAL.getx(feature)
end

ArchGDAL.getx(geom,i)
ArchGDAL.gety(geom,i)
ArchGDAL.getz(geom,i)
ArchGDAL.getpoint(geom,i)
ArchGDAL.getgeom(,1)

point = ArchGDAL.createpoint(1.0,2.0,3.0) # create a point geometry with x = 1.0,y = 2.0,z = 3.0
ArchGDAL.getz(point,0)
    
# Extract coordinates
GeoDataFrames.getgeometrycolumns(gdf)
coords = coordinates(gdf)

# Convert to dataframe
df = DataFrame(coords)

using GMT
pts = GMT.read(fn)


# using GeoJSON
# using DataFrames
# geojson_file="D:/ClimateExplorer/pressure/stations.geojson"
# jsonbytes = read(geojson_file) # read the geojson file as bytes
# fc = GeoJSON.read(jsonbytes) # parse the geojson file as a FeatureCollection object
# df = DataFrame(fc)

# using GeometryBasics

using GeoJSON
using DataFrames
geojson_file="/mnt/d/ClimateExplorer/pressure/stations.geojson"
jsonbytes = read(geojson_file) # read the geojson file as bytes
fc = GeoJSON.read(jsonbytes) # parse the geojson file as a FeatureCollection object
GeoJSON.coordinates(fc.geometry[1])
GeoJSON.coordinates(fc[1]) #|>Matrix
#using JSON3,DataFrames
json_array=GeoJSON.coordinates(fc[1])
# convert to DataFrame
df = DataFrame(x = json_array[1,:],y = json_array[2,:])


using GeoJSON
using DataFrames
geojson_file="/mnt/d/ClimateExplorer/pressure/stations.geojson"
jsonbytes = read(geojson_file) # read the geojson file as bytes
fc = GeoJSON.read(jsonbytes)
tmp=[]
df=[]
xc = []
yc = []
for x in eachindex(fc)
    arr = GeoJSON.coordinates(fc[x])
    #tmp = DataFrame(x = arr[1,:],y = arr[2,:])
    push!(xc,arr[1,:])
    push!(yc,arr[2,:])
end
DataFrame(hcat(xc,yc),:auto)
DataFrame([xc,yc],:auto)
DataFrame("X"=>xc,"Y"=>yc)

Plots.plot(Point(xc,yc))

using GeoJSON
using DataFrames
using GeometryBasics
geojson_file="/mnt/d/ClimateExplorer/pressure/stations.geojson"
jsonbytes = read(geojson_file) # read the geojson file as bytes
fc = GeoJSON.read(jsonbytes)
#cr = hcat([coordinates(x) for x in fc.geometry]...)
cr = DataFrame(hcat([(x) for x in fc.geometry]...),:auto)|>permutedims
nm = DataFrame(hcat([(x) for x in fc.Names]...),:auto)|>permutedims
ht = DataFrame(hcat([(x) for x in fc.ht]...),:auto)|>permutedims
df = hcat(nm,ht,cr,makeunique=true)
rename!(df,["Location","Height","X","Y"])


df|>typeof
# concatenate vectors of each element into matrix
matrix = hcat([vec(elem) for elem in df]...)
# create dataframe from matrix
df = DataFrame(matrix)

dcd = hcat(df)
dcd = reduce((left,right) -> 
      innerjoin(left,right,on = :x,makeunique=true),
      df)


for x in eachindex(fc)
    println(x)
end









#matrix = hcat([vec(elem) for elem in df]...)


geom = Point()
# Extract coordinates
coords = coordinates(geom)
# Print the coordinates
println(coords)

#maby https://github.com/JuliaGeo/LibGEOS.jl





function sday(m3s::Float64,area_km2::Float64)
    q_liters_sec = m3s * 1000.0  # Convert m^3/s to L/s
    q_liters_day = q_liters_sec * 86400  # Convert L/s to m3_day
    q_liters_sec_km2 = q_liters_day / area_km2  # Calculate specific discharge in L/s/km^2
    return q_liters_sec_km2
end

sday(2330.0,198735.0)

function specific_discharge(m3s::Float64,area_km2::Float64)
    q_liters_sec = m3s * 1000.0  # Convert m^3/s to L/s
    q_liters_sec_km2 = q_liters_sec / area_km2  # Calculate specific discharge in L/s/km^2
    return q_liters_sec_km2
end
specific_discharge(2330.0,198735.0)

pyq 2330.0 198735.0  
pysp 2330.0 198735.0  


#Leaflet.jl
#import Pkg; Pkg.add("Leaflet")
using Leaflet
using GeoJSON

# Load GeoJSON file
#geojson_data = GeoJSON.download("https://raw.githubusercontent.com/johan/world.geo.json/master/countries.geo.json")
geojson_data="/mnt/d/ClimateExplorer/pressure/stations.geojson"

fc = GeoJSON.read(read(geojson_data))

# Create Leaflet map
const L=Leaflet
#map = Leaflet.Map(provider = L.Provider.Stamen.TonerLite())
provider = L.OpenTopoMap()
# Add GeoJSON layer to map
geojson_layer = L.geojson(geojson_data)
addLayers!(map,geojson_layer)

# Set map view to focus on the GeoJSON data
fitBounds!(map,geojson_data)

map = Leaflet.Map(;)
# Display map
display(map)

#m = Leaflet.Map(; layers,provider,zoom=3,height=1000,center=[30.0,120.0]);

using EarthEngine,Leaflet

# Initialize EarthEngine
EarthEngine.Initialize()

# Set up an earthengine map
imcollection = filterDate(EarthEngine.ImageCollection("MODIS/006/MOD10A1"),"2021-07-15","2021-07-31")
image = mean(select(imcollection,"Snow_Albedo_Daily_Tile"))
blue_fluorite = ["#291b32","#2a1b34","#2b1b34","#2d1c36","#2f1c38","#301c39","#301d3a","#321d3b","#331d3d","#351d3f","#351e40","#371e41","#381e43","#3a1e45","#3b1f45","#3c1f46","#3e1f48","#3f1f4a","#401f4c","#42204d","#43204e","#44204f","#462051","#472052","#482054","#4a2056","#4a2157","#4c2158","#4e215a","#4f215b","#50215d","#52215e","#532160","#552162","#552263","#562264","#582265","#592267","#5b2268","#5c226b","#5e226c","#5f226e","#60226f","#622271","#632272","#642274","#662276","#672277","#692278","#6a227a","#6c227b","#6e227d","#6e237e","#6f247f","#702480","#712581","#722681","#732683","#742783","#752884","#762985","#772987","#792a87","#792b88","#7a2c89","#7b2c8a","#7c2d8a","#7d2d8c","#7e2e8d","#7f2f8d","#80308e","#813190","#823191","#833292","#843292","#863393","#863494","#873595","#893596","#8a3697","#8b3798","#8b3899","#8c389a","#8e399b","#8e3a9c","#8f3b9c","#8f3d9d","#8f3e9e","#903f9e","#90419e","#90439f","#9044a0","#9046a0","#9047a1","#9049a1","#914aa2","#914ca2","#914ca3","#914ea3","#9150a4","#9151a5","#9153a5","#9154a6","#9156a6","#9157a7","#9258a7","#9259a8","#925aa8","#925ba9","#925da9","#925faa","#9260ab","#9260ab","#9263ac","#9264ac","#9265ad","#9266ae","#9268ae","#9269ae","#926aaf","#926bb0","#926cb0","#926eb1","#926fb1","#9270b2","#9271b2","#9273b3","#9274b3","#9275b4","#9277b5","#9277b5","#9278b6","#927ab6","#927bb7","#927cb7","#927eb8","#927fb8","#9280b9","#9281ba","#9282ba","#9284bb","#9285bb","#9285bc","#9187bc","#9188bd","#918abd","#918bbe","#918cbf","#918dbf","#918ec0","#918fc0","#9191c1","#9092c2","#9094c2","#9094c2","#9095c3","#9096c3","#8f99c4","#8f9ac5","#8f9ac5","#8f9bc6","#8f9cc6","#8f9dc7","#8e9fc8","#8ea0c8","#8ea2c9","#8ea3c9","#8da5ca","#8da5ca","#8da6cb","#8da7cb","#8ca9cc","#8caacc","#8caccd","#8bacce","#8badce","#8baecf","#8ab0d0","#8ab2d0","#8ab2d1","#8ab4d1","#89b4d1","#89b5d2","#89b7d2","#88b8d3","#88bad4","#87bad4","#87bbd5","#86bdd6","#86bed6","#86c0d7","#85c0d7","#85c1d8","#84c3d8","#84c4d9","#83c5d9","#83c6da","#82c8da","#82c8db","#81cadc","#81cbdc","#80ccdd","#81cddd","#84cfdd","#85cfdd","#87d0dd","#8ad0de","#8dd1de","#8fd2de","#90d2de","#92d4de","#95d5de","#97d5de","#98d6de","#9bd7de","#9dd7df","#a0d8df","#a1d9df","#a2dadf","#a5dadf","#a7dbdf","#aadcdf","#abdddf","#acdde0","#afdfe0","#b1dfe0","#b3e0e0","#b4e1e0","#b7e2e0","#bae2e1","#bae3e1","#bee3e2","#c0e4e3","#c1e5e3","#c4e6e3","#c6e6e4","#c8e7e4","#cbe7e5","#cde8e5","#cee9e6","#d2e9e7","#d3eae7","#d5eae7","#d8ebe8","#d9ece8","#dcece9","#deedea","#dfeeea","#e2eeea","#e5efeb","#e6f0eb","#e9f0ec","#ebf1ed","#ecf2ed","#eff3ee","#f1f3ee"]
visParams = Dict(
    :min => 0,
    :max => 100,
    :palette => blue_fluorite,
)
map_id_dict = ee.Image(image).getMapId(visParams)
map_url = map_id_dict["tile_fetcher"].url_format

# Define a leaflet provider
ee_provider = Leaflet.Provider(
    map_url,
    Dict{Symbol,Any}(
        :maxZoom => 20,
        :attribution => """&copy; Openstreetmap France | {attribution.OpenStreetMap}"""
    ),
)

# And the map
map = Leaflet.Map(; provider=ee_provider,zoom=2)

# Notebooks may needs this to work.
WebIO.render(map)




AG.getcoords(pts[1])    
AG.getspatialref(pts[1])
pts|>first|>AG.getcoorddim
pts|>first|>AG.getpoint
pts|>first|>DataFrame
using Plots

x = pts|>first
AG.getpoint(x,0)

# create DataFrame
df = DataFrame(x=[x],y=[y])
plot(df.x,df.y,seriestype=:scatter)



    cr = DataFrame(hcat([(x) for x in fc.geometry]...),:auto)|>permutedims
    nm = DataFrame(hcat([(x) for x in fc.Names]...),:auto)|>permutedims
    ht = DataFrame(hcat([(x) for x in fc.ht]...),:auto)|>permutedims
    df = hcat(nm,ht,cr,makeunique=true)
    rename!(df,["Location","Height","X","Y"])
    metadata!(df,"filename",jsonfile,style=:note);
    return(df)
end



typeof(m)
typeof(vec_int)
parent(vec_int)


pt="/mnt/d/ClimateExplorer/tmean/tmean-stations.geojson"
nd = jread(pt)
DataFrames.metadata(nd)
nd

using ArchGDAL,GeoFormatTypes
const AG = ArchGDAL
myv = [nd[!,:X],nd[!,:Y]]
vec_int = map(x -> map(y -> round(y;sigdigits=0),x),myv)
vec_int = map(x -> map(y -> round(Int64,y),x),myv)
#vec_int = map(x -> map(y -> convert(Float32,y),x),myv)
#vec_int = map(x -> map(y -> convert(Int64,y),x),myv)
ArchGDAL.reproject(myv ,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
ArchGDAL.reproject(vec_int ,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
ArchGDAL.reproject(m,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))

# pts = [(x) for x in fc.geometry]
# using GeometryBasics
# GeometryBasics.Point(pts[1])
# ArchGDAL.reproject(pts[1],EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))

pts = AG.createpoint(nd[1,:X],nd[1,:Y])
ArchGDAL.reproject(pts,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))


using GeoJSON
using DataFrames


import ArchGDAL
dataset = ArchGDAL.read("points.shp")
layer = ArchGDAL.getlayer(dataset,0)
ArchGDAL.iterate(layer) |>first
ArchGDAL.iterate(layer) |>last

###reprojection works on shapes#########
import GeoFormatTypes as GFT
df = GDF.read("points.shp")
df.geometry = GDF.reproject(df.geometry,GFT.EPSG(4326),
    GFT.EPSG(25832))


df.geometry[1]


x=pt
function read_until_flag(file::IOStream,flag::String)
    line = readuntil(file,flag)
    return line[1:end-length(flag)]
end
file = open(x,"r")
while !eof(file)
    line = read_until_flag(file,"timeoff") 
    println(line)
end
close(file)


pts = GDF.read(gauge)
dfd = pts.data
plot(pts.geometry)

pts = GDF.read(gauge)
st = pts.geometry[1]
const AG = ArchGDAL
wkt = AG.toWKT(st)
# Convert geometry to WKT
geometry = pts.geometry[1]
wkt = ArchGDAL.toWKT(geometry)
# Remove Z-coordinate from WKT
wkt = replace(wkt," 0"=> "")
# Convert modified WKT back to geometry
geom = AG.fromWKT(wkt)
# Assign the new geometry back to the desired row
#replace(pts.geometry[1] => geom)




import GeoDataFrames as GDF
import ArchGDAL
pts = GDF.read(gauge)
pts.geometry = GDF.reproject(pts.geometry,GFT.EPSG(4326),
    GFT.EPSG(25832))
#new_geometries = []::Vector{IGeometry{wkbPoint}}
new_geometries = []
for i in 1:length(pts.geometry)
    geometry = pts.geometry[i]
    wkt = ArchGDAL.toWKT(geometry)
    wkt = replace(wkt," 0"=> "")
    new_geometry = ArchGDAL.fromWKT(wkt)
    push!(new_geometries,new_geometry)
end

# Create a new GeoDataFrame with XY coordinates
new_pts = GDF.DataFrame(geometry=new_geometries)
# Copy other columns from the original GeoDataFrame
for col in propertynames(pts)
    if col != :geometry
        new_pts[!,col] = pts[!,col]
    end
end
df = GDF.DataFrame(new_pts)
plot(df.geometry)
plot(pts.geometry)

df.geometry = convert(typeof(pts.geometry),df.geometry)

output_file = "output.geojson"
# Write the GeoDataFrame to a GeoJSON file
GDF.write(output_file,df)

@vv "geojson"

# kd = jread(output_file)
# kd = agjson(output_file)
#kd = GDF.read("points.shp")
kd = GDF.read(output_file)
plot(kd.geometry)
for i in 1:size(kd,1)
    #Plots.annotate!(kd.geometry[i],text(kd.Name[i],:center))
    Plots.annotate!(i,text(kd.Name[i],:center))
end



using Printf
fp = raw"D:\ClimateExplorer\precipitation\precipitation_climexp2509005.kml"
using GeoDataFrames
using Plots
plotlyjs()
kml_content = readlines(fp)
# Regular expression pattern to match coordinates in the KML file
coord_pattern = r"[-+]?\d+\.\d+,\s*[-+]?\d+\.\d+"
# Extract coordinates using regular expression
coordinates = filter(x -> occursin(coord_pattern,x),kml_content)
nam_pattern = r"<name>.*</name>"
nms = filter(x -> occursin(nam_pattern,x),kml_content)

stripped_strings = [strip(replace(s,
"</Name>"=>"",
r"\(DE\).*." => "",
"<name>"=>"")) for s in nms]

#last(stripped_strings,10)
formatted_strings = [titlecase(lowercase(s)) for s in stripped_strings[2:end]]
formatted_strings = [strip(replace(s,"</Name>"=>"")) for s in formatted_strings]
show(formatted_strings)
# Create an empty list to store points
points = []
# Process the coordinates and add them to the points list
for coord in coordinates
    parts = split(coord,",")
    x_coord = parse(Float64,(parts[1]))
    y_coord = parse(Float64,(parts[2]))
    push!(points,(x_coord,y_coord))
end

# Create a GeoDataFrame with the points
#using DataFrames
#df = DataFrame(X = [x for (x,_) in points],Y = [y for (_,y) in points])
gdf = DataFrame(X = [x for (x,_) in points],Y = [y for (_,y) in points],
    name = formatted_strings)
#gdf = hcat(df,formatted_strings)
#gdf = GeoDataFrame(df,crs="+proj=longlat +datum=WGS84")
# #add names
# nm = GeoDataFrames.read(fp)
# stripped_strings = [strip(replace(s,r"\(DE\)$" => "")) for s in nm[!,2]]
# #lowercase.(stripped_strings)
# formatted_strings = [titlecase(lowercase(s)) for s in stripped_strings]
# gdf.name=formatted_strings

Plots.scatter(
    gdf[!,:X],
    gdf[!,:Y],
    color=:blue,
    #label=gdf[!,:name],
    legend=false,
    title=basename(fp),
);

annotate!([(gdf[!,:X][i],gdf[!,:Y][i],
    text(gdf[i,:name],8,:bottom)) 
    for i in 1:size(gdf,1)])



using PlotlyJS




using ArchGDAL
using Plots
geojson_file="D:/ClimateExplorer/precipitation/pre_ogr.geojson"
dataset = ArchGDAL.read(geojson_file) # read the geojson file
layer = ArchGDAL.getlayer(dataset,0) # get the first layer
const AG = ArchGDAL

using GeoDataFrames
#geom = fc.geometry[1]
fc = GeoDataFrames.read(geojson_file)
pts=[]
for geom in fc.geometry
    wkt = ArchGDAL.toWKT(geom)
    # Remove Z-coordinate from WKT
    wkt = replace(wkt," 0"=> "")
    # Convert modified WKT back to geometry
    pt = AG.fromWKT(wkt)
    push!(pts,pt)
end

typeof(pts)
fc.geometry = convert(pts,typeof(fc.geometry))
plot(fc.geometry)


import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL
geojson_file="D:/ClimateExplorer/precipitation/pre_ogr.geojson"
pts = GDF.read(geojson_file)
pts.geometry = GDF.reproject(pts.geometry,GFT.EPSG(4326),
    GFT.EPSG(25832))

new_geometries = []
for i in 1:length(pts.geometry)
    geometry = pts.geometry[i]
    wkt = ArchGDAL.toWKT(geometry)
    wkt = replace(wkt," 0"=> "")
    new_geometry = ArchGDAL.fromWKT(wkt)
    push!(new_geometries,new_geometry)
end

# Create a new GeoDataFrame with XY coordinates
new_pts = GDF.DataFrame(geometry=new_geometries)
# Copy other columns from the original GeoDataFrame
for col in propertynames(pts)
    if col != :geometry
        new_pts[!,col] = pts[!,col]
    end
end
df = GDF.DataFrame(new_pts)

plot(pts.geometry)
output_file = "prec_station_25832.geojson"
# Write the GeoDataFrame to a GeoJSON file
of = joinpath(dirname(geojson_file),output_file)
GDF.write(of,df)
kd = GDF.read(of)
plot(kd.geometry)
for i in 1:size(kd,1)
    Plots.annotate!(kd.geometry[i],text(kd.Name[i],:center))
#    Plots.annotate!(i,text(kd.Name[i],:center))
end

kd.Name


#now a reader.
fl = CSV.read("D:/Wasim/Tanalys/DEM/Input_V2/meteo/pre_ce.txt",DataFrame;limit=4)
xc = fl[2,5:end]|>collect
yc = fl[3,5:end]|>collect

g = hcat(xc,yc)

pts=[]
for i in eachindex(g)
    println(i)
end

pts::ArchGDAL.IGeometry{ArchGDAL.wkbPoint} = []
for i in 1:size(g,1)
    #enumerate(g)    
    #pt = AG.createpoint(xc,yc,EPSG(25832))
    pt = AG.createpoint([xc[i],yc[i]])
    #pt = AG.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
    push!(pts,pt)
end

nd = DataFrame(geometry=pts,name=propertynames(fl)[5:end])

typeof(kd.geometry[1])

plot(nd.geometry)



function stprp(fn::String)
    fl = CSV.read(fn,DataFrame;limit=4)
    xc = fl[2,5:end]|>collect
    yc = fl[3,5:end]|>collect
    pts = ArchGDAL.IGeometry[]
    for i in 1:length(xc)
        pt = ArchGDAL.createpoint([xc[i],yc[i]])
        pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
        push!(pts,pt)
    end
    nd = DataFrame(geometry=pts,name=propertynames(fl)[5:end],xc=xc,yc=yc)
    return nd
end

od = stp("D:/Wasim/Tanalys/DEM/Input_V2/meteo/pre_ce.txt")
od = stp("D:/Wasim/Tanalys/DEM/Input_V2/meteo/temp.2021")
od = stprp("D:/Wasim/Tanalys/DEM/Input_V2/meteo/temp.2021")
plot(od.geometry)
for (i,pt) in enumerate(od.geometry)
    x = od.xc[i]
    y = od.yc[i]
    name = od.name[i]
    annotate!(x,y,text(name,8,:black,:center))
end
plot!()


AG.getpoint(od[!,1])

ArchGDAL.getx(od.geometry[1],0)

# for (i,pt) in enumerate(od.geometry)
#     pts = ArchGDAL.createpoint(pt)

#     x,y =  
#     [ArchGDAL.getpoint(x,0)[1] for x in pts],
#     [ArchGDAL.getpoint(x,0)[2] for x in pts]

#     name = od.name[i]
#     annotate!(x,y,text(name,8,:black,:center))
# end
   
#ArchGDAL.getpoint(od.geometry[5])[2]

od = stp("D:/Wasim/Tanalys/DEM/Input_V2/meteo/spec_1970_2017_pt2.txt")
plot(od.geometry)
for (i,pt) in enumerate(od.geometry)
    x = od.xc[i]
    y = od.yc[i]
    name = od.name[i]
    annotate!(x,y,text(name,8,:black,:center))
end
plot!()


wa.stplot("D:/Wasim/Tanalys/DEM/Input_V2/meteo/temp.2021")



###

#try to rasterize in julia
using Rasters,ArchGDAL,Plots,Dates,Shapefile
using Rasters.LookupArrays

input = 'D:/Wasim/main/db_table.csv'
shapefile_name = "D:/Bodendaten/buek200_2020/merged_buek200.shp"
bks = Shapefile.Handle(shapefile_name).shapes[1:end]

#reducer: a reducing function to reduce the fill value for all geometries that cover or touch a pixel down to a
#single value. The default is last. Any that takes an iterable and returns a single value will work,including
#custom functions. However,there are optimisations for built-in methods including sum,first,last,minimum,
#maximum,extrema and Statistics.mean
bks[1]|>plot
bks[end]|>plot

rs = rasterize(first,bks[1]; 
    shape = :polygon,
    res=60,missingval=0,fill=1,
    boundary=:touches,progress=true)


d = bks[end]
d.parts

using GeoDataFrames
bks = GeoDataFrames.read(shapefile_name)
select!(bks,[:geometry,:Symbol])

rs = rasterize(first,bks; 
    shape = :polygon,
    res=60,missingval=0,fill=1,
    boundary=:touches,progress=true)

plot(rs)
Rasters.values(rs)

heatmap(bks.geometry,bks.Symbol)

plot(bks.geometry[1:5])
@vv "reproject"

# "reprojection of shapes"
# note: this is inplace!
using ArchGDAL
const AG = ArchGDAL
AG.reproject(bks.geometry,EPSG(25832),
    ProjString("+proj=longlat +datum=WGS84 +no_defs"))

# #pts = typeof(bks.geometry)[]
# for geom in bks.geometry
#     pt = AG.reproject(geom,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
# #    push!(pts,pt)
# end

plot(bks.geometry[end])
plot(bks.geometry[4:8])
plot(bks.geometry[50:end])

typeof(bks.geometry)
typeof(pts)



raw"D:\Wasim\regio\out\rc200\x2n"|>cd

x = @gl "sb05"
r = Raster(x)
flags = Dict(:tr => [1000,1000],:r => :near)
#r2 = Rasters.warp(r[t=1],flags)
Rasters.resample|>methods
r2 = Rasters.resample(r[t=1];res=tuple(1000,1000))

tuple(1000,1000)|>typeof

r = readras(x)
r2 = Rasters.resample(r[t=1];res=tuple(1000,1000))

const AG = ArchGDAL
r = AG.read(x)
f1 = AG.getlayer(r,0)
pt = AG.reproject(f1,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs")) 


using Rasters,ArchGDAL,Plots,Dates,Shapefile
shapefile_name = "D:/Bodendaten/buek200_2020/merged_buek200.shp"
bks = Shapefile.Handle(shapefile_name).shapes[1:end]
rs = rasterize(first,bks[1]; 
    shape = :polygon,
    res=60,missingval=0,fill=1,
    boundary=:touches,progress=true)
plot(rs)



using ArchGDAL
const AG = ArchGDAL
pt="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/wasim_pest/app/input/str_smoothed.tif"
# Open the raster dataset (replace 'x' with your raster file path)
r = AG.readraster(pt)

# Check if the raster dataset was successfully opened
if r isa ArchGDAL.RasterDataset
    # Select the band you want to reproject (e.g.,band 1)
    band_index = 1
    f1 = AG.getband(r,band_index)
    #f1 = AG.getlayer(r)
    
    # Define the target spatial reference (e.g.,EPSG 25832)
    #target_srs = EPSG(25832)
    target_srs = EPSG(4326)
    
    # Reproject the raster band to the target spatial reference
    pt = AG.reproject(f1,target_srs)
    
    # Save the reprojected raster to a new file (optional)
    # AG.write("path/to/output/reprojected_raster.tif",pt)
    
    println("Raster reprojected successfully.")
else
    println("Failed to open the raster dataset.")
end



function reproject_netcdf(input_file::AbstractString,
    output_file::AbstractString,
    source_srs::AbstractString,target_srs::AbstractString)
    cmd = `cmd.exe /c gdalwarp -s_srs $source_srs -t_srs $source_srs -of NetCDF -dstfile $output_file $input_file`
    run(pipeline(cmd))
end

x
reproject_netcdf(x,"test.nc","EPSG:25832","EPSG:4326")

cmd = `cmd.exe /c gdalwarp -s_srs EPSG:25832 -t_srs EPSG:4326 -of NetCDF $x test.nc`
run(pipeline(cmd))
wa.pw()

cp(x,raw"C:/OSGeo4W64/sb.nc")

using Conda
#Conda.add("gdal")
Conda._get_conda_env()

using RCall
function repro_tr(   
    input_file::AbstractString,
    target_srs::Int64)
    tr = rimport("terra")
    rin = tr.rast(input_file,drivers="NETCDF")
    @rput rin
    R"""terra::crs(rin) <- 'EPSG:25832'"""
    n = rcopy(R"rin")
    #tr.crs(rin) = "EPSG:25832"
    # R"
    # require(terra)
    # rin = rast($input_file,drivers='NETCDF')
    # terra::crs(rin) = 'EPSG:25832'
    # ot = project(n,crs='EPSG:$target_srs')
    # "
    # #out = tr.project(in,"EPSG:"*string(target_srs))
    #out = rcopy(R"ot")
    out = tr.project(n,"EPSG:"*string(target_srs))
    return out
end

z = repro_tr(x,4326) #RObj


function repro_tr(   
    input_file::AbstractString,
    output_file::AbstractString,
    target_srs::Int64)
    tr = rimport("terra")
    rin = tr.rast(input_file,drivers="NETCDF")
    # tr.NAflag(rin)=-9999
    # rin = tr.flip(tr.trans(
    #     tr.rev(rin)),
    #     direction = "vertical")
    
    @rput rin

    R"""
    terra::crs(rin) <- 'EPSG:25832'
    terra::NAflag(rin) <- -9999

    rin <- terra::flip(terra::trans(
    terra::rev(rin)),direction = 'vertical')
    
    """
    # @rput rin
    # R"""terra::crs(rin) <- 'EPSG:25832'"""
    
    n = rcopy(R"rin")
    out = tr.project(n,"EPSG:"*string(target_srs))
    tr.writeCDF(out,output_file,overwrite=true)
    #return out
    println("$output_file done!")
end

repro_tr(x,"ot.nc",4326)
rx=Raster("ot.nc")
rx=Raster(x)
plot(rx)


function nctc(   
    input_file::AbstractString,
    output_file::AbstractString,
    target_srs::Int64;mask = true)
    #@rput mask=T
    
    if mask
        rmask = "TRUE"
    else
        rmask = "FALSE"
    end

    R"""
    rin = terra::rast($input_file)
    terra::crs(rin) <- 'EPSG:25832'
    terra::NAflag(rin) <- -9999

    rin <- terra::flip(terra::trans(
        terra::rev(rin)),direction = 'vertical')
    #rin[rin <= 0] <- NA #masking
    mask = as.logical($rmask)
    ifelse(mask,yes = {
        rin[rin <= 0] <- NA
},no = rin)

    out = terra::project(rin,paste0("EPSG:",$target_srs))
    writeCDF(out,$output_file,overwrite=T)
    """
    return Raster(output_file)
    println("$output_file done!")
end

z = nctc(x,"ot2.nc",4326)
plot(z)
rpr(z)
agheat("ot2.nc";step=120)
agheat(x;step=100)

flags = Dict(:tr => [1000,1000],:r => :near)
r2 = warp(z,flags)


pwd()
x = @gl "sb1"
xr=pyimport("xarray")
plt=pyimport("matplotlib.pyplot")
ds = xr.open_dataset(x)
s = pystr(ds.keys()) 
tv = split(s)[end-5]
#plt.plot(ds[tv].isel(t=0)) #wrong
ds[tv].isel(t=0).transpose().plot() #right
plt.show() #show the plot

x="tempsinn100.2017.nc"
dx = xr.open_dataarray(x)
dx.values|>typeof #das ist julia array.

maskval = float(5)
# dx = ds[tv].isel(t=0).transpose()
# dx.where(dx.values>maskval).plot(cmap="cividis")
xr.open_dataarray(x,mask_and_scale=true).transpose().plot(cmap="cividis")
plt.show()

x=nconly("va")|>last
py"""
from xarray import open_dataset
from matplotlib.pyplot import show
maskval = float(0)
dx = open_dataset($x,mask_and_scale=True).isel(t=0).transpose().to_array()
dx.where(dx.values>maskval).plot(cmap="cividis")
show()
"""

function xrplot(x::AbstractString;maskval=0,lyr=0)
    x = nconly(x)|>last
    py"""
    from xarray import open_dataset
    from matplotlib.pyplot import show
    maskval = float($maskval)
    dx = open_dataset($x,mask_and_scale=True).isel(t=$lyr).transpose().to_array()
    dx.where(dx.values>maskval).plot(cmap="cividis")
    show()
    """
end

xrplot("wind";maskval=7.25)


x = nconly(".")|>last

function frx(x::AbstractString;maskval=0,lyr=0)
    x = nconly(x)|>last
py"""
ncf = $x;
ncf = ncf if ncf.endswith('.nc') else print('not an NetCDF grid!\ntry rioplot',ncf) & exit(0)
lyr=$lyr
maskval = float($maskval)
import xarray as xr
from matplotlib import pyplot
import numpy as np
from pandas import DataFrame
dx=xr.open_dataset(ncf)
m=list(dx.variables)[-1]
print('cf attrs:',dx[m].attrs)
df=dx[m].quantile(q=[0,.05,.1,.25,.5,.75,.9,.95,1]).to_dataframe()
print('fieldmean:',dx[m].mean().to_numpy())
print(list(dx.variables),'\n',df,'\n','-'*18)
idx = list(dx.indexes._coord_name_id)
lng = ''.join([i for i in idx if (i.startswith('lon') or i.endswith('x'))])
lat = ''.join([i for i in idx if (i.startswith('lat') or i.endswith('y'))])
tim = ''.join([i for i in idx if (i.startswith('t'))])
if idx == ['x','y','t']:
    print('this is a wasim grid... transposing now...','\n','-'*18)
    dx = dx.transpose().to_array().isel(t=lyr)
elif idx == ['x','y',tim]:
    dx = dx.rename_dims({tim:'t'}).transpose('x','y','t')
    dx = dx.to_array().isel(t=lyr)
else:
    dx = dx.rename_dims({tim:'t',lng:'x',lat:'y'})
    print('renamed coord_name_ids: ',idx,'to',dx.dims,'!!!')
    dx = dx.isel(t=lyr).to_array()
y=dx.where(dx.values>maskval)
y=y if not np.all(np.isnan(y.data)) else print('all data of selection is NaN\n exiting now from',ncf) & exit(0)
print('fieldmean of selection:',y.mean().to_numpy())
df=y.variable.quantile(q=[0,.05,.1,.25,.5,.75,.9,.95,1])
ar=np.array([0,.05,.1,.25,.5,.75,.9,.95,1])
df=DataFrame(data=df,index=ar).reset_index().rename({'index':'q',0:m},axis=1)
print('Quantiles of selection:\n',df,'\n','-'*18)
pyplot.rcParams['figure.figsize']=(10,9)
y.plot(cmap='turbo')
pyplot.tight_layout()
title=(ncf.rsplit('/')[-1]+' Layer: '+str(lyr))
pyplot.title(title)
pyplot.show()
"""
end

#p.fig.savefig("test.png")


frx(x)


Plots.gr()
od = stp("D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt")
wa.fzplot(;dirpath="D:/Wasim/regio/rcm200/v12/")
plot!(od.geometry)
for (i,pt) in enumerate(od.geometry)
    x = od.xc[i]
    y = od.yc[i]
    name = od.name[i]
    annotate!(x,y,text(name,8,:black,:center))
end
plot!()



using JSON

#x="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/saale_s2"
# cmd = `wsl.exe -e "tail +7 $x|cut -f 5| jq -s '{minimum:min,maximum:max,average:(add/length),
# median:(sort|if length%2==1 then.[length/2|floor]else
# [.[length/2-1,length/2]]|add/2 end),sum:add}'" `
#run(cmd)

#qq saale_wolfsmuenster.txt
x = """
{"minimum": 1.71,"maximum": 326,
"average": 17.528997919804915,
"median": 11.1,"sum": 244371.76000000033 
}
"""

JSON.parse(x)|>DataFrame

v="D:/Wasim/Tanalys/DEM/Grundwasser/gw2022/wegfurth.wats"
readdlm(v)



########join two GeoDataFrames ########
using GeoDataFrames
using DataFrames

dfa = GeoDataFrames.read("a.csv")
dfb = GeoDataFrames.read("b.csv")

X = crossjoin(dfa, dfb, makeunique=true)
subset!(X, [:geometry, :geometry_1] => (a, b) -> intersects.(a, b))


#########try to crop a shpfile## ->> see rst.crop_geodf ! ########
S="D:/Wasim/Tanalys/DEM/brend_fab/in3/cir.shp"
using GeoDataFrames
@time setup() #25s
cdof(S)
gd = GeoDataFrames.read(S)
#dropmissing!(gd,5)
plot(gd.geometry)
r = Raster("D:/Wasim/Tanalys/DEM/brend_fab/in3/fab.intern.sl.nc")
@edit @vv "crop"
#rc = Rasters.crop(r,gd.geometry)
#wintree()
"Shapefile"|>wa.vgjl
using Shapefile
sh = Shapefile.Handle(S)

plot(r)
sh.shapes[40:400]|>plot!
sh.shapes[500:3000]|>plot
sh.shapes[end-3000:end]|>plot!

rc = Rasters.crop(r;to=sh.shapes[end])
rc = Rasters.crop(r;to=gd)
plot(rc)

src = Rasters.crop(gd.geometry;to=r)
typeof(gd)|>cb
using ArchGDAL, Rasters
# Load the GeoDataFrame and Raster
geodf = ArchGDAL.read(S)
#raster = Rasters.read("raster.tif")
# Get the extent of the Raster
raster_extent = ArchGDAL.extent(r)
# Crop the GeoDataFrame with the Raster extent
cropped_geodf = ArchGDAL.intersect(geodf, raster_extent)
# Save the cropped GeoDataFrame
ArchGDAL.write("cropped_geodf.geojson", cropped_geodf)

v1 = raster_extent[1]
#v1 = range(v1[1],v1[2])
v1 = [(v1[1],v1[2])]
v2 = raster_extent[2]
v2 = [(v2[1],v2[2])]
Tp = Tuple(v1,v2)
raster_polygon = ArchGDAL.createpolygon()

raster_polygon = ArchGDAL.createpolygon(r)

#ext = extend(tr)
xmin = dims(r, (X, Y))[1]|>first
xmax = dims(r, (X, Y))[1]|>last
ymin = dims(r, (X, Y))[2]|>first
ymax = dims(r, (X, Y))[2]|>last

crop_poly = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]]
Rasters.crop(gd;to=crop_poly)

cropped_geodf = ArchGDAL.intersect(geodf, crop_poly)

sub=Rasters.mask(sh; with=r)
gd[!, :geometry][1]|>plot

gd[!, :geometry][1]|>typeof|>cb

r = replace_missing(r,missing)
wa.vgjl("polygonize")

raster_path = "D:/Wasim/Tanalys/DEM/brend_fab/in3/fab.intern.sl.nc"
r = ArchGDAL.readraster(raster_path)
lyr=1
dx = ArchGDAL.getband(r,lyr)
x = ArchGDAL.getgeotransform(r)[2]
y = ArchGDAL.getgeotransform(r)[6]
#x = Raster(dx,dims=(Y,X))
const AG=ArchGDAL
AG.getgeometry(dx)
ArchGDAL.intersection(g1, g2)

using Conda
Conda.add("geopandas")
using PyCall
gp = pyimport("geopandas")


fn="D:/Wasim/Tanalys/DEM/brend_fab/in3/fab.art-bfid"
cdof(fn)
@time setup() #25.026s

agcont(fn)
agcont2(fn;step=150)
using GeoDataFrames
br = GeoDataFrames.read("fab-bk.shp")
plot(br.geometry[1:end])
for i in 1:length(br)
    annotate!(br.geometry[i].x, br.geometry[i].y, text(br[:Symbol][i], 10))
end
#plot!(br.geometry[1])
step=150
br.geometry|>first
for i in 1:step:size(br,1)
    #plot!(br.geometry[i],color=:black)
    annotate!(br.geometry[i],
    Plots.text(string(br.Symbol[i]), 7, :red, :center, 
    halign=:center, rotation=-35.0))
end

for i in 1:step:size(dx, 1)
    for j in 1:step:size(dx, 2)
        value = round(dx[i, j]; digits=roundto)
        color = isnan(value) ? :white : :black
        annotate!(j, i, Plots.text(string(value), 7, color, :center, 
            halign=:center, rotation=-35.0))
    end
end


od = br
begin
p = plot(od.geometry);
for (i, pt) in enumerate(od.geometry)
    midpt = ArchGDAL.centroid(od.geometry[i])
    x = ArchGDAL.getx(midpt, 0)
    y = ArchGDAL.gety(midpt, 0)
    #name = od.Symbol[i]
    name = od.Legende[i]
    name = split(name,",")[1]
    name = split(name,":")[1]
    name = split(name," ")|>last
    annotate!(x, y, text(name, 8, :black, :bottom, :left))
end
display(p)
end 


###lookup cords.
lk="D:/Wasim/Tanalys/DEM/brend_fab/out/c8/cdx/sb1_fab.2012.nc"
mr = readras(lk)
dims(mr)
DimPoints(mr)

ar = DimPoints(mr)|>collect
ar[:,:,1]
#df = DataFrame(ar[:,1]|>collect,:auto)
df = DataFrame(ar[:,860])
#df = DataFrame(ar[:,:,1],:auto)
rename!(df,["X","Y","val"])
plot(df.val)
lookup(mr,:X,:Y)|>DataFrame
x,y = (lookup(mr,:X)|>collect,lookup(mr,:Y)|>collect)

"""
fn_prefix::AbstractString,
    array,    wkt::AbstractString,
    transform,    file_format::String
"""
function write_raster(fn_prefix::AbstractString,
    array,
    wkt::AbstractString,
    transform,
    file_format::String)
    # transponse array back to columns by rows
    array_t = permutedims(array, [2, 1])

    width, height = size(array_t)

    # Define extension and driver based in file_format
    file_format == "tif" ? (ext = ".tif"; driver = "GTiff") :
    (ext = ".asc"; driver = "AAIGrid")

    file_format == "tif" ? (options = ["COMPRESS=LZW"]) :
            (options = [])

    # Append file extention to filename
    fn = string(fn_prefix, ext)

    # Create raster in memory
    ArchGDAL.create(fn_prefix,
    driver = ArchGDAL.getdriver("MEM"),
    width = width,
    height = height,
    nbands = 1,
    dtype = eltype(array_t),
    options = options) do dataset
    band = ArchGDAL.getband(dataset, 1)
    # Write data to band
    ArchGDAL.write!(band, array_t)

    # Write nodata and projection info
    ArchGDAL.setnodatavalue!(band, -9999.0)
    ArchGDAL.setgeotransform!(dataset, transform)
    ArchGDAL.setproj!(dataset, wkt)
    # Copy memory object to disk (necessary because ArchGDAL.create
    # does not support creation of ASCII rasters)
    ArchGDAL.destroy(ArchGDAL.copy(dataset,
        filename = fn,
        driver = ArchGDAL.getdriver(driver),
        options = options))
    end
    nothing
end

using ArchGDAL,Rasters
#r = Raster("D:/Wasim/Tanalys/DEM/brend_fab/in3/fab.art-bfid")
geojson_file="d:/ClimateExplorer/pressure/stations.geojson"
cdof(geojson_file)
dataset = ArchGDAL.read(geojson_file) # read the geojson file
layer = ArchGDAL.getlayer(dataset,0) # get the first layer

pwd()
const AG = ArchGDAL
image = fill(UInt8(127), (150, 250))
AG.create(
    "example.tif",
    driver = AG.getdriver("GTiff"),
    width = 500,
    height = 300,
    nbands = 1,
    dtype = UInt8) do raster
    AG.write!(
        raster,
        image, # image to "burn" into the raster
        1, # update band 1
        30:180, # along (window) xcoords 30 to 180
        50:300 # along (window) ycoords 50 to 300
    )
end
# import ArchGDAL
r = Raster("example.tif")


cd("d:/Wasim/Tanalys/DEM/brend_fab/out/c8/loc3/spin")
using GDAL_jll
rp=@gl "sb1"
# list information about a raster dataset
run(`$(gdalinfo_path()) $rp`)
# convert raster data between different formats
run(`$(gdal_translate_path()) -of COG input.asc output.tif`)

# list information about an OGR-supported data source
run(`$(ogrinfo_path()) path/to/vector-file`)

@time @time_imports setup() 
#20.201747 seconds nice!
agheat(rp;step=150,roundto=1)
wa.agheat(rp;step=110,roundto=false)
wa.agheat(rp;step=91,roundto=3)
agcont(rp)
agcont2(rp)
wa.agcont2(rp;step=99,roundto=false)
rp=@gl "temper"
wa.agcont2(rp;step=99,roundto=4)


cd("D:/Wasim/regio/rcm200/v12/")
fn = "rcm.art-bfid"
fo = "out.tif"
#run(`$(gdal_translate_path()) -of COG input.asc output.tif`)
run(`$(gdal_translate_path()) -of COG $fn $fo`)
r=Raster(fo)
plot(r)

v = dfonly("sb")
pv = map(cmk.cloudplot2,map(fread,v[1:2]))
pv[1]
pv[2]
cmk.tsp(fread(v[1]))
cmk.tsp(fread(v[2]))
v[1]|>fread|>dfp

#check wit and dep grids
cd("D:/Wasim/regio/rcm200/v12/")
wit = "rcm.wit"
dep = "rcm.dep"
wa.agheat(wit;step=25,roundto=1)
wa.agheat(dep;step=15,roundto=1)
wa.agcont2(dep;step=5,roundto=1)
cdu()
dvec = rglob(r"dep$")
@edit descr(r)
@doc wa.descr
wa.descr(dvec[2];missval=0)
wa.descr(dvec[3])
a = []
for i in dvec
    rs = Raster(i)
    r = replace_missing(rs, 0)
    rd = wa.stats(r)
    rd.name .= i
    #rd = hcat(wa.stats(r),i)
    #push!(a,DataFrame(rd,:auto))
    push!(a,rd)
end
da = reduce(vcat,a)
import CairoMakie as Mke
#Mke.scatter(da[!,:mean],da[!,:std])
Mke.scatter(tovec(da,1),tovec(da,3),
    legend=tovec(da,:name))


import CairoMakie as Mke
# Assuming `da` is your DataFrame and `tovec` is a function that converts a DataFrame column to a vector
x = tovec(da, 1)
y = tovec(da, 3)
#labels = tovec(da, :name)
lab = replace.(String.(da[!,:name]),".\\"=>"","\\"=>"_")

fig = Mke.Figure()
ax = Mke.Axis(fig[1, 1])
# Create a scatter plot with labels
for i in 1:length(x)
    Mke.scatter!(ax, [x[i]], [y[i]], label=lab[i])    
end
Mke.Legend(fig[1, 2], ax)
fig

scatter(da[!,:mean],da[!,:sd])
lab = replace.(String.(da[!,:name]),".\\"=>"","\\"=>"_")
pie(lab,da[!,:sd])

sort(da,:max)
#maxdepth
bar(da[!,:max],legend=false,
    xticks=(1:nrow(da),lab),rotation=45)


xa = cmk.tsp
xa = cmk.tsbar
"D:/Wasim/regio/out/rc200/x3/cl/qbasrcm.x3.2005"|>dfr|>xa
"D:/Wasim/regio/out/rc200/x3/cl/qbasrcm.x3.2005"|>dfr|>cmk.cloudplot2

#using FileIO 
ds = stack(da)
fig = Mke.Figure()
ax = Mke.Axis(fig[1,1])
# Mke.barplot!(ax,
#             ds.variable,
#             ds.value,
#             color=ds.name)


x=raw"L:\04-Phil\Geo1data\prj-Efre-Daten\GIS-Daten@Server K\hydrosheds\hybas_lake_eu_lev07_v1c.shp"
hsd = gdf.read(x)
plot(hsd.geometry[101:end])
n="D:/Wasim/regio/rcm200/v12/catchment.shp"
g = gdf.read(n)
#xm = mask_trim(ras,g.geometry)

names(hsd)
@rsubset hsd



vgjl("dep ")
"D:/Wasim/regio/rcm200/"|>cd
v = rglob(r"dep$")
a=[]
@edit descr(r)
for i in v
    #df = wa.stats(Raster(i,missingval=-9999))
    r = Raster(i) #,missingval=Float64(-9999)
    r = replace_missing(r, 0)
    m = mean(r) # get the mean for each band
    n = minimum(r) # get the minimum for each band
    x = maximum(r) # get the maximum for each band
    d = median(r) # get the median for each band
    s = std(r) # get the standard deviation for each band
    arr=[m,n,x,d,s]'
    df = DataFrame(arr,:auto)
    nm=["mean", "min", "max", "median", "sd"]
    rename!(df,nm)
    df.nm .= i
    push!(a,df)
end
a = reduce(vcat,a)
sort(a,:max,rev=true)
sort!(a,:mean,rev=true)

wd = wa.rstats(r"[.]wit$")



a = Raster(v[2])
plot(a)
describe(a)

@doc wa.agheat(v[1];step=25,msk=1)
wa.agheat(v[1];step=25,msk=1)

wa.agheat(dep;step=15,roundto=1)
wa.agcont2(dep;step=5,roundto=1)


k=raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\pre_genRE\yrsum2014.tif"
agheat(k)



###excel to geojson
import GeoDataFrames as GDF
import XLSX
fn="C:/Users/chs72fw/Desktop/OneDrive - Universität Würzburg/divers/Bodendaten_Zusammengefasst.xlsx"
ojs="C:/Users/chs72fw/Desktop/OneDrive - Universität Würzburg/divers/Bodendaten.geojson"
begin
    df = DataFrame(XLSX.readtable(fn, "Ergebnis"))
    # Create a Point geometry for each row
    geometry = [ArchGDAL.createpoint(
        (df.RECHTSWERT[i],df.HOCHWERT[i]
            )) for i in 1:nrow(df)]
    df.geometry = geometry
    # Convert all columns to if type is Any
    for col in names(df)
        if eltype(df[!, col]) == Any
            try
                df[!, col] = string.(df[!, col]) #parse.(Float64, df[!, col])
            catch e
                println("Failed to convert column $col: $e")
            end
        end
    end
    # Write the GeoDataFrame to a file
    GDF.write(ojs, df) #wrks!
end
#tryparse.(String, df[!, 1])
df


fn=raw"C:\Users\chs72fw\Desktop\OneDrive - Universität Würzburg\divers\GRDC_Stations_selection.xlsx"
ojs="C:/Users/chs72fw/Desktop/OneDrive - Universität Würzburg/divers/grdc.geojson"
sheet = "station_catalogue"

df = DataFrame(XLSX.readtable(fn, sheet))
#names(df)|>println
# Create a Point geometry for each row
gm = [ArchGDAL.createpoint(
    (df.long[i],df.lat[i]
    )) for i in 1:nrow(df)]

df.geometry = gm
dropmissing!(df, :country)
df = filter(x->x.country=="DE",df)

begin
    # Convert all columns to if type is Any
    for col in names(df)
        if eltype(df[!, col]) == Any
            try
                df[!, col] = string.(df[!, col]) #parse.(Float64, df[!, col])
            catch e
                println("Failed to convert column $col: $e")
            end
        end
    end
    # Write the GeoDataFrame to a file
    GDF.write(ojs, df) #wrks!
end



########
import GeoDataFrames as GDF
ojs = "https://cdn.jsdelivr.net/npm/world-atlas@2/countries-10m.json"
gd = GDF.read(ojs)
#filter(row -> occursin(r"G", row[:name]), gd)
@rsubset gd begin
    :name=="Germany"
end
#@subset gd @byrow
filter(row -> occursin(r"Ger", row[:name]), gd).geometry|>
plot

filter(row -> occursin(r"Geor", row[:name]), gd).geometry|>z->
plot(z,fillcolor=:transparent,title="Georgia")

filter(x-> occursin(r"German", x[:name]), gd).geometry|>X->
plot(X,fillcolor=:transparent,title="")

@subset gd @byrow begin
    :name=="Germany"
end

sel = filter(x-> occursin(r"German", x[:name]), gd)
dir="C:/Users/chs72fw/Desktop/OneDrive - Universität Würzburg/divers/"
isdir(dir) ? println("yes") : println("no")
@doc GDF.write(joinpath(dir,"ger.geojson"), sel;layer_name="Germany")
import GeoFormatTypes as GFT
GDF.write(joinpath(dir,"ger.geojson"), 
    sel;layer_name="Germany",
    crs=convert(WellKnownText, EPSG(4326)))
    #crs=EPSG("EPSG:4326"))
    
GDF.write(joinpath(dir,"ger.shp"), sel;layer_name="Germany")
GDF.write("tst.geojson", sel;layer_name="Germany") #error
Shapefile.write("tst.shp",sel) #wrks!
Shapefile.write(joinpath(dir,"ger.shp"),sel) #wrks!


inpath="D:/Wasim/regio/rcm200/v0/"
bsns=gdf.read("D:/Wasim/regio/rcm200/v0/basins-v0.shp")
gdf.write("D:/Wasim/regio/rcm200/v0/basins-v0.geojson",bsns)
gdf.write("D:/Wasim/regio/rcm200/v0/basins-v0.kml",bsns)


z=nread("D:/Wasim/regio/out/rc200/x12/mfile")
z=fread("D:/Wasim/regio/out/rc200/x12/mfile")
z=waread2("D:/Wasim/regio/out/rc200/x12/mfile")
select!(z,Not(Cols(r"Wolfsmuenster-qoutjl_Wolfs")))
heat(z[!,(1:10)])
dfm(z)


#smc = crop(sm.geometry,cmt.geometry)
# Load the GeoDataFrame and Raster
geodf = ArchGDAL.read("D:/Wasim/regio/rcm200/v10/bk200_ezg.shp")
# Get the extent of the Raster
raster_extent = ArchGDAL.extent(r)
r = rst.readras("D:/Wasim/Tanalys/DEM/brend_fab/out/c10/loc/sb1_fab.2015.nc")
xmin = dims(r, (X, Y))[1]|>first
xmax = dims(r, (X, Y))[1]|>last
ymin = dims(r, (X, Y))[2]|>first
ymax = dims(r, (X, Y))[2]|>last
crop_poly = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]]
crop_poly = ArchGDAL.createpolygon(crop_poly)

gdf.within(sm.geometry[end],crop_poly)
gdf.intersection(sm.geometry[end],crop_poly)
new = Vector{ArchGDAL.IGeometry}(undef, length(sm.geometry))
# for i in 1:length(sm.geometry)
#     new[i] = gdf.intersection(sm.geometry[i], crop_poly)
# end
np = [gdf.intersection(geom, crop_poly) for geom in sm.geometry]
npw = filter(x -> !ArchGDAL.isempty(x), np)

plot!(r)



using ArchGDAL
using DataFrames

"""
r = rst.readras(file)
xmin = dims(r, (X, Y))[1]|>first
xmax = dims(r, (X, Y))[1]|>last
ymin = dims(r, (X, Y))[2]|>first
ymax = dims(r, (X, Y))[2]|>last
crop_poly = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]]
crop_poly = ArchGDAL.createpolygon(crop_poly)
xdf = crop_geodf(sm, crop_poly)
edf = reduce(GeoDataFrames.union,xdf.geometry)
"""
function crop_geodf(sm::DataFrame, crop_poly::ArchGDAL.IGeometry{ArchGDAL.wkbPolygon})
    # Perform intersection
    np = [gdf.intersection(geom, crop_poly) for geom in sm.geometry]
    
    # Filter out empty geometries
    npw = filter(x -> !ArchGDAL.isempty(x), np)
    
    # Get the indices of the non-empty geometries
    indices = findall(x -> !ArchGDAL.isempty(x), np)
    
    # Create a new GeoDataFrame with the non-empty intersections
    new_gdf = DataFrame(geometry = npw)
    
    # Copy the attributes from the original GeoDataFrame
    for col in names(sm)
        if col != :geometry
            new_gdf[!, col] = sm[!, col][indices]
        end
    end
    
    return new_gdf
end

xdf = crop_geodf(sm, crop_poly)
edf = reduce(gdf.union,xdf.geometry)
plot(edf;fillcolor=:transparent,title="Union of cropped polygons")

xx = DataFrame(geometry = [crop_poly])
plot(r)
plot!(xx.geometry;fillcolor=:transparent,title="Cropped polygon")


