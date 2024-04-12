
using Plots          # generate simple plots
using Rasters
z = Raster("/mnt/d/Wasim/Goldbach/revision/v3/hgeoggb_stack.2021.nc";missingval=0)
plot(z[t=3],xlabel = "Longitude", ylabel = "Latitude")
#plot_title = string("title")

x=z[t=3]

##mask and plot...
zm = x .< 300
k = Rasters.mask(x; with=zm)
plot(k)



## Modules 
using NCDatasets     # open and manipulate NetCDFs
using Plots          # generate simple plots
using Dates          # to work with dates and time indices
# Open the file and show the metadata 
nc = NCDataset("/mnt/d/Wasim/Goldbach/revision/v3/hgeoggb_stack.2021.nc")
println(nc) 
 
# Read dimensions
lon = nc["x"][:]
lat = nc["y"][:]
time = nc["t"][:]
# depth = nc["depth"][:] NA

# Load all data of a variable (here rvar)
rvar = nc["hgeoggbstack"][:]

##check
findall(rvar[:,:,1] .> 0)[1] 
findall(rvar[:,:,1] .> 0)[1:50] 

# Load the attributes
#rvar.attrib #NA

# # Define day and depth level of interest
# date = DateTime(2022, 12, 15, 12)
# depth_level = 50

# Extract the indices 
# time_indices = findall(date .== time)[1]   
# depth_indices = findall(depth_level .<= depth)[1]

# Generate the plot 
heatmap(lon, lat, rvar[:,:, findall(rvar[:,:,1] .> 0)]', 
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string("title"))


pt="/mnt/d/Wasim/regio/out/lowres/c6/ez2/tsoilrcm_stack.2016.nc"
nc = NCDataset(pt)
println(nc) 
 
# Read dimensions
lon = nc["x"][:]
lat = nc["y"][:]
time = nc["t"][:]

#rvar = nc["tsoilrcmstack"][:]
rvar = nc["tsoilrcmstack"][:]
gr()
Plots.heatmap(lon, lat, rvar[:,:,end]', 
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string("tsoilrcmstack"))



pt="/mnt/d/Wasim/regio/out/lowres/c6/ez2/sb05rcm.2016.nc"

ds = NCDataset(pt,"r")
sbval = ds["sb05rcm"][:,:,:]
close(ds)

# Plots.heatmap(
#         sbval[1, :, :][:], 
#         sbval[:, 1, :][:], 
#         sbval[:, :, 1][:], 
#         xlabel = "Longitude", ylabel = "Latitude", 
#         plot_title = string("sb05rcm"))

rplot(pt)
nc = NCDataset(pt)
println(nc) 
 
# Read dimensions
lon = nc["x"][:]
lat = nc["y"][:]
time = nc["t"][:]

rvar = nc["sb05rcm"][:] #|> dropmissing
#msk=(rvar[:,:,1]').>0
#o=rvar[msk]
∛9

typeof(lon)
typeof(rvar[:,:,1])
Plots.plot(rvar[100:end,50:end,1])
Plots.plot(rvar[100:end,50:end,:])

Plots.plot([lon,lat,
rvar[1,:,1]
])



gr()
Plots.heatmap(lon, lat, rvar[:,:,1], 
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string("sb05rcm"))

# rvar = nc["tsoilrcmstack"][:]
# gr()
# Plots.heatmap(lon, lat, rvar[:,:,end]', 
#         xlabel = "Longitude", ylabel = "Latitude", 
#         plot_title = string("tsoilrcmstack"))

using NCDatasets
pt="/mnt/d/remo/cordex/eobs/tx_ens_mean_crop.nc"
nc = NCDataset(pt)
println(nc) 
 
# Read dimensions
lon = nc["longitude"][:]
lat = nc["latitude"][:]
time = nc["time"][:]
date = DateTime(2010, 12, 15, 0)
# Load all data of a variable (here tx)
tx = nc["tx"][:]
#tx = nc["tx"]
Dict(nc.attrib)
# Load the attributes
xx = nc["tx"]
ad = Dict(xx.attrib)
# Extract the indices 
time_indices = findall(date .== time)[1]
# Generate the plot ##array as to be []'
Plots.heatmap(lon, lat, tx[:,:,time_indices]', 
        xlabel = "Longitude", 
        ylabel = "Latitude", 
        yflip = true,
        plot_title = string(ad["long_name"],
        " on\n",
        date))


# nope.---
# xr = Raster(pt;tx=time_indices)
# #Dim{:t}=2;
# z = Raster("/mnt/d/Wasim/Goldbach/revision/v3/hgeoggb_stack.2021.nc";
#         missingval=0)
# z[Dim{:t}(Rasters.Where(x -> x >= 1 && x < 4))] |>Plots.plot
# lons, lats = map(parent, dims(z, (X, Y)))
# #z[t=2]'|>Plots.plot
# #z[t=2]|>Plots.plot
# Plots.heatmap(lons, lats, z[t=2]', 
#         xlabel = "Longitude", 
#         ylabel = "Latitude")
#         yflip = true,
#         plot_title = string(ad["long_name"],
#         " on\n",
#         date))

pt="/mnt/d/Wasim/regio/out/lowres/c8/thetrcm_stack.2003.nc"
nc = NCDataset(pt)
(nc.dim)

#getindex(g::NCDatasets.Groups,groupname::AbstractString)
getindex(nc,"y")
dimnames(nc)
dimsize(nc)
#st = parse.(string,nc.group)
nc.group["Variables"]
#nc.group|>eltype
Dict(nc.attrib)
#Dict(nc.group)
#DataFrame(nc.group)

# Load the attributes
xx = nc["t"]
#nc|>parent
lon = nc["x"][:]
lat = nc["y"][:]
println(nc)
# Extract the indices 
#time_indices = findall(date .== time)[1]
tx = nc["thetrcmstack"][:]
time_indices = 3
# Generate the plot ##array as to be []'
Plots.heatmap(lat,lon, tx[:,:,time_indices]', 
        xlabel = "Longitude", 
        ylabel = "Latitude", 
        yflip = true,
        plot_title = string(basename(pt),
        " on\n",
        date))

Plots.plot(tx[:,:,time_indices]',legend=false)

Plots.plot(tx[:,50,time_indices]',legend=false)


using NCDatasets, BenchmarkTools, Statistics
ds = NCDataset(pt)
var = "thetrcmstack"
@btime mean(ds[var])            # takes 10.180 s
@btime mean(ds[var][:,:])       # 2.369 ms, faster
close(ds)

ds = Raster(pt)
@btime mean(ds)                 #956.704 μs (3 allocations: 48 bytes) fastest
contourf(ds[t=(3)])
#Plots.plot(ds[:,:,(ds.dims)[3]|>last])

Plots.plot(ds[:,:,Int.((ds.dims)[3][9])])

Plots.plot(ds[:,:,Int.((ds.dims)[3][9])])

ds[:,:,(ds.dims)[3]|>first]
vgjl("index")

length((ds.dims)[3])
getindex((ds.dims)[3])

(ds.dims)|>eltype
ds.name
ds.metadata

zm = ds .> 0
k = Rasters.mask(ds; with=zm)
plot(k)


"/mnt/d/Wasim/regio/out/lowres/c8"|>cd
rr = readallras(".")
#map(x->Rasters.name(x),rr)
#names(rr[5])

a = string.(map(Rasters.name,rr))
b = filter(s->occursin(r"^ts",s),string.(map(Rasters.name,rr)))
typeof(indexin(b,a))
#Integer.(indexin(b,a))
ii = indexin(b,a)
rr[first(ii)]|>Plots.plot
tras = rr[first(ii)]
#tras[t=At(2)]

c=Int.((tras.dims)[3])
vgjl("Where")
tras[Dim{:t}(Rasters.Where(x -> x >= 1))] |>Plots.plot
tras[Dim{Rasters.name(tras.dims)[end]}(Rasters.Where(x -> x >= 6))] |>Plots.plot



Plots.plot(ds[:,:,Int.((ds.dims)[3][9])])



###climate
## Modules 
using NCDatasets     # open and manipulate NetCDFs
using Plots          # generate simple plots
using Dates          # to work with dates and time indices
# Open the file and show the metadata 

di="/mnt/d/remo/cordex/eobs"
fn="qq_ens_spread_0.1deg_reg_v26.0e.nc"
fn=joinpath(di,fn)

fn="/mnt/d/remo/cordex/eobs/qq_ens_spread_crop.nc"
nc = NCDataset(fn)
println(nc)
 
# Read dimensions
lon = nc["longitude"][:]
lat = nc["latitude"][:]
time = nc["time"][:]
#  rad = nc["qq"][:] #memory errors!
# Load the attributes
#rad.attrib
# Define day and depth level of interest
date = DateTime(2021, 12, 15, 00)
# Extract the indices 
time_indices = findall(date .== time)[1]
# Generate the plot 
#nc["qq"][:,:,5]

heatmap(lon, lat, nc["qq"][:,:,time_indices]', 
        xlabel = "Longitude", ylabel = "Latitude",
        plot_title = string("",nc["qq"].attrib["long_name"],
        "\non ", date))


cd(di)
pt=glob("tg")[1]
nc = NCDataset(pt)
println(nc)
 
# Read dimensions
lon = nc["longitude"][:]
lat = nc["latitude"][:]
time = nc["time"][:]
#println(nc)
plot(time, nc["tg"][5,5,:],
label=nc["tg"].attrib["units"],
title=nc["tg"].attrib["long_name"],
#yaxis=:log
)


#climateplot on time with Rasters.jl timeseries
pt="/mnt/d/Wasim/regio/out/rc200/v5/scnrun/tsoilrcm__stack.2012.nc"
fn=pt
#z=Raster(fn)
z=readras(fn)
z[X=4,Y=4]|>plot
z[X=1,Y=1]|>plot
z[X=Rasters.Center(),Y=1]|>plot

Rasters.Center

facets(z)




r=Raster(raw"C:\Users\chs72fw\Documents\EFRE_GIS\radolan\17-22\pre_pro2021.nc")
"time series"|>vgjl
"Ti"|>vgjl

r[:,:,5]'|>Plots.plot
r.dims
r[Ti=Near(DateTime(2021,6,30))] |> Plots.plot
r[Ti=(2)]|>Plots.plot

#r[Dim{:easting}=Rasters.Center(),Dim{:northing}=1,Ti=(2)]|>plot

"Dim{"|>vgjl

r[Dim{:easting}=5]

# rnew = rebuild(r; dims = (X(name(r.dims[1])), Y(name(r.dims[2])), 
#     Ti(name(r.dims[1]))))
# rnew = rebuild(r; dims = (X,Y,Ti))
# rnew.dims


using NCDatasets
pt=raw"C:\Users\chs72fw\Documents\EFRE_GIS\radolan\17-22\pre_pro2021.nc"
nc = NCDataset(pt)
#println(nc)
 
# Read dimensions
lon = nc["easting"][:]
latl = nc["northing"][:]
time = nc["time"][:]
#println(nc)
nc.attrib
plot(time, nc["pre"][2,2,:],
label=nc["pre"].attrib["units"],
title=nc["pre"].attrib["long_name"])

plot(time, r[15,10,:])

rvar = nc["pre"][:]|>skipmissing
Plots.heatmap(lon, latl, rvar[:,:,end]',
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string("title"))

rvar[:,:,4]'|>skipmissing


x="d:/remo/qm/tas-reorder.nc"
nc = NCDataset(x)
nc|>Dict        #yes

time = nc["time"][:]
nc.attrib
v="tas"

md = nc|>Dict   
(md,:auto)

#md|>collect|>x->x[3]
dict = nc|>Dict   
mykeys = keys(dict)
string.(mykeys)

v = string.(mykeys)|>lastbefore
#string.(md)
typeof(md)

DataFrame(md|>first)
DataFrame((md|>last))

plot(time, nc[v][end,end,:],
label=nc[v].attrib["units"],
title=nc[v].attrib["long_name"])


nc.dim|>collect|>first
nc.dim|>collect|>second|>last


x="pre_qdm_result.nc"
function ftsp(x::AbstractString)
        nc = NCDataset(x);
        #nc.attrib
        dict = nc|>Dict   
        mykeys = keys(dict)
        println(string.(mykeys))
        time = nc["time"][:]
        v = filter(x->!occursin(r"time|lon|lat|x|y",x),string.(mykeys))|>first
        plot(time, nc[v][end,end,:],
        label=nc[v].attrib["units"],
        title=nc[v].attrib["long_name"])        
end
"D:/remo/qm/"|>cd
ftsp("pre_qdm_result.nc")



df = hcat(time, nc[v][end,end,:])
df = DataFrame(df,:auto)
dropmissing!(df)


kk=nc[v][end,end,:]|>skipmissing|>collect
dt = time|>skipmissing
df = DataFrame(hcat(dt,kk),:auto)

typeof(time)
typeof(nc[v][end,end,:])

#indexin(nc[v][end,end,:])
#parse(df.x1,Date)
#DateTime(df.x1)
#Date(df.x1)
datetime_vector = coalesce.(time, missing)
df = hcat(nc[v][end,end,:],datetime_vector)
df = DataFrame(df,:auto)
# Parse dates in x2 column
df.x2 = Date.(string.(df.x2),"yyyy-mm-ddTHH:MM:SS")
# df.x2 = Date.(df.x2, "yyyy-mm-ddTHH:MM:SS")
# Create a line plot
plot(df.x2, df.x1, xlabel = "Date", ylabel = "Value", legend = false)


#@df df plot(df.date,cols(df.pre))

function nctodf(x::AbstractString)
        nc = NCDataset(x);
        #nc.attrib
        dict = nc|>Dict   
        mykeys = keys(dict)
        #println(string.(mykeys))
        v = filter(x->!occursin(r"time|lon|lat|x|y",x),string.(mykeys))|>first
        time = nc["time"][:]
        datetime_vector = coalesce.(time, missing)
        #df = hcat(nc[v][end,end,:],datetime_vector)
        xm = nc[v]|>size|>first
        xm = Int(round(median(1:xm);digits=0))
        ym = nc[v]|>size|>second
        ym = Int(round(median(1:ym);digits=0))
        df = DataFrame(
                v=>nc[v][xm,ym,:],      #x, y, indices
                "date"=>datetime_vector)
        # df = DataFrame(
        #         v=>nc[v][end,end,:],      #x, y, indices
        #         "date"=>datetime_vector)
        DataFrames.metadata!(df, "filename", x, style=:note);        
        #df.date = Date.(string.(df.x2),"yyyy-mm-ddTHH:MM:SS") #not needed
               # plot(time, nc[v][end,end,:],
                # label=nc[v].attrib["units"],
        # title=nc[v].attrib["long_name"])        
end

df = nctodf(x)
dy = yrsum(df)

df=nctodf("simp.nc")
bardf(df)

"nc"|>glob

df=nctodf("obsh.nc")
df=nctodf("tas_delta.nc")

x="D:/remo/cordex/wgs/utm/gb_cdx_qm/tas-obs1980-1981.nc"
x|>nctodf|>bardf

nc = NCDataset(x)
dict = nc|>Dict   
mykeys = keys(dict)
println(string.(mykeys))
v = filter(x->!occursin(r"time|lon|lat|x|y",x),string.(mykeys))|>first
time = nc["time"][:]
datetime_vector = coalesce.(time, missing)
df = DataFrame(
        v=>nc[v][50,50,:],      #x, y, indices
        "date"=>datetime_vector)
DataFrames.metadata!(df, "filename", x, style=:note);    

bardfm(df)

nc[v]|>first

r=Raster(x)
plot(r[Ti(4)])

xm = nc[v]|>size|>first
xm = Int(round(median(1:xm);digits=0))
ym = nc[v]|>size|>second
ym = Int(round(median(1:ym);digits=0))
df = DataFrame(
        v=>nc[v][xm,ym,:],      #x, y, indices
        "date"=>datetime_vector)







        "6421"|>vgctl
        "6525"|>vgctl
        "radiation_24"|>vgctl

        "[landuse_table]"|>vgctl





        usf="D:/Wasim/regio/rcm200/v2/rcm.use2"
        rx=readras(usf)
        plotlyjs()
        plot(rx)
        
        usf="C:/Users/chs72fw/Documents/EFRE_GIS/Landcover/20181016_Deutschland_LC_clip_for_Ullmann/v2/DFD_LULC_DE_2014_v1_clip_AOI_Franken_Domain.tif"
        
        usf=raw"C:\Users\chs72fw\Documents\EFRE_GIS\Landcover\20181016_Deutschland_LC_clip_for_Ullmann\lcproj.tif"
        usf="C:\\Users\\chs72fw\\Documents\\EFRE_GIS\\Landcover/20181016_Deutschland_LC_clip_for_Ullmann/v2/LULC_DE_2014_nbg_200m.tif"
        rx=Raster(usf)
        typeof(rx)
        #dirname(usf)|>cd
        plot(rx)
        
        data = rx.data
        counts = count(data)
        
        summary(rx)
        M=rx.data|>collect
        filter(x->x .> 0,M)
        filter(x->x .> 0 && x .< 4,M)
        
        rsub = filter(x->x .> 0 && x .< 4,rx.data)
        "Where"|>vgjl
        #rx[rx.data(Rasters.Where(x -> x .> 0))] 
        
        usf="D:/Wasim/regio/out/rc200/r3-cl/Soil_Temperature_Stack.nc"
        rx=readras(usf)
        plot(rx[t=3])
        
        "D:/Wasim/regio/out/rc200/r3-cl/sb05rcm_1100.mit.nc"|>readras|>plot



x = usf
using NCDatasets
nc = NCDataset(x)
dict = nc|>Dict   
mykeys = keys(dict)
println(string.(mykeys))
v = filter(x->!occursin(r"time|lon|lat|x|y|spatial_ref",x),string.(mykeys))|>first
time = nc["time"][:]
datetime_vector = coalesce.(time, missing)
df = DataFrame(
        v=>nc[v][50,50,:],      #x, y, indices
        "date"=>datetime_vector)
DataFrames.metadata!(df, "filename", x, style=:note);    

###########
x="D:/remo/cordex/wgs/utm/proj_rsds_hist+rcp85_utm.nc"
nc = NCDataset(x)
# Read dimensions
lon = nc["x"][:]
latl = nc["y"][:]
time = nc["time"][:]
#rvar = nc["rsds"][:]
#println(nc)
nc.attrib

rvar = nc["rsds"][:]
gr()
Plots.heatmap(lon, latl, rvar[:,:,end]', 
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string("rsds"))
     
        
dict = nc|>Dict   
mykeys = keys(dict)
string.(mykeys)
v = string.(mykeys)|>lastbefore
DataFrame(md|>first)
plot(time, nc[v][end,end,:],
label=nc[v].attrib["units"],
title=nc[v].attrib["long_name"])
        
adf=DataFrame(hcat(time, nc[v][end,end,:]),:auto)

ind=1234
heatmap(lon, latl, nc[v][:,:,ind]', 
        xlabel = "Longitude", ylabel = "Latitude",
        plot_title = string("",nc[v].attrib["long_name"],
        "\non ", time[ind]))


r=Raster("D:/remo/cordex/wgs/utm/proj_rsds_hist+rcp85_utm.nc")
plot(r[Ti=234])

@code_llvm globdf(r"qg")

#Using @code_lowered, you can examine the intermediate form of the 
#code after certain optimizations and transformations have been applied, helping you understand how the Julia compiler operates on your code before 
#it gets converted into lower-level representations like LLVM IR.

@code_lowered globdf(r"qg")

#This provides insights into how the code is optimized and typed by the Julia compiler.
@code_typed globdf(r"ad")
@code_typed readras(x)




str = raw"D:\Wasim\regio\out\v9\evaprcm.2017.nc"
str = "D:/remo/cordex/wgs/utm/gb_cordex/pre_gb.nc"
nc = NCDataset(str)
dict = nc|>Dict   
mykeys = keys(dict)
println(string.(mykeys))
v = filter(x->!occursin(r"time|lon|lat|x|y|spatial_ref",x),string.(mykeys))|>first
time = nc["time"][:]
#time = nc["t"][:]
#time = nc[v][:]
lon = nc["x"][:]
latl = nc["y"][:]
#time = nc["time"][:]   

varname = filter(x->!occursin(r"^t|^y|^x|^spatial_ref",x)    
            ,string.(mykeys))|>first
            
rvar = nc[varname][:]

Plots.heatmap(lon, latl, rvar[:,:,end]', 
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string(varname))
        # ,
        # yflip = true)

v = varname
Plots.heatmap(nc[v][:,:,550]',
        xlabel = "Longitude", ylabel = "Latitude", 
        #xticks = string.(lon), yticks = string.(latl),
                legend_title = string(varname))

datetime_vector = coalesce.(time, missing)
df = DataFrame(
        v=>nc[v][
                Int.(median(1:length(lon))),
                Int.(median(1:length(latl))),:],      #x, y, indices
        "date"=>datetime_vector)

DataFrames.metadata!(df, "filename", x, style=:note);    

Plots.plot(df[!,:date],df[!,1])

Plots.plot(datetime_vector,
        nc[v][
                Int.(median(1:length(lon))),
                Int.(median(1:length(latl))),:],
                        plot_title = string(varname),
                        legend=false,
                        yflip = true)

dk = yrsum(df)
Plots.plot(
        dk[!,:year],
        dk[!,2],
        #dk[!,Not(Cols(:year))]
        )


Plots.heatmap(nc[v][:,:,end]',
        xlabel = "Longitude", ylabel = "Latitude", 
                legend_title = string(varname))
