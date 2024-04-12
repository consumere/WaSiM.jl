pt = split(file)[1]
m = match(r".*[.]",basename(file))
outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"
println("loading...\n",pt,"\nand save it to ",outfile,"...\n")

using StatsPlots;
Plots.gr()
using DataFrames, CSV, Dates

function loaddf(path::AbstractString)
    ms=["-999","-9999","lin","log","LIN","LOG"]
    df = CSV.read(path,DataFrame,
    missingstring=ms,
    delim="\t",comment="-",
    silencewarnings=false,
    ntasks=4,downcast=true,	    
    normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
end

df = loaddf(path)


function reader(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    dfs = DataFrame[]
    for file in files
        file_path = joinpath(path, file)
        if isfile(file_path) && endswith(file, ext)
            df = loaddf(file_path)
	    fact=0.66
	    plotsize = (1600*fact,800*fact)
	    s = propertynames(df)[1:end-1]
	    p1 = @df df plot(:date, cols(s),size=plotsize,yaxis=:log)
            push!(dfs, p1)
        elseif isdir(file_path)
            dfs_in_subdir = loaddf(file_path, ext)
            #dfs = vcat(dfs, dfs_in_subdir)
        end
    end
    return dfs
    #return df
end

fact=0.66
plotsize = (1600*fact,800*fact)
s = propertynames(df)[1:end-1]
@df df plot(:date, cols(s),size=plotsize)


using Plots; unicodeplots()
lineplot(df[:,:date],df[:,1]) 
lineplot!(df[:,:date],df[:,2]) 


@df df plot(:date, cols(s), yaxis=:log,size=plotsize)

@df df corrplot(cols(s),size=plotsize)  





x,y = df[!,1],df[!,2]
density(x,legend = :topleft,size=plotsize)
density!(y,legend = :topleft,size=plotsize)

plot(qqplot(df,cols(s)), corrplot(df, cols(1:2)))

using Distributions 


plot( 
qqplot(df[!,1],df[!,2], qqline = :fit), 
qqplot(Cauchy,df[!,2]), 
qqnorm(df[!,2], qqline = :R) 
)

@df df qqplot(cols(s), legend = :topleft,size=plotsize)





out = pline(file)
println("saving plotly plot to",outfile,"...")
savefig(out,outfile)
println("done! ...")



	
	# function pline(path::AbstractString)
	    # ms="-9999"
	    # df = CSV.read(path,DataFrame,
	    # missingstring=ms,
	    # delim="\t",comment="-",
	    # silencewarnings=false,                                         
	    # normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
	    # df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
	        # df=df[:,Not(1:3)]
	    # nrows=size(df)[2]-1
	    # st=[]
	    # for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
	    # p = make_subplots(rows=nrows, cols=1, 
	    # shared_xaxes=true, 
	    # shared_yaxes=false,
	    # vertical_spacing=0.05,
	    # )
	    # for i in 1:nrows;
	            # add_trace!(p, 
	                        # scatter(x=df.date, y=df[:,i],
	            # name=st[i]),   row=i,     col=1);
	    # end
	# end

function ll(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    for file in files
        file_path = joinpath(path, file)
        if isfile(file_path) && endswith(file, ext)
           println(file_path)
        end
    end
end

function lg(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    v=[]
    for file in files
        file_path = joinpath(path, file)
        if isfile(file_path) && endswith(file, ext)
           println(file_path)
	   push!(v,file_path)
        end
    end
    return(v)
end

using Printf

  #  cwd = pwd() ? path : Nothing

function rec(ext::AbstractString)
    cwd = pwd() 
    osize = 0
    for (root, dirs, files) in walkdir(cwd)
     for file in files
     if isfile(file) && endswith(file, ext)
         nm=joinpath(root, file)
	 osize = stat(nm).size
	 #nm=basename(file)
	 @printf("%-30s %15.2f MB\n","$(nm):",osize/1024^2);
	 #@printf("%-20s\n","$(nm):")
	 #@printf("%-40s %15.2f MB\n","$(cwd):",osize/1024^2);
     end
    end 
end 
end 


function ldf(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    dfs = DataFrame[]
    outname = []
    for file in files
        file_path = joinpath(path, file)
        if isfile(file_path) && endswith(file, ext)
           p1 = loaddf(file_path)
	   push!(dfs, p1)
	   m = match(r".*[.]",basename(file_path))
	   nm = contains(basename(file_path),".") ? string(m.match,"png") : basename(file_path)*".png"
	   push!(outname,nm)
        end
    end
    return(dfs)
    return(outname)
end



@df md[2] plot(:date, :_1,size=plotsize,yaxis=:log)
@df md[2] plot(:date, cols(propertynames(md[2])[1:end-1]),size=plotsize,yaxis=:log) 

for i in 1:length(md);print(first(md[i]));end

for i in 1:length(md);
o=(md[i])
p=@df o plot(:date, cols(propertynames(o)[1:end-1]),
size=plotsize,yaxis=:log) 
savefig(p,outname)
;end



fact=0.66
plotsize = (1600*fact,800*fact)

function splot(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    for file in files
        file_path = joinpath(path, file)
        if isfile(file_path) && endswith(file, ext)
           df = loaddf(file_path)
	   m = match(r".*[.]",basename(file_path));
	   nm = contains(basename(file_path),".") ? string(m.match,"png") : basename(file_path)*".png" ;
	   p=@df df plot(:date, cols(propertynames(df)[1:end-1]),size=plotsize,yaxis=:log) ;
	   Plots.savefig(p,nm);
	   println("saved to $nm !");
        end
    end
end

splot(pwd(),"10")

Matrix(df[!,Not(:date)])|>cornerplot

M = df[!,Not(:date)]
corrplot(M)


myv=ct("wind")


df = loaddf(myv[end])
M = Matrix(df[!,Not(:date)])
corrplot(M)
@df df plot(:date, cols(propertynames(df)[1:end-1]),yaxis=:log,title=myv[end]) 

replace(myv, r"<(\w+)>"=>s"\1")

broadcast(length,myv) 
broadcast(contains("nc"),myv) 
myv[Not(broadcast(contains("nc"),myv))]    ##opposite of contain match

v=ct("wind")
v[Not(broadcast(contains("nc"),v))]

p=v[Not(broadcast(contains("nc"),v))][]       #vect -> string, [] needed!

M = v[Not(broadcast(contains("nc"),v))][]|>loaddf|>Matrix
corrplot(M[:,Not(end)])

df = v[Not(broadcast(contains("nc"),v))][]|>loaddf
@df df plot(:date, cols(propertynames(df)[1:end-1]),yaxis=:log) 


v=ct("soi")
v=ct("x")
s = v[Not(broadcast(contains(r"nc$"),v))][]
df = loaddf(s)
@df df plot(:date, cols(propertynames(df)[Not(end)]),yaxis=:log,title=basename(s) )

broadcast(x ->length(x), v ) 

z = v[Not(broadcast(x->occursin(r"nc$|xml|put",x),v))]
a=[]
for i in z; df=loaddf(i); push!(a,df);end    
broadcast(x ->size(x),a)


z  = outerjoin(a[2], a[3], on = :date,makeunique=true) 

cs=propertynames   
dfs = a #vect of dfs
nn = cs(dfs[3])[Not(end)]
rename(dfs[3][:,Not(end)], string.(nn) .* "-v3") 

string(nn) != string.(nn)        #!!!!
propertynames(df)|>string|> x -> replace(x,"_"=>"-") |> x -> replace(x,r"-|:"=>"")



getproperty(df,propertynames(df)[2]) .^.5   
replace(string(nn), "_"=>"")






@df df violin(string.(:_1), :tot_average, linewidth=0)
@df df boxplot!(string.(:_1), :tot_average, fillalpha=0.75, linewidth=2)
@df df dotplot!(string.(:_1), :tot_average, marker=(:black, stroke(0)))

@df df boxplot(string.(year.(df.date)), :tot_average, fillalpha=0.75, linewidth=2,legend=false) 

@df df violin(string.(year.(df.date)), :tot_average, linewidth=0,legend=false);
@df df boxplot!(string.(year.(df.date)), :tot_average, fillalpha=0.75, linewidth=2,legend=false);
@df df dotplot!(string.(year.(df.date)), :tot_average, marker=(:black, stroke(0)),legend=false)


@df df violin((month.(df.date)), :tot_average, linewidth=0.01,legend=false)


@df df violin(string.(month.(df.date)), :tot_average, linewidth=0.01,legend=false);
@df df boxplot!(string.(month.(df.date)), :tot_average, fillalpha=0.75, linewidth=0.25,legend=false);
@df df dotplot!(string.(month.(df.date)), :tot_average, fillalpha=0.75,marker=(:black,stroke(1)),legend=false)


@df df dotplot(string.(month.(df.date)), :tot_average, fillalpha=0.75,marker=(:black,stroke(1)),legend=false)

@df df dotplot(string.(month.(df.date)), :tot_average, alpha=0.5,marker=(stroke(0.0023)),legend=false,yaxis=:log)  



for i in eachindex(a)
    x = a[i]
    @df x plot!(:date, cols(propertynames(x)[Not(x)]),yaxis=:log)
end


using Rasters,Plots
pwd()
pt="/mnt/d/temp/saale/output/jan23/coarse-pest/coarse-pest/pestout/qgfiles"
cd(pt)
file="/mnt/d/temp/saale/output/jan23/coarse-pest/coarse-pest/pestout/qgfiles/qb__smf180_0700.sum.nc"
lyr=1
ts=read(Raster(file,missingval=-9999))
x = ts[t=lyr]


x=trim(x)

x=read(Raster(file,missingval=0))
contourf(x; dpi=300, size=(800, 400))





pt="/mnt/d/temp/julia_env/Soil_Temperature_Stack.nc"
ts=read(Raster(pt,missingval=0))
x = ts[t=2]
plot(x)


rp(pt,3)

using Shapefile

using Rasters, Plots, Dates, Shapefile

rl="/mnt/d/Relief_DGMs/FABDEM/tpi2.tif"
rl = "/mnt/d/Relief_DGMs/FABDEM/franken_selection/N50E011_FABDEM_V1-0.tif"
dm="/mnt/d/Relief_DGMs/FABDEM/franken_fabdem_domain.shp"



rl = "/mnt/d/Relief_DGMs/FABDEM/franken_selection/N50E011_FABDEM_V1-0.tif"
raster = Raster(rl; missingval=0)


cropped = Rasters.crop(raster, polygon)

plot([raster, cropped])

shapefile_name=dm
using Rasters, Plots, Dates, Shapefile
using Rasters.LookupArrays


shapefile_name="/mnt/d/Relief/domain_uf.shp"
shapefile_name="/mnt/d/Wasim/Goldbach/glo/ggb_catchment.shp"
shapefile_name="/mnt/d/Wasim/Goldbach/glo/ggb_catchment_4326.shp"
shp = Shapefile.Handle(shapefile_name).shapes[1]
shp = Shapefile.Handle(shapefile_name)






rl = "/mnt/d/Relief_DGMs/FABDEM/ufra30_fabdem.tif"
raster = Raster(rl)
crs(raster)



cropped = crop(raster; to=shp)
plot(cropped)
contourf(cropped)




function rcp(rl::AbstractString,shapefile_name::AbstractString)
    #rl ="/mnt/d/Wasim/Goldbach/revision/v5/ei__ggb.2021.nc"
    rl ="/mnt/d/Wasim/Goldbach/revision/v5/tsoilggb_stack.2021.nc"
    
    raster = Raster(rl;missingval=0,crs=EPSG(25832))

    shapefile_name="/mnt/d/Wasim/Goldbach/glo/ggb_catchment.shp"
    shp = Shapefile.Handle(shapefile_name)
    #reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)

    cropped = crop(raster[t=2]; to=shp)
    contourf(cropped)
end

wr = Raster("/mnt/d/Relief_DGMs/FABDEM/tpi2.tif";missingval=0)

raster[t=3]|>crs

pr = reproject(raster[t=2],EPSG(4326))

using ArchGDAL
dataset = ArchGDAL.read(shapefile_name)
layer = ArchGDAL.getlayer(dataset, 0)

ArchGDAL.getx(geom, i)

using Rasters

output_crs = ProjString(proj4string="+proj=utm +zone=18 +datum=WGS84 +units=m")
output_crs = EPSG(4326)

raster_in = Raster("/mnt/d/Relief_DGMs/FABDEM/tpi2.tif";missingval=0)
raster_out = Rasters.reproject(raster_in, output_crs)


using Plots;
A = raster_in


flags = Dict(
    #:tr => [2.0, 2.0],
    :tr => [500,500],
    :r => :near,
)

plot(warp(A, flags))




flags = Dict(


    :tr => [100,100],
    :r => :near,
)


flags = Dict(
    :s_srs => "epsg:25832",
    :t_srs => "epsg:4326",
    :tr => [100,100],
    :r => :near,
)

plot(warp(A, flags))




wr = read(Raster("/mnt/d/temp/saale/output/v3/spin/out.tif";missingval=-9999))
plot(wr)

flags = Dict(
    :tr => [100,100],

    :r => :near,
)
A=wr[t=1]
A|>plot
plot(warp(A,flags))



cd("/mnt/d/Wasim/regio/rcm/radolan")

flist = ct("2016-")


rl= "/mnt/d/Wasim/regio/rcm/radolan/projected.nc"
dm="/mnt/d/temp/saale/in_mf/fabdem/smf-ezg.shp"

rl = flist[4]

polygon = Shapefile.Handle(dm).shapes[1]
raster = Raster(rl; missingval=0,)
crs(raster)

xr = read(Raster(rl;crs=EPSG(25832)))
xr[Ti=4]|>plot
xr[Ti=150][:,:,1]|>plot

xr[1,:,1]|>plot
xr[1,:,:]|>plot

rl="/mnt/d/temp/saale/out_30m/v1/tsoilsmf_stack.2016.nc"
xr = read(Raster(rl;crs=EPSG(25832),missingval=0))
xr[t=5]|>plot

xr[t=20,cname="RdBl"]|>plot

xr[t=20, c=:thermal]|>plot



fnt=rl
cl=colormap("Purples",logscale=true,1000)
z = read(Raster(fnt,missingval=0))
plot(z[t=2];c=cl)
plot(z[t=3];c=cgrad(:thermal))


plot(z[t=20];c=:oslo)


cropped = Rasters.crop(xr[t=5], polygon)

d = [z[t=5];z[t=8]]
flags = Dict(    :tr => [100,100], :r => :near)

plot(d;c=:oslo)



	# using ArchGDAL
	# spatialref = ArchGDAL.importEPSG(25832)
	# ArchGDAL.toPROJ4(spatialref)
	# const AG = ArchGDAL
	# p_WGS_84 = AG.getfeature(x, 0) do feature
	    # AG.getgeom(feature, 0) do geom
	        # plot(geom; fa=0.1, title="UTM")
	    # end
	# end
	# layer=AG.readraster(x)
	# layer
	# p_WGS_84 = AG.getfeature(layer, 0) do feature
	    # AG.getgeom(feature, 0) do geom
	        # plot(geom; fa=0.1, title="UTM")
	    # end
	# end
	# l1 = ArchGDAL.getlayer(layer, 0)
	# ds=AG.read(x)
	# l1 = ArchGDAL.getlayer(ds,1)
	# Rasters.reproject
	# xx = Rasters.reproject(AG.importEPSG(4326),z)

xx = Rasters.reproject(EPSG(4326),z)                                                                    
#ERROR: ArgumentError: Cannot reproject from:    
	x="/mnt/d/Wasim/regio/out/utm_v6/station/qd__rcm_1700.sum.nc"
	x="tsoilrcm_stack.2017.nc"
	using Rasters,Plots
	z = read(Raster(x;crs=EPSG(25832),missingval=0))
	plot(z[t=1];c=cgrad(:thermal))
	# xx = Rasters.reproject(EPSG(4326),z[t=1])
	# xx = Rasters.reproject(EPSG(4326),z) error
	# xx = Rasters.reproject(EPSG(25832),z) #geht, aber bringt ja nix.
	xx = replace_missing(z,0) #masking 
	write("test.nc",xx) #works, but still transposed in R/Python
	

missingmask(z)|>plot
z	|>plot

z.metadata

v = z.dims[3] 
broadcast(x->trunc(Int,x),v)

e = length(z.dims[3])
z[t=collect(2:e)]|>plot

s = z[t=collect(2:e)]
plot(s;title="TS",c=cgrad(:thermal))    #title is wrong set


sx="/mnt/d/Wasim/regio/rcm/lulc_rcm.shp"
using Shapefile
sh=Shapefile.Handle(sx)
sh|>plot


Plots.backends()      
Plots.backend(:pythonplot)
sh|>plot


include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/functions.jl")
ct("png") 
methods(ldf)


title!("Trigonometric functions") #add
xz = z[t>1]
xz	|>plot
w=missingmask(z[t=4])
	

resample(xx, 200; 4326, :bilinear)
A=z[t=2]
nn="/mnt/d/Wasim/regio/rcm/radolan/2016-11-wgs.nc"
B=read(Raster(nn;missingval=-9999))
xb = B[Ti=5]

resample(A; to=xb)


import GeoDataFrames as GDF
x="/mnt/d/Wasim/Goldbach/out/v2/tou.shp"
df = GDF.read(x)
df
import GeoFormatTypes as GFT
df.geometry = reproject(df.geometry, GFT.EPSG(25832), GFT.EPSG(4326))
reproject
import ArchGDAL as ag
ag.reproject
df.geometry = ag.reproject(df.geometry, GFT.EPSG(25832), GFT.EPSG(4326))


##11.03.
"/mnt/d/Wasim/regio/out/lowres/c6/c8"|>cd
v = ct("tem")
rp(v[5],1)

fl=readdir()
fl[broadcast(x->occursin(Regex("nc"),x),fl)]

fl[broadcast(x->occursin(r"hu*nc",x),fl)]
fl[broadcast(x->occursin(r"^t?nc",x),fl)]
fl[broadcast(x->contains(x,r"^t?nc"),fl)]

rstring = "\\([^)]*\\)"
#match(Regex(rstring), "Dominion Diamond Corporation (DDC) ")
fl[broadcast(x->occursin(Regex("[nc]$"),x),fl)]
df=loaddf(v[1])
#dfs = ldf(pwd(),"te")
dfs = loadso(pwd(),"te")

#outerjoin(df_left, df_right, on=:id, order=:left)
all=[]

for i in 1:length(dfs)-1
    all[i] = innerjoin(dfs[i],dfs[i+1],on = :date,makeunique=true)
end


s = propertynames(df)[Not(end)]
combine(df, s[1] => sum, s[2])

#gd = groupby(df, month.(df.date));
##metaprogramming...
@view df[2, 2]

i=15
df = innerjoin(dfs[i],dfs[i+1],on = :date,makeunique=true,order=:left)
propertynames(df)
s = propertynames(df)[Not(end)] #masks last column == date     #[1:end-1]
@df df density(cols(s), legend = :topleft)

#like dcast
zz = stack(df, Not([ :date]))
#s = propertynames(zz)[2:end]
@df zz density(:value, legend = :topleft)

CSV.read(str,DataFrame)
str = v[1]
occursin(r"^\s*(?:#|$)", str)

matchall(r".nc",fl)

m=match(str,"t")
#replace(str, r"\+\p{U}" => x -> lowercase(x[2]))

using CSV
FilePathsBase.readable(v[1], allowcomments=true, commentmark='-')
FilePathsBase.readable(v[1])


"/mnt/d/Wasim/regio/out/c8"|>cd
cpal("t")
stackplot("tsoil")

#####16.03.2023###########
include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/functions.jl")
cd("/mnt/d/Wasim/regio/out/c8")
p=0#clear!(:p)
ll()
ras=readallras(pwd(),"t")
ras[5]|>plot
#ras[5]|>plots_heatmap
dfs=loadalldfs(pwd())

r=ras[7]
#Rasters.Key(r)
Rasters.getfield(r,1)
Rasters.NCD_DIM_MAP
dims(r)
#plot(r)
dims(r)[3][1]
lons, lats = map(parent, dims(r, (X, Y)))
#val = map(parent, dims(r, (Dim{:t})))


using Rasters, NCDatasets
fn = download("https://www.globsnow.info/swe/archive_v3.0/L3B_monthly_biascorrected_SWE/NetCDF4/201805_northern_hemisphere_monthly_biascorrected_swe_0.25grid.nc", "201805_northern_hemisphere_monthly_biascorrected_swe_0.25grid.nc")
# Get the crs with NCDatasets, and wrap it in `WellKnownText`:
wkt = WellKnownText(Dataset(fn)["crs"].attrib["spatial_ref"])
# Choose the :swe layer from the netcdf file
swe = Raster(fn; name=:swe, crs=wkt, mappedcrs=wkt)
# Use typemin(Int32) as the missing value.
swe = replace_missing(swe, typemin(Int32))
# And resample to EPSG(4326)
resampled = resample(swe, 1; crs=EPSG(4326))

contourf(resampled)


x = ct("tsoil")
fn = x[2]
ds = Dataset(fn)
wkt = EPSG(25832)
# Choose the :swe layer from the netcdf file
swe = Raster(fn; name=:tsoilrcmstack, crs=wkt, mappedcrs=wkt)
# Use typemin(Int32) as the missing value.
swe = replace_missing(swe, -9999)
swe = swe[t=5]
#swe|>plot

Rasters.transform
#ArgumentError: (2, 1) is not a valid permutation of dimensions 1:3

# And resample to EPSG(4326)
resampled = resample(swe, 1; crs=EPSG(4326))
contourf(resampled)

using ArchGDAL, GeoFormatTypes
const AG = ArchGDAL
ArchGDAL.reproject([[118, 34], [119, 35]], ProjString("+proj=longlat +datum=WGS84 +no_defs"), EPSG(25832))
ArchGDAL.reproject([[118, 34], [119, 35]], EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
m =[[118, 34], [119, 35]]
m = [lats.data,lons.data]
#m = broadcast(x->Integer(x),m)
oo = ArchGDAL.reproject(m, EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
#m = map(parent,swe)
dataset= ArchGDAL.readraster(fn)
#mb = ArchGDAL.getband(m,3)
AG.getdriver(dataset)
AG.nraster(dataset)
gt = AG.getgeotransform(dataset)
p = AG.getproj(dataset)
#AG.toPROJ4(AG.importWKT(p))
band = AG.getband(dataset, 3)

#rr = Raster(dataset[:, :, 3];Dims=[X,Y,Ti])
#rr = Raster(dataset[:, :, 3])

xm = ArchGDAL.reproject(mb, EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))

#x = Rasters.reproject(r,target=ProjString("+proj=longlat +datum=WGS84 +no_defs"))
x = Rasters.reproject(r,target=ProjString("+proj=longlat +datum=WGS84 +no_defs"))


# A=r
# classes = <=(15) => 10,
#           15..25 => 20,
#           25..35 => 30,
#           >(35) => 40
# #classified = classify(A, classes; others=0, missingval=0)
# classes = >=(4)=> 10, <(7)=> 500
# classified = classify(A ,classes;)
# plot(classified; c=:magma)

source = ArchGDAL.importEPSG(25832)
target = ArchGDAL.importEPSG(4326)
ArchGDAL.createcoordtrans(source, target) do transform
    point = ArchGDAL.fromWKT("POINT (1120351.57 741921.42)")
    println("Before: $(ArchGDAL.toWKT(point))")
    ArchGDAL.transform!(point, transform)
    println("After: $(ArchGDAL.toWKT(point))")
end



using Rasters, NCDatasets,Plots
#fn = "/mnt/d/remo/cordex/eobs/tx_ens_mean_crop.nc"
fn="/mnt/d/Wasim/regio/rcm/radolan/pre_pro.nc"
# Get the crs with NCDatasets, and wrap it in `WellKnownText`:
wkt = WellKnownText(Dataset(fn)["crs"].attrib["spatial_ref"])
# Choose the :swe layer from the netcdf file
swe = Raster(fn; name=:pre_pro, crs=wkt, mappedcrs=wkt)
#swe = RasterStack(fn; name=:pre_pro, crs=wkt, mappedcrs=wkt)
# Use typemin(Int32) as the missing value.
swe = replace_missing(swe, typemin(Int32))
# And resample to EPSG(4326)
resampled = resample(swe[Ti=2], 1; crs=EPSG(4326))
contourf(resampled)

r = Rasters.slice(swe,1)
plot(r)
r = swe[Ti=2]
plot(r[:,:,1])

r_sub_temp = swe[:, :, 1:10]
#r_out = reproject(r_sub_temp, EPSG(4326))

swe[Ti=2]|>plot
swe[:,:,2]|>plot
swe[:,1,2]|>plot

#swe.Where(Ti=20)
###
subset_i = 1:127
subset_j = 1:102
subset_k = 1:183
# create the subset view
subset = view(swe, subset_i, subset_j, subset_k)
size(subset)
#plot(subset[Ti=3])

dx = view(swe, :, :, 3)
plot(dx)

r = dx
# Rename dimensions
r = rename_dims(r, [:X, :Y, :Ti])
# Check new dimension names
dims(r)

swe[Ti(80)]|>plot

x = swe[Ti=At(DateTime(2016))]
plot(x)
ras = x
lon = lookup(ras, :northing) # if X is longitude
lat = lookup(ras, Y) # if Y is latitude

Rename(1 => :x, 3 => :y)

A=swe
A[Ti=1:3:12] |> plot

dims(A)
ncrename -O -v "easting","x" $a
ncrename -O -v "northing","y" $a
#ncrename -O -d old_dim_name,new_dim_name -v old_dim_name,new_dim_name input_file.nc output_file.nc
ncrename -O -d "easting","x" -v "easting","x" $a
ncrename -O -d "northing","y" -v "northing","y" $a
ncrename -O -d "variable","pre_pro" -v "variable","pre_pro" $a
cdx ${a%nc}|head
cdo -v divc,24 $a div.$a
cp -v $a $a.tmp
mv -v div.$a $a
xrcrds $a 524002.7 5537726 -1



v = Rasters.NCD_DIM_MAP
fn="/mnt/d/Wasim/regio/rcm/radolan/gzfiles/radolan_raw.nc"
#ras=Raster(fn)
#ras[Ti(3)]|>plot
wkt= WellKnownText("+proj=stere +lat_0=90 +lat_ts=90 +lon_0=10 +k=0.93301270189 +x_0=0 +y_0=0 +a=6370040 +b=6370040 +to_meter=1000 +no_defs")
#wkt = WellKnownText(Dataset(fn)["crs"].attrib["spatial_ref"])
# Choose the :swe layer from the netcdf file
swe = Raster(fn; name=:pre, crs=wkt, mappedcrs=wkt)
#swe[Ti(3)]|>plot
#swe = RasterStack(fn; name=:pre_pro, crs=wkt, mappedcrs=wkt)
# Use typemin(Int32) as the missing value.
swe = replace_missing(swe, typemin(Int32))
# And resample to EPSG(4326)
resampled = resample(swe[Ti=2], 1; crs=EPSG(4326))
#resampled = resample(swe, 1; crs=EPSG(4326))
contourf(resampled)


fn="/mnt/d/Wasim/regio/rcm/radolan/pre_pro.nc"
#wkt = WellKnownText(Dataset(fn)["crs"].attrib["spatial_ref"])
ras = Raster(fn; name=:pre_pro,crs=EPSG(25832), mappedcrs=EPSG(25832))
#ras = Raster(fn; name=:pre_pro)


qf=ct("out")
qf
lplot(qf[3])
lplotf(qf[3])

###############df stuff##############
#https://blog.devgenius.io/julia-dataframes-accc01afbaf8
df=readdf("qgesrcm.c8.2016")
#show(df, allcols=true)
describe(df, :mean, :std, :min, :median, :max)
getfield(df, :colindex)
#Rename columns
rename!(df,:_4 => :wlf , :_6 => :mittelsinn, :_11 => :sw,:_14 => :wue)
names(df)
#Filter the data-frame
df3 = df[df[:, :tot_average] .> 1.95, :]
last(df.date,6)

getindex(df.date,"2016-01-01")

df[df[:, :date] .> df.date[50], :]
df[df[:, :date] .> "2016-02-01", :]
##but
using DataFramesMeta, Dates
@rsubset(df, month(:date) == 4)
@subset(df, month.(:date) .== 4)
#The difference is that @rsubset works by row and @subset works on whole columns.
@subset(df, :wlf .> mean(:wlf))
@chain df begin
    @subset year.(:date) .== 2016
end

using StatsPlots
gr(size = (900, 950))

@df df scatter(:tot_average
              , :wlf
              , title="scatter plot"
              , xlabel = "tot_average"
              , ylabel = "wlf"
              , marker = (:bullet)    #(:star5)
              #,color=year.(df.date)    #:tot_average
              ,color=month.(df.date)    #:tot_average
              ,markeralpha = 0.7
              ,legend = :none
)


#interactive use -not working in wsl?
using StatsPlots, Interact
using Blink
w = Window()
body!(w, dataviewer(df))

df = loaddf("dou")
describe(df, :mean, :std, :min, :median, :max)


finaldf = df
finaldf[!,:year]=year.(finaldf[!,:date]) ;
#ve = combine(groupby(finaldf,:year), nrow => :count,  :wlf => sum)
#ve.year .= year.(ve.date)
#rename!(ve, :year => :year, Symbol(names(ve)[end]) => :sim)

getfield(finaldf, :colindex)
ve = combine(groupby(finaldf,:year), Symbol(names(finaldf)[end-3]) => :sim,
Symbol(names(finaldf)[end-2]) => :obs)
#ve.sim ./= ve[:, end-1] .* 100
#@show ve


xd = df
rename!(xd, Symbol(names(xd)[1]) => :sim,Symbol(names(xd)[2]) => :obs)
@df xd marginalhist(:sim, :obs)
out = nse(xd[!,:sim], xd[!,:obs])


println("NSE from $simfile:\nPredictions: $newname\nTarget: $(names(finaldf)[end-1])\nNSE: $out")


getfield(df, :colindex)

@df df corrplot(cols(1:4), grid = false)

marginalkde(xd[!,:sim], xd[!,:obs])
marginalkde(df[!,2], xd[!,1])

@df xd andrewsplot(:sim, cols(1:4), legend = :topleft)

xd[!,:year]=year.(xd[!,:date]) ;
@df xd andrewsplot(:year, cols(1:4), legend = :topleft)

@df xd andrewsplot(:year, cols(1:), legend = :topleft)


function aplot(df::DataFrame)
    df[!,:year]=year.(df[!,:date]) ;
    s = propertynames(df)[Not(end-1:end)]
    o = DataFrames.metadata(df)|>collect
    ti = "AndrewsPlot of "*basename(o[1][2])
    @df df andrewsplot(:year, cols(s), legend = :topleft,title=ti)
end

aplot(df)

z=ct("qg")
dx = readdf(z[1])
aplot(dx)


Matrix(dx[!,Not(:date)])|>cornerplot

using MultivariateStats
X = convert(Matrix, dx[:, 1:4])
M = fit(MDS, X'; maxoutdim=2)

plot(M, group=df.sim)




#fn="/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/proj_stereo.wkt"
fn="/mnt/d/Wasim/regio/rcm/radolan/gzfiles/full/tar11/rado_2011_raw.nc"
#readlines(fn)
#Dataset.read(fn)
# Get the crs with NCDatasets, and wrap it in `WellKnownText`:
#wkt = WellKnownText(Dataset(fn)["crs"].attrib["spatial_ref"])
wkt = WellKnownText("+proj=stere +lat_0=90 +lat_ts=90 +lon_0=10 +k=0.93301270189 +x_0=0 +y_0=0 +a=6370040 +b=6370040 +to_meter=1000 +no_defs")
# Choose the :swe layer from the netcdf file
swe = Raster(fn; name=:pre, crs=wkt, mappedcrs=wkt)
swe[Ti(3)]|>plot

# Use typemin(Int32) as the missing value.
#swe = replace_missing(swe, typemin(Int32))
# And resample to EPSG(4326)
resampled = resample(swe, 1; crs=EPSG(4326))


x = swe[Ti(3)];
EPSG(x)
v = resample(x, 1; crs=EPSG(4326))


fn="/mnt/d/Wasim/regio/out/lowres/c8/v2/thetrcm_stack.2016.nc"
z = Raster(fn; name=:thetrcmstack,missingval=-9999, crs=EPSG(25832), mappedcrs=EPSG(25832))

zs=z[t=3]

resampled = resample(z, 1; crs=EPSG(4326))



using ArchGDAL

# open input raster dataset
input = ArchGDAL.readraster(fn)

# specify output file name
output = "output.tif"

#layer
ds = ArchGDAL.getband(input,3)
# specify output SRS
srs = "+init=epsg:4326"

# reproject input to output
#x = reproject(ds, output, srs)

using ArchGDAL, GeoFormatTypes
zs

t = [[118, 34], [119, 35]]

t2 = zs.dims[1],zs.dims[2]

lons, lats = map(parent, dims(zs, (X, Y)))

ArchGDAL.reproject(
      zs.dims[1],zs.dims[2],
      ProjString("+proj=longlat +datum=WGS84 +no_defs"),
      EPSG(2025)
  )


  options = [
    "SRC_SRS=EPSG:25832",
    "DST_SRS=EPSG:4326",
    "FORMAT=GTiff",
    "RESAMPLE=BILINEAR"
]

flags = Dict(
    :tr => [2.0, 2.0],
    :r => :near,
)

#"FORMAT"=>"GTiff",
flags = Dict(
"SRC_SRS"=>"EPSG:25832",
"DST_SRS"=>"EPSG:4326",
"RESAMPLE"=>"BILINEAR")

flags = Dict(
    :tr => [2.0, 2.0],
 :r =>:bilinear)

b = (warp(zs, flags))

plot(zs)


using ArchGDAL
const AG = ArchGDAL
# create a vector of lons and lats
lons = [-80.1918]
lats = [25.7617]

# create a WKT string for the point
lons, lats = map(parent, dims(zs, (X, Y)))
#wkt = "POINT ($(lons[1]) $(lats[1]))"

wkt = []
for i in lons
    for j in lats
        pt = AG.fromWKT("POINT ($(i) $(j))")
        push!(wkt, pt)
    end
end

wkt
ot = ArchGDAL.reproject(
          wkt,      EPSG(25832),      EPSG(4326)  )


# create an abstract geometry object from the WKT string
#geom = AG.fromWKT(wkt)
#ProjString("+proj=longlat +datum=WGS84 +no_defs"),
ArchGDAL.reproject(
          geom,      EPSG(25832),      EPSG(4326)  )


x = AG.readraster(fn)
AG.getproj(x)
AG.getgeom(x)
xb = AG.getband(x,4) 
xb[:,1]

point = AG.createpoint(lons, lats)

pll = AG.point(lons, lats)


using ArchGDAL
const AG = ArchGDAL
# create a 57960-element Vector{Any} of x-y coordinates
v = wkt
# create an AG.LineString object from the vector
#linestring = AG.createlinestring(v)

# convert it to a WKT string
wkt = AG.toWKT(linestring)



using GDAL
#input = GDAL.open(fn, GDAL.GA_ReadOnly)
fn="/mnt/d/Wasim/regio/out/lowres/c8/v2/temprcm.2016.nc"
input=GDAL.open(fn)
output = "/mnt/d/Wasim/regio/out/lowres/c8/v2/temp.tif"
# specify warp options
options = [
    "SRC_SRS=EPSG:25832",
    "DST_SRS=EPSG:4326",
    "FORMAT=GTiff",
    "RESAMPLE=BILINEAR"
]

# warp input to output
GDAL.warp(input, output; options=options)


using Rasters,Plots
raster = Raster(fn)
replace!(raster, -9999 => 0)
plot(raster[t=2])



x="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/met0/wind_1970.txt"
df = DataFrame(CSV.File(x; missingstring="-9999",
                        skipto=6,
                        limit=typemax(Int),
                        comment="#",
                        stringtype=String)) #|>dropmissing
df=loaddf(x)
denseplot(df)


file = open(x,"r")
while !eof(file)
    #line = readuntil(file, '\n') # your stuff
    line = Base.read(file, '-') # your stuff
    print(line)
end
close(file)


function read_until_flag(file::IOStream, flag::String)
    line = readuntil(file, flag)
    return line[1:end-length(flag)]
end

file = open(x,"r")
while !eof(file)
    line = read_until_flag(file, "1971") 
    println(line)
end
close(file)


cd("/mnt/d/Wasim/ecad/rr")
x="rr_dom.txt"
#df = fread("rr_dom.txt",nrows=365,fill = T,skip="20000101",na.strings = '3e+33',select = 1:5)
#delim="\t",comment="-",
df = DataFrame(CSV.File(x; missingstring="3e+33",
                        skipto=3000,
                        limit=typemax(Int),
                        ntasks=4,downcast=true,	  
                        comment="#",
                        stringtype=String)) #|>dropmissing

#rr_dom.txt lines: 52111 

#parsing will try to detect the most consistent delimiter on the first 10 rows of the file
#hdr = CSV.File(x, header=3,skipto=50000,delim="\t",comment="#") |> DataFrame

#
df = CSV.File(x, header=1, skipto=51500,
ntasks=4,downcast=true,missingstring="3.0e33",
delim=" ",	  
ignorerepeated=true, #whether repeated (consecutive/sequential) delimiters should be ignored while parsing;
comment="#") |> DataFrame

propertynames(df)[1]
df.date = Date.(string.(df[!,1]),"yyyymmdd");
df
df=df[:,Not(1)]
metadata!(df, "filename", x, style=:note);
describe(df)

using PlotlyJS
dfplotjs(df,logy=false,fact=0.6)


df = CSV.File(x, header=1, #skipto=51500,
select = [1,6,9],
ntasks=4,downcast=true,missingstring="3.0e33",
delim=" ",	  
ignorerepeated=true, #whether repeated (consecutive/sequential) delimiters should be ignored while parsing;
comment="#") |> DataFrame
df.date = Date.(string.(df[!,1]),"yyyymmdd");
df=df[:,Not(1)]
metadata!(df, "filename", x, style=:note);
describe(df)
#propertynames(df)[2]=["A"]
rename!(df, Dict(propertynames(df)[1] => "A"))
rename!(df, Dict(propertynames(df)[2] => "B"))
#replace!(df, Dict(3.0e33 => 0))
#DataFrames.replace

start=1980
stop=1999
#df2 = df.year(start:stop)
#df[in(start:stop).(df.year)]
df2=df[in(start:stop).(year.(df.date)),:]
dfp(df2)
aplot(df2) #adds year col
dfp(df2[!,1:3])
dfplotjs(df2,logy=false,fact=0.6)


###streu
"/mnt/d/Wasim/streu"|>cd


###################lots of Rastersubsets ######################
#rr,nms=readallras(".")
rr[10]|>plot
#rr[10]|>name
#map(string,nms)
#map(s->occursin(r"ts",s),map(string,nms))
filter(s->occursin(r"^ts",s),map(string,nms))
map(name,rr)
filter(s->occursin(r"^ts",s),string.(map(name,rr)))

map(s->occursin(r"^ts",s),string.(map(name,rr)))
#mapreduce(x->x^2, +, [1:3;])
#xx = rr[map(s->occursin(r"^ts",s),string.(map(name,rr)))]
#typeof(rr)


a = string.(map(name,rr))
b = filter(s->occursin(r"^ts",s),string.(map(name,rr)))
typeof(indexin(b,a))
#Integer.(indexin(b,a))
ii = indexin(b,a)
rr[first(ii)]|>plot
only(ii) #first(ii)          #Return the one and only element of collection x

tras = rr[first(ii)]
tras[t=At(2)]
Not(tras[Band(1)])
keys(tras)

Rasters.metadata(tras)
#tras.Dim(:t)
tras[:,:,Not(1)]|>plot
tras[:,:,[2,3]]|>plot
using NCDatasets
xx = ct("tsoi")
ds=NCDataset(only(xx))
#v = getindex(ds::NCDataset,varname::AbstractString)
ds
ds.dim|>Dict|>first
getindex(ds,"tsoilstrstack")
#k = tras.dims[end]
k = map(Integer,tras.dims[end][2:end])
k = broadcast(Integer,tras.dims[end][2:end])

z=tras[t=first(k)+1]
#z=vcat(z,tras[t=first(k)+2])
#z=collect(z,tras[t=first(k)+2])
z|>plot



#for i in range(1,last(size(r)));println(describe(r[t=i]));end 
#for i in range(1,(size(r)[end]));println(describe(r[t=i]));end
function descr(r::Raster)
    nm=name(r);
for i in 1:last(size(r));
    printstyled("$nm t $i \n",color=:green)
    describe(r[:,:,i])
end
end

#descr(r)
#lastdim
#d=(name(r.dims)[end])
#for i in range(1,(size(r)[end]));describe(r[d(i)]);end #nope
#describe(r[(r.dims)[end]=1])

# r = readras(r"tso")
# r[:,:,[2,3]]|>plot
# for i in range(1,(size(r)[end]));
#     sm=r[:,:,i];
#     describe(sm);end









A=tras
A[X(Between(1,2))]
A = DimArray([1 2 3; 4 5 6], (X(10:10:20), Y(5:7)))
using DimensionalData
A[X(Rasters.Between(15, 25)), Y(Rasters.Between(4, 6.5))]
#k=(map(Integer,k))[2:end]
tras[:,:,k]|>plot
tras[Band=k]|>plot
#tras[..tras[Band(1)]]
#|>plot
tras[Dim{:t}(Rasters.Between(2,3))] |>plot
tras[Dim{:t}(Rasters.Between(2,end))] |>contourf

#subset by layer and geog.
tras[
Dim{:t}(Rasters.Between(2,end)),
X(Rasters.Between(570500.003, 600000)),
Y(Rasters.Between(5.58e6, 5.59e6))
] |>contourf

#A[X(Where(x -> x > 15)), Y(Where(x -> x in (19, 21)))]
tras[Dim{:t}(Rasters.Where(x -> x >= 1))] |>plot
tras[Dim{:t}(Rasters.Where(x -> x < 1))] |>plot


"/mnt/d/Wasim/streu/out/coarse/v2/"|>cd
i="hgeostr_stack.2005.nc"
r=read(Raster(i,missingval=0,mappedcrs=EPSG(25832)));
rn = r[Dim{:t}(Rasters.Where(x -> x >= 1 && x < 4))]
rn|>contourf
rn|>plot
rn = r[t=(Rasters.Where(x -> x >= 10 && x < 13))]
plot(rn)

using Shapefile
shp = "/mnt/d/Fernerkundungsdaten/glo30/eubasins.shp"
ezg = Shapefile.Handle(shp)
ezg.crs
#EPSG(raster)
ezg|>plot
r = Raster("/mnt/d/Fernerkundungsdaten/glo30/sglo30_basin.tif",missingval=-32768)
plot(r)
#plot!(r)
#poly = ezg
#close(r)
using ArchGDAL 
const AG = ArchGDAL 
r=AG.readraster("/mnt/d/Fernerkundungsdaten/glo30/sglo30_basin.tif")
typeof(r)
#xr=AG.getband(r)
#AG.indexof(r)
using GeoArrays
geoarray = GeoArray(r)
geoarray|>plot

#Rasters.GeoFormat.(r)
#Rasters.GeoInterface.

r = tras[Dim{:t}(Rasters.Where(x -> x >= 1))]
mask_trim(raster, poly) = trim(mask(raster; with=poly); pad=10)
#ezg.shapes[end]
msk=mask_trim(r, ezg.shapes[20:end])
plot(msk)

"/mnt/d/Wasim/streu/out/f1"|>cd
rr = readallras(pwd(),"rain")
map(name,rr)
contourf(rr[end])

"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/pestout"|>cd
rr = readallras(pwd(),"sum")
z=trim(rr[end])
contourf(z)


shp="/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2/long.shp"
ezg = Shapefile.Handle(shp)
ezg.crs
#EPSG(raster)
ezg|>plot
msk=mask_trim(z, ezg.shapes[1:end])
plot(msk)
mean(msk)

a=CSV.read("/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2/baselist",
DataFrame,header=3,downcast=true,normalizenames=true,delim=" ")

#comment="#",#header=3,#)
b=CSV.read("/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2/xlimlist",
DataFrame,header=3,
downcast=true,	    delim=" ",    normalizenames=true)
#comment="#",delim="\t")
#[a,b]
#merge(a,b)
map(propertynames, [a,b])
c=innerjoin(a, b, on = :_,makeunique=true)
c|>dropmissing
c[!,1]
#a[!,1]
#b[!,1]


#lnk="/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/ppold/brp01/qgkofab.p1.2016"
#lnk="/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/ppold/brp01/qges-merge"
lnk="/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/pestout/brp01/qgkofab.p1.2016"

function readqk(x::AbstractString)
    ms=["-9999","lin","log"]
    df = CSV.read(x,DataFrame,missingstring=ms,
    skipto=4, ntasks=4,
    #type = Int,
    delim="\t",comment="-", silencewarnings=true,
    normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
    df.DD = map(Int, df.DD)
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
    metadata!(df, "filename", x, style=:note);
end

z=readqk(lnk)
dfplot(z)
dfp(lnk)
#df = readdf(lnk)
df = readdf(lnk)


s = propertynames(df)[Not(end)]
@df df plot(:date,cols(s[1:end-4]),legend = :topright)

vibx(lnk)
vio(lnk)
dfp(df)

@less dfplot(df)   
@edit dfplot(df)   
Base.uncompressed_ast(methods(dfplot).ms[1]).code 

#run(`grep --color=auto -rIin --include=*jl -e sd`)
#run(`find -maxdepth 1 -name "*$file_ending" -type f -not -regex '.*.\(png\|svg\|html\|ftz\|ftz_0\|log\|list\|nc\|xml\|sh\|grd\|yrly\)$' | xargs -I% grep --color=auto -C2 -HIinw -e "$regex" %`)

"/mnt/c/Users/Public/Documents/Python_Scripts/julia"|>cd
vgg(["readras","jl"])
vgg(["allr","jl"])
#run(`ls -l`)
file="climate.jl"
xx="s"
run(`grep --color=auto -C2 -HIinw -e $xx $file`)
read(`grep --color=auto -C2 -HIinw -e $(Regex(regex)) $file`, String)

non_matching_extensions = ["png", "svg", "html", "ftz", "ftz_0", "log", "list", "nc", "xml", "sh", "grd", "yrly"]
file_ending="jl"
# && !in(basename(file) |> lowercase |> split(".")[end] #&&!non_matching_extensions
files = function filter(file -> endswith(file, file_ending),readdir())
cmd
    `grep --color=auto -C2 -HIinw -e allr climate.jl`
p = run(cmd, stdout=Pipe, stderr=Pipe)
if length(stderr) > 0
    print(stderr)
else
    print(stdout)
end

using Base: Pipe


ending="jl"
using Base: Pipe

function vgg(regex::AbstractString, ending::AbstractString)
    #find_cmd = `find -maxdepth 1 -name "*$ending" -type f -not -regex '.*.\(png\|svg\|html\|ftz\|ftz_0\|log\|list\|nc\|xml\|sh\|grd\|yrly\)$'`
    #files = run(find_cmd,stdout)
    files = filter(file -> endswith(file, file_ending),readdir())
    if stderr.readerror != nothing
        print(stderr)
        return
    end
    #files = strip(String(files))
    #if length(files) == 0
    if files == nothing
        println("..no match found...abort")
        return
    end
    #for file in split(files, '\n')
    for file in files
        cmd = `grep --color=auto -C2 -HIinw -e "$regex" $file`
        output = run(cmd, stdout)
        if length(stderr) > 0
            print(stderr)
        elseif length(output) > 0
            print(output)
        end
    end
end


vgg("ras","jl")


function redir_out(cmd::Cmd, out::Union{IO, AbstractString})
    p = run(cmd, stdout=out)
    wait(p)
    return out
end


#grep -rIHn -E 'LIN|LOG' --include='qgko*'
function vg2(regex::AbstractString, ending::AbstractString)
    cmd = `grep --color=always -C2 -rIHn -E "$regex" --include="*.$ending"`
    run(cmd)
end
cd("/mnt/c/Users/Public/Documents/Python_Scripts/julia")
vg2("readallras","jl")


# function vgg(args::Vector{String})
#     if length(args) < 2
#         println("... usage: vgg <regex> <file_ending>")
#         return 1
#     end
#     file_ending = args[2]
#     #non_matching_extensions = ["png", "svg", "html", "ftz", "ftz_0", "log", "list", "nc", "xml", "sh", "grd", "yrly"]
#     #regex = r"\b$(args[1])\b"  # Use \b to match the whole word
#     regex = args[1]
#     # Find files with the given file ending and exclude files with the non-matching extensions
#     #files = filter(file -> endswith(file, ".$file_ending") && !in(lowercase(splitext(basename(file))[2]), non_matching_extensions),
#     #               readdir())


function vgg(regex::AbstractString, ending::AbstractString)
    files = filter(file -> endswith(file, file_ending),readdir())
    # Iterate over the files and grep for the regex
    for file in files
        cmd = `grep --color=always -C2 -rIHnE "$regex" $file`
        run(cmd)
    end
end
##liber vg2<-vgg

vgg("Thread","jl")

function vgrep(regex, file_ending)
    # list files that start with "qgko" and end with file_ending
    #files = filter(file -> startswith(file, "qgko") && endswith(file, file_ending), readdir())
    files = filter(file -> endswith(file, file_ending), readdir())
    # loop over each file
    for file in files
        # open the file and read each line
        #xlines = readlines(file)
        #filter(z -> true, xlines) |> (x -> for i in 1:length(x) println("$i\t$(x[i])") end)

        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                # check if the line matches the regex
                if occursin(Regex(regex), line)
                    # print the file name, line number and line content
                    #println("$file:$(f.lineno):$line") <-nope
                    #m=match(regex, line)
                    #m=count(_ -> true, line) #das zählt die linechars
                    println("$file: $counter:\t $line")
                end
            end
        end
    end
end

vgrep("Base.Threads","jl")

file="climate.jl"
regex="Quan"
open(file) do f
    for line in eachline(f)
        #filter(z -> true, line) |> (x -> for i in 1:length(x) println("$i\t$(x[i])") end)

        #m=count(_ -> true,line)
        #count(line)
        
        println("$file: $m : $line")
        end
    end

lines = readlines("kge.jl")
for i in 1:length(lines)
  println("$i\t$(lines[i])")
end

# Mit der filter-Funktion
lines = readlines("kge.jl")
filter(line -> true, lines) |> (x -> for i in 1:length(x) println("$i\t$(x[i])") end)

# Mit der filter-Funktion
filter(line -> occursin(r"Base", line), readlines(open("kge.jl"))) |> println

#wc -l in julia:
open(file) do f
    println(count(_ -> true, eachline(f)))
end
# or count lines without opening the file
println(count(_ -> true, eachline("test.txt")))

a = ['a', 'b', 'c', 'b', 'd', 'a'];
b = ['a', 'b', 'c'];
indexin(a, b)

open(file) do f
    for line in eachline(f)
        # print the file name, line number and line content
        println("$line")
        println(line.lineno)
    end
end

# print the current line number in code
println("This is line number @__LINE__")

# open a file and get a File object
f = open(file)
end
for line in eachline(f)
    println("file.txt:$(f.lineno):$line")
end
# close the file
close(f)


"/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2"|>cd
#rplot("fab.slp") 
z=Raster("fab.slp") 
plot(z)

"fab.dhm"|>Raster|>plot

lk="fab.dhm"
r=AG.readraster(lk)
typeof(r)
#xr=AG.getband(r)
#AG.indexof(r)
using GeoArrays
geoarray = GeoArray(r)
geoarray|>plot

geoarray = GeoArray(AG.readraster(lk))
geoarray|>plot
rplot(lk,1)

function gplot(r::AbstractString)
    geoarray = GeoArrays.GeoArray(ArchGDAL.readraster(lk))
    geoarray|>plot
end

lk="fab.slp"
gplot(lk)

p()

lk="fab.art_bfid"
rplot(lk,1)

#To make a case-insensitive regex in Julia, you can use the Regex constructor with the i flag as an argument. For example:
r = Regex("hello", "i") # matches "hello", "Hello", "HELLO", etc.
#Alternatively, you can use the @r_str string macro with the i flag after the pattern. For example:
r = r"hello"i # same as above
#The i flag makes the whole regex case-insensitive. If you want to make only part of the regex case-insensitive, you can use the (?i) and (?-i) mode modifiers inside the pattern. For example:
r = r"(?i)hello(?-i) world" # matches "hello world", "Hello world", "HELLO world", etc. but not "hello World" or "Hello World"
#You can also use character classes to specify both lowercase and uppercase letters. For example:
r = r"[hH]ello [wW]orld" # same as above


as=filter(x -> occursin(r"flux"i,x),vv) 
dfs = []
for i in as;push!(dfs,readdf(i));end
k=[]
for i in dfs; oo=innerjoin(i,i+1,on="date");push!(k,oo);end


vv=rglob("qges")
as=filter(x -> occursin(r"hamon"i,x),vv) 
dfs = []
for i in as;push!(dfs,readdf(i));end
k=[]
for i in dfs; oo=innerjoin(dfs[1],i,on="date",makeunique=true);push!(k,oo);end 

odf=DataFrame[]
for i in dfs; 
    push!(odf,innerjoin(dfs[1],i,on="date",makeunique=true))
end 

out = reduce(vcat, odf) 

# rearrange columns using select
#df2 = select(df, :c, :a, :b)
#df2 = select(out, filter(x -> !occursin(r"date"i,x),names(out)) )

df2 = select(out, filter(x -> !occursin(r"date"i,x),names(out)),:date ) 
dfplotjs(df2;logy=true,fact=0.8)
try fdf(df2) catch; 
    xdf(df2)
    #fdf2(df2)
    println("hi, there was an error in fdf, so i took xdf!") end




# rearrange columns using select!
select!(df, :c, :a, :b)


s = propertynames(df)[occursin("date",names(df))]

filter(x -> occursin(r"date"i,x),names(df)) 


try
    x = parse(Int, "not a number")
catch e
    println("An error occurred: ", e)
finally
    println("This runs anyway")
end


fdf(df)

To use plotly themes in Julia, you can use the `layout` 
attribute of the `Plot` object to specify the `template` 
property. The `template` property can be set to one of the
 built-in themes, such as `"plotly"`, `"plotly_white"`, 
 `"plotly_dark"`, `"ggplot2"`, `"seaborn"`, `"simple_white"`, 
 or `"none"`. For example:

```julia
using PlotlyJS
trace = scatter(;x=1:4, y=[10, 15, 13, 17])
layout = Layout(;template="plotly_dark")
plot(trace, layout)
```

This will produce a scatter plot with the `"plotly_dark"` theme applied.

You can also create and register your own custom themes using the `Template` object and the `PlotlyJS.register_template` function. For example:

```julia
using PlotlyJS
# create a custom theme with a blue background and white text
my_theme = Template(;layout=attr(backgroundcolor="blue", font_color="white"))
# register the theme with a name
PlotlyJS.register_template("my_theme", my_theme)
# use the theme in a plot
trace = scatter(;x=1:4, y=[10, 15, 13, 17])
layout = Layout(;template="my_theme")
plot(trace, layout)
```

This will produce a scatter plot with the custom theme applied.
For more details on using plotly themes in Julia, see ¹ and ³.

Quelle: Unterhaltung mit Bing, 26/03/2023(1) Layout in Julia - Plotly. https://plotly.com/julia/reference/layout/ Zugegriffen 26/03/2023.
(2) Getting started with plotly in Julia. https://plotly.com/julia/getting-started/ Zugegriffen 26/03/2023.
(3) Theming and templates in Python - Plotly. https://plotly.com/python/templates/ Zugegriffen 26/03/2023.



r|>plot

"/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2"|>cd
"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/c1/spin13/"|>cd
fdd()
dd()

dfs = loadalldfs(pwd())
dfs[end]|>denseplot

z=map(x->DataFrames.metadata(x)|>collect,dfs)
z=map(x->basename(x[1][2]),z)
filter(x->startswith("q",x),z)

y=dfs[2]
only(DataFrames.metadata(y))[2]

basename(only(DataFrames.metadata(y))[2])

filter(x->startswith("clo",
basename(only(DataFrames.metadata(x))[2])),dfs)
#collect(basename(only(DataFrames.metadata(x))[2])),

nms = map(x->basename(only(DataFrames.metadata(x))[2]),dfs)

filter(n->occursin(r"qo",n),
map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
)

# findall("al",z[2])
# map(x->occursin(r"wol",x),z)
|>only

getf("qout")
#set other backend
plotlyjs()

homreg()
cdb()
cd("regio/control")
vgrep("50m","ctl")
"/mnt/d/Wasim/main"|>cd
vgrep("art","ctl")


"/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/"
vgrep("Unterwei","txt")

df=readmeteo("brendpegel_utf8.txt")

"/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2"|>cd
rs=Raster("fab.fzs";missingval=-9999)
gr()
contourf(rs)
Plots.plot(rs;c=cgrad(:matter))

#colorschemes: https://docs.juliaplots.org/latest/generated/colorschemes/

cs= [:default		,
      :blues		,
      :bluesreds		,
      :darkrainbow		,
      :darktest		,
      :grays		,
      :greens		,
      :heat		,
      :lightrainbow		,
      :lighttest]

Plots.plot(rs;c=cgrad(cs[2]))
Plots.plot(rs;c=cgrad(cs[4]))

plotlyjs()
files=nconly("te")
rplot(files[2])

vgrep("Nordheim","txt")


###geht nicht, zwischen Distributions zu wechseln
#ENV["R_HOME"] = "\\wsl.localhost/Ubuntu-18.04/usr/lib/R"
using Pkg
Pkg.rm("RCall")

using PyCall
Pkg.gc()
using PkgCleanup


nk="/mnt/d/Wasim/Tanalys/DEM/brend_fab/in2/fab.intern.sl.nc"
rplot(nk)



"/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/"|>cd
df=readdf(r"wind")
df = yrmean(df)
PlotlyJS.plot(PlotlyJS.contour(
    x=df.year,
    y=df.Wasserkuppe,
    z=df.Wasserkuppe))  

    PlotlyJS.plot(PlotlyJS.contour(
        x=df.year,
        y=df.year,
        z=df.Wasserkuppe))  
    
    PlotlyJS.plot(PlotlyJS.contour(
        x=df.year,
        z=df.year,
        y=df.Wasserkuppe))  


cd("/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m6")
df = yrmean(readdf(r"qges"))
df = yrsum(readdf(r"qges"))
ln = Symbol.(filter(x->!occursin("year",x),names(df)))
last(ln)

fig=PlotlyJS.plot(
    PlotlyJS.contour(
    x=df.year,
    y=df.year,
    z=df[!,last(ln)])) 
Layout(yaxis_autorange="reversed")) 

o=readras(r"vap") 
replace_missing!(o,0)
rplot(o,1)

#####plotting raster raw data!
z=o[:, :, 1].data
# x=o[:, 1, :].data
# y=o[1, :, :].data
#z = z[z .>0]
#z = z .>0
y = range(1,size(z)[1])
x = range(1,size(z)[2])
PlotlyJS.plot(
    PlotlyJS.contour(
    x=x,
    y=x,
    z=z)  )


dfyrs(readdf(r"so_wi"))
dfyrs(readdf(r"qges"))
dfpjs(readdf(r"qges"))


for i in ln
    println(i)
end


myr = readras(r"qi")
rplot(myr,1)


Layout (yaxis_autorange="reversed"))
PlotlyJS.Lax


using GeoJSON, Plots
# Read in GeoJSON file
#geojson = GeoJSON.read("/mnt/d/ClimateExplorer/pressure/pa_ogr.geojson")
vgjl("GeoJSON")
vgjl("GeoDataFrames")
using GeoDataFrames
fn="/mnt/d/ClimateExplorer/pressure/pa_ogr.geojson"
fn="/mnt/d/ClimateExplorer/pressure/stations.geojson"
gdf=GeoDataFrames.read(fn)
names(gdf)
# Extract coordinates and properties
coords = [point for point in gdf.geometry]
#c = map(parent,coords)
#parent(coords)|>DataFrame
#coords|>DataFrame
# Create DataFrame
df = DataFrame(x=[coord[1] for coord in coords], y=[coord[2] for coord in coords])
Plots.plot(gdf)

Plots.plot(gdf.geometry,legend=true) #labels=gdf.Names
#annotate!(gdf.geometry, gdf.Names,:bottomleft)
Plots.annotate!(gdf.Names)
# Add label

pt = Rasters.read(fn)
plot(pt)
Rasters.points(pt)

lab = gdf[!,:Name]
desc = gdf[!,:description]

# Create plot
scatter(
    [coord[1] for coord in coords], 
    [coord[2] for coord in coords], 
    label=[x for x in lab])
#        color=desc, label=lab)

Plots.plot(coords)
Plots.plot(gdf.geometry,label=gdf.Name)
#Plots.scatter(coords)

using ArchGDAL
# Initialize DataFrame
df = DataFrame(x=[], y=[])
# Loop over points and append coordinates to DataFrame
points = [Point(1.0, 2.0), Point(3.0, 4.0), Point(5.0, 6.0)]
for point in points
    x, y = point
    push!(df, [x, y])
end

scatter(
    [x for x in coords], 
    [coord[2] for coord in coords], 
    label=[x for x in lab])








using RDatasets
using PlotlyJS
##vioplot plotly
function violin_side_by_side()
    # requires RDatasets and DataFrames
    tips = RDatasets.dataset("reshape2", "tips")
    parts = zip(
        ("Female", "Male"),
        ("positive", "negative"),
        ("#bebada", "#8dd3c7"),
        (1.0, -0.6)
    )
    traces = GenericTrace[]
    for (sex, side, color, pointpos) in parts
        sub_tips = tips[tips[!,:Sex] .== sex, :]
        sub_traces = PlotlyJS.violin(sub_tips,
            group=:Day,
            x=:TotalBill, y0=(df) -> df[1, :Day],
            side=side, orientation="h",
            marker=attr(line=attr(width=2, color=color), symbol="line-ns"),
            line_color=color,
            hoveron="points+kde", text=(df) -> "Sample length $(size(df, 1))",
            scalemode="count", scalegroup=sex, legendgroup=sex, name=sex,
            points="all", jitter=0, pointpos=pointpos,
            span=[0],
            box_visible=true, meanline_visible=true,
            showlegend=false,
        )
        sub_traces[1][:showlegend] = true
        PlotlyJS.append!(traces, sub_traces)
    end
    # TODO: make the layout
    layout = PlotlyJS.Layout(
        hovermode="closest", violinmode="overlay",
        title="Total bill distribution<br><i>scaled by number of bills per gender",
        legend_tracegroupgap=0, violingap=0, violingroupgap=0,
        yaxis=attr(showgrid=true, categoryarray=["Thur", "Fri", "Sat", "Sun"]),
    )
    PlotlyJS.plot(traces, layout)
end

violin_side_by_side()



###########waread###########
pt="/mnt/d/ClimateExplorer/precipitation/pre_ce.txt"
#df = readdf(pt)
function waread(x::AbstractString)
    ms=["-9999","lin","log","-9999.0"]
    df = CSV.read(x,
    DataFrame,
    missingstring=ms,
    ntasks=8,
    skipto=6,
    limit=typemax(Int),
    delim="\t",
    silencewarnings=false,
    normalizenames=true,drop=(i, nm) -> i == 4) #|> dropmissing
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
    metadata!(df, "filename", x, style=:note);
end
df = waread(pt)
dy = yrsum(df)
# barp(dy)
# @less barp(dy)
vgjl("subset")

dy = filter(:year => x -> x > 1950, dy) 
ln = Symbol.(filter(x->!occursin("year",x),names(dy)))
@df dy Plots.plot(:year, cols(ln),
    legend = false, 
    seriestype=:bar)

@df dy Plots.plot(:year, cols(ln),
    legend = false, 
    seriestype=:scatter)

@df dy andrewsplot(cols(ln),
    legend = false)

homreg()
cd("../control")
vg("QQ","ctl")
vg("qq","ctl")
vg("_rr","ctl")
vg("ts-0","ctl")


"/mnt/d/temp/saale/out_30m/v1-new"|>cd

dfs,ncs = readalloutput()
filterplot("qg",dfs)
filterplot("rad",ncs)
filterplot("sb05",ncs)
filterplot("sb05",dfs)
filterplot!("ts",dfs)

filterplot("qbas",dfs)
filterplot("SCNRA",dfs)
filterplot("nenbal",dfs)
filterplot("bal",dfs)
filterplot("cantemp",dfs)
filterplot("vapo",dfs)
getnames(dfs)
tdifnc()

#names(ncs[1])

map(x->name(x),ncs)


x="/mnt/d/temp/saale/in_mf/fabdem/smf.ezg"
x="/mnt/d/temp/saale/saale_25/saale.ezg"
x="/mnt/d/temp/saale/saale.ezg"
r=Raster(x)
plot(r)

filterplot("vapo",ncs)
filterplot("wurz",ncs)

dfpjs(r"q")
"sns.heatmap"|>vgpy
"heatmap"|>vgpy
"df.corr"|>vgpy




vgjl("GeoDataFrames")
using GeoDataFrames
fn=raw"C:\Users\chs72fw\Downloads\basin.geojson"
gdf=GeoDataFrames.read(fn)
names(gdf)
# Extract coordinates and properties
coords = [point for point in gdf.geometry]
#c = map(parent,coords)
#parent(coords)|>DataFrame
#coords|>DataFrame
# Create DataFrame
df = DataFrame(x=[coord[1] for coord in coords], y=[coord[2] for coord in coords])
Plots.plot(gdf)

Plots.plot(gdf.geometry,legend=true) #labels=gdf.Names
#annotate!(gdf.geometry, gdf.Names,:bottomleft)
Plots.annotate!(gdf.Names)
# Add label


pt = Rasters.read(fn)
plot(pt)
Rasters.points(pt)

##GRDC dataset
fn="c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/GRDC/2022/main/main_daily_runoff.nc"
runoff=Raster(fn)
plot(runoff[:,:,1])
#rs = runoff[Ti > DateTime(2000,1,1)]
#rs = runoff[Ti=Near(DateTime(2000))] 
rs = runoff[Ti=Rasters.Between(DateTime(2000),DateTime(2020))] 
plot(rs,yaxis=:log)


##earth engine...

using Conda
Conda.add("earthengine-api", channel="conda-forge")
#add EarthEngine
using EarthEngine
Initialize()
ENV["PYTHON"]


vgjl("shapefile")

#using Rasters, Plots, Dates 
#using Rasters.LookupArrays
using Shapefile, Plots
x=raw"L:\04-Phil\Geo1data\prj-Efre-Daten\GIS-Daten@Server K\hydrosheds\hybas_lake_eu_lev07_v1c.shp"
shapefile_name=x
shp = Shapefile.Handle(shapefile_name).shapes[1]
plot(shp)




using CairoMakie, DataFrames

df = DataFrame(year = [2014, 2015, 2016], C160 = [2.24795e6, 2.24183e6, 2.27062e6], C200 = [3.1661e6, 3.14731e6, 3.16366e6], C240 = [1.63031e6, 1.60538e6, 1.61424e6], C320 = [1.9535e6, 1.93716e6, 1.95553e6], C624 = [0.0, 0.0, 0.0], C626 = [5.28363e6, 5.23311e6, 5.22309e6], C640 = [5.22712e6, 5.19747e6, 5.20459e6], C680 = [3.04017e6, 3.06468e6, 3.13001e6], C772 = [1.10017e6, 1.15139e6, 1.19938e6], C800 = [2.67579e6, 2.72284e6, 2.79059e6], C840 = [1.7312e6, 1.7147e6, 1.73182e6], C876 = [2.25112e6, 2.23787e6, 2.25873e6], C880 = [3.80737e6, 3.8741e6, 3.95725e6], C920 = [3.74305e6, 3.75907e6, 3.80505e6])

dfm = stack(df)
dfm.variable = string.(dfm.variable)

fig = Figure()
ax = Axis(fig[1,1])
barplot!(ax,
    dfm.year,
    dfm.value,
    color=dfm.year#,
    #stack=dfm.year
)
ax.xticklabelrotation = π / 4
ax.xticklabelalign = (:right, :center)
fig



using Shapefile

function cropper(rl::AbstractString,shapefile_name::AbstractString)
    """
    klappt noch nicht so richtig...terra ftw!
    """
    rl ="D:/Wasim/regio/out/rc200/x3/cl/tsoilrcm__stack.2005.nc"   
    raster = Raster(rl;missingval=0,crs=EPSG(25832))
    #raster = raster[t=2]
    
    shapefile_name="D:/Wasim/regio/rcm200/v9/ezp2.shp"
    shp = Shapefile.Handle(shapefile_name)
    #reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)
shp.rcm
    cropped = crop(raster[t=2]; to=shp.shapes[1:5])
    cropped = crop(raster[t=2]; to=shp.shapes[3])
    contourf(cropped)
    plot!(shp.shapes[3];c=:white,alpha=0.4)
end

using GeoDataFrames 
const GDF   = GeoDataFrames
dshp = GDF.read(shapefile_name)
plot(dshp.geometry)
dshp.rcm


#for (i,pt) in enumerate(shp.shapes)
for (i,pt) in enumerate(shp.geometry)
    name = "Basin"*string(shp[!,2][i])
    Plots.annotate!(shp.geometry[i],Plots.text(name, 8, :black, :bottom, :left))
end
Plots.plot!()

@vv ("rlibrary")
using RCall
tr = @rlibrary("terra")
#...tbc


using Pkg
installed_packages = Pkg.installed()
package_names = String[]
for (package, version) in installed_packages
    push!(package_names, package)
end
println(package_names)
##load ALL installed pkgs.
for package in package_names
    eval(Meta.parse("using $package"))
end

@code_typed debuginfo=:source +(1,1)

@code_typed debuginfo=:source waread("x")

fieldnames(df)

lk=raw"C:\Users\chs72fw\Downloads\df_Wurzburg_DEU_Clima_SIunit.csv"
df = CSV.read(lk,
    DataFrame,
    limit=typemax(Int))
    #missingstring="-9999",
    #ntasks=8,
    #skipto=6,
    #delim="\t",
    #silencewarnings=false,
    #normalizenames=true,drop=(i, nm) -> i == 4) #|> dropmissing
#,"yyyy-mm-dd"
#dt = Date.(reduce(hcat,select(df,2:5)|>Matrix|>permutedims))
#dt = Date.(select(df,2:5)|>Matrix)

dt = DateTime.(df.Column1,"yyyy-mm-dd HH:MM:SS+00:00")
#df.date .= Date.(df.Column1,"yyyy-mm-dd HH:MM:SS+00:00")
df.date .= dt
df|>names|>printstyled
dfp(df[!,Cols(r"RH|date")])
dfp(df[!,Cols(r"glob_hor_rad|date")])
dfm(df[!,Cols(r"glob_hor_rad|date")])
dfm(df[!,Cols(r"wind|date")];log=true)
#dfm(df[!,7:ncol(df)])

#polygon reproject polygon 
using GeoDataFrames
const gdf = GeoDataFrames
n="D:/Wasim/regio/rcm200/v8/ezp2.shp"
g = gdf.read(n)
plot(g.geometry, fillcolor=false)

#gl = Rasters.reproject(g.geometry,EPSG(4326))
#Array{ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}, 1}
pol = ArchGDAL.IGeometry[]
for i in 1:size(g,1)
    polygon = g.geometry[i]  # assuming the first geometry is a polygon
    src = EPSG(25832)
    dst = EPSG(4326)
    reprojected_polygon = ArchGDAL.reproject(
        polygon, src, dst) #;order=:compliant
    push!(pol,reprojected_polygon)
end
# first(pol)
# plot(pol, fillcolor=false)
#gp = g
#gp.geometry .= pol #<-breaks
gp = gdf.DataFrame(nm=g.rcm,geometry=pol)
plot(gp.geometry, fillcolor=false)

#get noch nciht

using ArchGDAL,Rasters,NCDatasets
g = gdf.read("D:/Wasim/regio/rcm200/v13/cmtv13.shp")
pol = ArchGDAL.IGeometry[]
for i in 1:size(g,1)
    polygon = g.geometry[i]  # assuming the first geometry is a polygon
    src = EPSG(25832)
    dst = EPSG(4326)
    reprojected_polygon = ArchGDAL.reproject(
        polygon, src, dst) #;order=:compliant
   coords = GeoInterface.coordinates(reprojected_polygon)
        ptc = coords[1] #muss sein
        ptc_tuples = [tuple(pt...) for pt in ptc]
        prev = [(Y,X) for (X,Y) in ptc_tuples]
        #reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
    reversed_polygon = ArchGDAL.createpolygon(prev)
    push!(pol,reversed_polygon)
end




using Rasters
import ArchGDAL
using Shapefile
using Statistics
using DataFrames
function jlzonal(shapefilepath, rasterpath;agg=mean)
    shp = Shapefile.Table(shapefilepath)
    raster = Raster(rasterpath)
    #Float64.(zonal(agg, raster; of=shp, boundary=:center))
    Int.(zonal(agg, raster; of=shp, boundary=:center))
end

shapefilepath = "D:/Wasim/regio/rcm200/v11/bk1000_cropv11.shp"
rasterpath = "D:/Wasim/regio/rcm200/v11/rcm.art-bfid"

@time julia_zonal = jlzonal(shapefilepath, rasterpath;agg=maximum)
shp = Shapefile.Table(shapefilepath)
DataFrame(id=julia_zonal,bk=shp.BODART_GR)

rasterpath = "D:/Wasim/regio/rcm200/v11/rcm.art1000"
rasterpath = "D:/Wasim/regio/rcm200/v11/rcm.slp"
r = Raster(rasterpath)
@doc agheat(rasterpath;step=50,roundto=0)
wa.agheat(rasterpath;step=50,roundto=0)
wa.agcont2(rasterpath;step=33)
agcont(rasterpath)
plot(r)

digits(123) #3-element Vector{Int64}:
@doc digits

digits(1234;base=2) #11-element Vector{Int64}:
digits(11;base=2) #binary

cd("D:/Wasim/regio/coarse-pest/v5/eobs/")
myv=ct("qg")
df=dfr(first(myv))
#Matrix(df[!,Not(:date)])|>cornerplot #breaks repl
Matrix(df[!,2:3])|>cornerplot
M = df[!,Not(:date)]
#corrplot(Matrix(M)) #breaks repl
Matrix(df[!,2:3])|>corrplot

@df df[!,2:3] corrplot(cols(1:2))

#https://www.julia-vscode.org/docs/stable/userguide/formatter/
#Ctrl+Shift+P, F1 Show Command Palette
#Shift+Alt+F Format document
#Shift+Alt+← Shrink selection
using RCall
@rlibrary("terra")
r = rast(@nco "sb")
tr = RCall.rimport("terra")
tr.summary(r)

#gr()
@time setup()
fn=raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\merit_hydro\elevtn.tif"
fn=raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\merit_hydro\rivwth.tif"
r = Raster(fn)
plot(r)
cdof(fn)
op()

pt=raw"D:\Wasim\main_basin.geojson"
using GeoDataFrames
const gdf = GeoDataFrames
g = gdf.read(pt)
g.geometry|>extrema

fn=raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\saale\saale.tif"
#r = RasterStack(fn)
#r = Raster(fn,key=:wth_1)
#r = Raster(fn,key=:saale)
r = Raster(fn)
plot(r,title="",xlab="",ylab="",legend=false)
#describe(r)

fn=raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\chirps_global.nc"
fn=raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\worldclim.nc"
r = Raster(fn)
plot(r,title="",xlab="",ylab="",legend=false)
plot(r[Ti=3],title="",xlab="",ylab="",legend=false)

@vv "Near"
r[Y(Near(45.6)),X(Near(12))]|>bar


a = Iterators.accumulate(+, [1,2,3,4]);
foreach(println, a)
collect(a)
(c, rest) = Iterators.peel(a); # peel off the first element
collect(rest)

Iterators.drop(a,1)|>collect
Iterators.take(a,1)|>collect

(b, rest) = Iterators.peel("abc"); # peel off the first element
collect(rest)


166006 files in directory
/home/ubu/.julia/conda/3:                           8.55 GB

##conda bin:
/home/ubu/.julia/conda/3/condabin/conda
/home/ubu/.julia/conda/3/condabin/conda env list
/home/ubu/.julia/conda/3/condabin/conda env --help
/home/ubu/.julia/conda/3/condabin/conda env remove -p "/home/ubu/.julia/environments/v1.8/.CondaPkg/env"
/home/ubu/.julia/conda/3/condabin/conda env list
/home/ubu/.julia/conda/3/condabin/conda doctor
/home/ubu/.julia/conda/3/condabin/conda clean -p #875.2 MB
condasize()


import WaSiM as wa
"D:/Wasim/regio/out/rc200/x22/cl5/"|>wa.towsl|>cd
wa.qba()
wa.dfp("qba")

#perl -pe "s/{/{\n'file:'\t$i\n/; s/},/}/;s/\'/\"/g"|

pnc () 
{ 
    for i in *"$1"*nc;
    do
    cx $i | perl -lne "
    s/{/{\n\"file\":\t\"$i\",\n/;
    s/.$/,/g;
    print if /{/../}/" | 
        perl -0pe "s/,\n\Z/\n}\n/;s/\'/\"/g" | 
        jq -c; done | jq -s '.' > $1.json
    echo "$1.json created"
}

jq -s '.' *.json > all.json

cx $i | perl -lne "
s/{/{\n'file:'\t'$i',\n/;
s/.$/,/g;
s/(.*),/$1}/; 
;print if /{/../}/"
#s/,(?=[^,]*$)/}/;

cx $i | perl -lne "
s/{/{\n'file:'\t'$i',\n/;
s/.$/,/g;
print if /{/../}/" | perl -pe "s/,$/}/;"

#works. last cmd back to double quotes
cx $i | perl -lne "
s/{/{\n\"file\":\t\"$i\",\n/;
s/.$/,/g;
print if /{/../}/" | perl -0pe "s/,\n\Z/\n}\n/;s/\'/\"/g" | jq 



da = """1.5 5.5 4 9 9"""
#cumsum This
da = "1.5 5.5 4 9 9"
da_array = parse.(Float64, split(da))
cumulative_sum = cumsum(da_array)


using CSV,DataFrames
df = CSV.read(IOBuffer(da), DataFrame, delim=" ", header=false, transpose=true)