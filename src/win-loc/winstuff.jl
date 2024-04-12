#winstuff
#julia
#from cmd:
#C:/Users/chs72fw/AppData/Local/Programs/Julia-1.9.0-rc1/bin/julia.exe
#C:/Users/chs72fw/AppData/Local/Programs/Julia-1.9.0-rc1/bin/julia.exe --startup-file=no --color=yes --threads auto -q
@time_imports setup()
#cd("D:/Wasim/Goldbach/revision/")
"D:/Wasim/regio/out/"|>cd
vec = rglob("so_Cap")
dfs=readall(vec)
#map(describe,dfs)
dfpall(dfs)

lal = map(yrsum,dfs)
pt="D:/Wasim/regio/out/lowres/capil/cpyr.txt"


dy = reduce((left, right) -> 
outerjoin(left, right, on = :year,makeunique=true,renamecols = uppercase => lowercase), lal)
writedf(pt,dy)
#DataFrames.metadata(df)|>only|>last|>basename

df=readdf("D:/Wasim/regio/out/lowres/capil/so_Capillar_Uprise.c2.2017")

function writewa(file::AbstractString, df::DataFrame)
    dout = df
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 24
    #df = select!(df,Symbol.(filter(x->!occursin("date",x), names(df))))
    dout = select(df, Not(:date))
    dout[!,Cols([:YY,:MM,:HH,:DD],1:end)]
    #cls = propertynames(df)|>sort|>reverse
    #df = df[!,cls[2:end]] 
    CSV.write(file, dout, transform = (col, val) -> something(val, missing),delim="\t")  
    nothing
end


df.YY = map(x ->year(x),df.date)
df.MM = map(x ->month(x),df.date)
df.DD = map(x ->day(x),df.date)
df[!, "HH"] .= 24
#select!(df, [:c, :b, :a]) # reorder columns by name
#select!(df, Not(:a)) # drop column :a
#sort!(df) # sort all columns lexicographically
#sort!(df, [order(:a, rev=true), order(:b)]) # sort column :a in descending order and column :b in ascending order
#cls = propertynames(df)|>sort|>reverse

#df = select!(df,Symbol.(filter(x->!occursin("date",x), names(df))))
select!(df, Not(:date))
df[!,Cols([:YY,:MM,:HH,:DD],1:end)]

df[!,cls[2:end]] 

df=readdf("D:/Wasim/regio/out/lowres/capil/so_Capillar_Uprise.c7.1981")
writewa("tst",df)



# sed -iE '/# time*/ d' $pt
# sed -iE '/# mode*/ d' $pt


# #all
#create_sysimage([:StatsPlots, :DataFrames, :CSV, :Dates],
using PackageCompiler
pt,nm=(pwd(),"winsys_all.so")
pts = joinpath(pt,nm)
create_sysimage(
	[:Printf, :DataFrames, :CSV, :Statistics, :Dates, :Rasters, :Distributions, :StatsPlots, :PlotlyJS],
sysimage_path=pts)



	dfs=readall(".")
	
	lal = map(yrsum,dfs)
	"D:/Wasim/Goldbach/revision/"|>cd
	vec = rglob("so_Cap")
	dfs=readall(vec)
	dfs
	vec
	vec = filter(x->!occursin("txt",x),vec)
	dfs=readall(vec)
	vec
	vec = filter(x->!occursin("old",x),vec)
	vec = filter(x->!occursin("test",x),vec)
	dfs=readall(vec)
	vec
	vec = filter(x->occursin("2021",x),vec)
	dfs=readall(vec)
	vec
	vec = filter(x->occursin("fab",x),vec)
	dfs=readall(vec)
	map(yrsum,dfs)
	dfs|>getnames



###############climRasters.
r=Raster("D:/wslconda/cmip6/out.nc")
se = r[X(Rasters.Near(200)),Y(Rasters.Near(45))]   
se|>plot  

se = r[Ti=Near(DateTime(2030))] 
se|>plot  

tm=DateTime(2050,05,05)
r[Ti=At(tm)] 

##get time axis
r[X(Rasters.Near(50)),Y(Rasters.Near(45))]|>plot

###new
cd("D:/Wasim/Tanalys/DEM/brend_fab/out/m7")
r"qges"|>dfp
r"qb"|>dfp!
plotlyjs()
using PlotlyJS

r"qb"|>readdf|>dfpjs

r"rgex"|>dfp
df = r"rad"|>readdf
describe(df)

"rad"|>facets

"rad_fab_1200.mit.nc"|>rp


cdb

my=rglob("rgex")

dfpall(my)

dfs = loadalldfs(my)
map(x->describe(x),dfs)
map(x->yrmean(x),dfs)

dfs|>last|>dfpjs

dfp(dfs[end])
bardfm!(dfs[3])

cd("D:/Wasim/Tanalys/DEM/brend_fab/out/pp/out/")


function tree(cwd::AbstractString, prefix=" ")
    paths::Vector{String} = []
	#cwd = abspath(cwd)
	cwd = relpath(cwd)
    for (looproot, dir, filenames) in walkdir(cwd)
        for relpath in dir
            push!(paths, joinpath(looproot,relpath))
        end
    end
	if length(paths)==0
		printstyled("no subfolders present in \n"*pwd(),color=:red)
		return
	end
	ap=abspath(cwd)
	printstyled("$ap\n",color=:red)
    for relpath in paths
		prefix = " "
				relpath = replace(relpath, r"[^\\]"=> ":",count=1)
				###helper to fast cdinto....
				relpath = replace(relpath, r"[^\\]*."=>"-- ",count=1)
				relpath = replace(relpath, "\\"=> "/")
		println(relpath)
	end
end

tree("../")
tree(pwd())
pwd()|>typeof


tree(raw"D:\Wasim\Tanalys\DEM")

abspath("../../../out/m7")|>cd
dfp(r"qges")
dfp(r"qbas")
dfp(r"rgex")
dfp(r"sb1")
dfp!(r"sb05")

r"sb05"|>lplot
r"sb1"|>dfl!
r"qbas"|>dfl!
r"prec"|>dfl!

mp=raw"D:\Wasim\Tanalys\DEM\Input_V2\meteo"

df = waread(joinpath(mp,"radiation_ce24.txt"))

df[:,Cols(r"WUR",:date)]|>dfp
df[:,Cols(r"BADKISSINGEN",:date)]|>dfp!

vgjl("shp")

using Pkg
Pkg.add(["Shapefile"])
using Shapefile
fn="D:/Wasim/regio/rcm200/ezg.shp"
shp = Shapefile.Handle(fn)

tr = readras("D:/Wasim/regio/out/rc200/v2/temprcm.2012.nc")
cut = mask(tr; with=shp.shapes[1:10])
p=Plots.plot
cut|>p

(shp)|>p

"D:/Wasim/regio/out/rc200/v2/precrcm.2012.nc"|>rp


# pcks=("DataFrames",  "CSV",  "Statistics",
#   "Dates",  "Rasters",  "StatsPlots",  "Distributions",
#   "Printf")
#create_sysimage([Symbol.(pcks)],sysimage_path=pts)

"D:\\Wasim\\Docker\\JuliaImagesWin\\"|>cd  
pt,nm=(pwd(),"winsys_13-05.so")
pts = joinpath(pt,nm)


#all loaded packs:
create_sysimage(sysimage_path=pts,incremental=true)


#######make image ##############
using PackageCompiler
using Plots, StatsPlots, DataFrames, CSV, Dates
pt,nm=(pwd(),"winsys_ts_plots.so")
pts = joinpath(pt,nm)
#nope... 
create_sysimage([:Plots, :StatsPlots, :DataFrames, :CSV, :Dates],
	#sysimage_path=pts,incremental=false)
	sysimage_path=pts,incremental=true)


ENV
s="D:/remo/cordex/eobs/nt.txt"
ds = CSV.read(s,DataFrame;header=false,delim=" ",
	ignoreemptyrows=true,transpose=true,stripwhitespace=true)|>dropmissing
rename!(ds,1=>Symbol.("date"))
dubs = ds[findall(x -> count(==(x), ds.date) > 1, ds.date),:]
grp = groupby(dubs,:date)
v = [size(yx,1) for yx in grp]|>unique
println("unique counts: $(v)")	

@time @time_imports setup()
pt="D:/remo/cordex/eobs/tx_ens_mean_crop.wa_ctlp"
dubplot(pt)


s="D:/remo/cordex/eobs/ts2"
ds = CSV.read(s,DataFrame;header=false,
		ignoreemptyrows=true,stripwhitespace=true)|>dropmissing
#rename!(ds,1=>Symbol.("date"))
#m = map(x->(x=>Int64),[:YY,:MM,:HH,:DD])
rename!(ds,[:YY,:MM,:DD,:HH,:hh])

join(string.(ds[!,1:3]),"-")
join.(string.(eachcol(ds)),"-")
join.((eachcol(ds)),"-")


df = ds
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
dubs = df[findall(x -> count(==(x), df.date) > 1, df.date),:]
grp = groupby(dubs,:date)
v = [size(yx,1) for yx in grp]|>unique
println("unique counts: $(v)")	




@time @time_imports setup()
cd("D:/remo/cordex/eobs/")
ps = lpro("D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt") #crds dr.

x="D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
fl = CSV.read(x,DataFrame;limit=4)
xc = fl[2,5:end]|>collect
yc = fl[3,5:end]|>collect
pts = ArchGDAL.IGeometry[]
for i in 1:length(xc)
	pt = ArchGDAL.createpoint([xc[i],yc[i]])
	pt = ArchGDAL.reproject(pt,EPSG(25832),
		ProjString("+proj=longlat +datum=WGS84 +no_defs"))
	push!(pts,pt)
end
od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)

plot(od.geometry)
od.geometry

ps.geometry[2]|>ArchGDAL.getx

df = []
for (i, pt) in enumerate(ps.geometry)
	x = ArchGDAL.getx(ps.geometry[i], 0)
	y = ArchGDAL.gety(ps.geometry[i], 0)
	nm= ps.name[i]
	#push!(df,  DataFrame([x, y, nm, [:x,:y,:name]]))
	push!(df, Dict(:x=>x,:y=>y,:name=>nm))
end
df = DataFrame(df)

Rasters.points(df.x,df.y)

# Convert observations to points
pnts = collect((ArchGDAL.getx(o,0), ArchGDAL.gety(o,0)) for o in (od.geometry) if !ismissing(o))  

st=Raster("D:/remo/cordex/eobs/qq_v28.nc")
# Extract values from raster
dou = collect(extract(st, pnts;atol=.1)|>skipmissing)
dx = DataFrame(dou)
GeoDataFrames.getgeometrycolumns(dx)
k = dx.qq[1]
plot(k, title=dx.geometry[1])
typeof(k)
#make station df
#df = DataFrame(Rasters.lookup(k,1),k.data,dx.geometry[1],[:date,:qq,:geometry])
df = DataFrame([Rasters.lookup(k,1),
	Float64.(k.data)],[:date,:qq])

wa.dubplot(df)
wa.hydromon(df;col=2)


##all in loop
cd("D:/remo/cordex/eobs/")
od = wa.lpro("D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt")
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
st = Raster("D:/remo/cordex/eobs/qq_v28.nc")
#st = Raster("D:/remo/cordex/eobs/qq_ens_spread_crop.nc") #very wrong
# Extract values from raster
dou = collect(extract(st, pnts;atol=0.1))
dx = DataFrame(dou)|>dropmissing
GeoDataFrames.getgeometrycolumns(dx)
k = dx.qq[1]
plot(k, title=dx.geometry[1])
typeof(k)
#make station df loop
df = []
for i in 1:nrow(dx)
	k = dx.qq[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,:qq]))
end
df

ds = innerjoin(df..., on= :date, makeunique=true)
#nn = map(y->join(y[1:5],y[end],"_") ,string.(od.name))
#nn = map(y->hcat(y[1:5],y[end]) ,string.(od.name))
nn = map(y->(y[1:5]*y[end]) ,string.(od.name))
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)

wa.writewa("qq_ext.wa",ds)
op()


# for i in enumerate(names(ds)[2:ncol(ds)])
# 	#print(i)
# 	rename!(ds,i[2]=>nn[i[1]])
# end
# for x in names(ds)
# 	if startswith(x,"qq")
# 	newname=replace(x,"_"=>"C", count=1)
# 	rename!(df,Dict(x=>newname))
# 	end
# end

data = string.(od.name)
# Ersetze bestimmte Muster in den Zeichenketten mithilfe von Vektoroperationen
map(x->replace(x,     r"ß" => "ss",
	"\\x" => "ue",
    r"\/" => "_",    r"_" => "-",
	r"-" => "",
    r"," => "_",    r"\xc4" => "Ae",
    r"\xd6" => "Oe",   
	r"\xdc" => "Ue",
    r"\xe4" => "ae",    r"\xf6" => "oe",
    r"\xfc" => "ue",    r"\xdf" => "ss"),data)
map(x->replace(x, 	"\xfc" => "ue",    "\xdf" => "ss",
	r",.*" => "",r"-.*" => ""),data)

#v28 test.
raw"D:\remo\cordex\eobs\v28" |>cd
#rh 

od = wa.lpro("D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt")
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
rasterpath= @nco "hu"
st = Raster(rasterpath;lazy=true) #7GB 
dou = collect(extract(st, pnts;atol=0.1)) #takes to long.
dx = DataFrame(dou)|>dropmissing
GeoDataFrames.getgeometrycolumns(dx)
k = dx.qq[1]
plot(k, title=dx.geometry[1])
typeof(k)
#make station df loop
df = []
for i in 1:nrow(dx)
	k = dx.qq[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,:qq]))
end
df

ds = innerjoin(df..., on= :date, makeunique=true)

nn = map(x->replace(x, 	"\xfc" => "ue",    "\xdf" => "ss",
r",.*" => "",r"-.*" => ""),string.(od.name))
#nn = map(y->(y[1:5]*y[end]) ,string.(od.name))
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)

wa.writewa("qq_ext.wa",ds)
op()


pt="D:/Wasim/regio/rcm200/pest/v11/gcloud/output/qkfiles/qgkorcm.p11.2016.10-06-23_13_04.txt"
cd("D:/Wasim/regio/rcm200/pest/v11/gcloud/")
cp(pt,"qg_best.txt")
df = dfr("qg_best.txt")
ob = dfr("D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt")
"D:/Wasim/regio/rcm200/pest/v11/gcloud/output/p11/route.txt"
readlines("D:/Wasim/regio/rcm200/pest/v11/gcloud/output/p11/route.txt")

k=("D:/Wasim/regio/rcm200/pest/v11/gcloud/output/p11/route.txt")
#run(`cmd.exe /c more $k`)
run(`powershell -noprofile get-content $k`)
run(`powershell -noprofile gc $k`)
run(`pwsh -NoProfile -NonInteractive -Command get-content $k`)
[split(line, '\t') for line in readlines(k)]
ENV |> grep(r"pyt"i)
ENV |> grep(r"pyt"i)|> values|>first|>U->split(U,";")

dx = mall(selt(df,:C22),selt(ob,:Wolfsmuenster))
ftp(dx)
savefig("qq-pest.png")


###remap v4 ezg
using RCall
tr=rimport("terra")
tr.terra_version()
vgr("a::classify")
#from,to
d=Dict(800=>840,680=>294,320=>240,200=>144,640=>144,920=>144,560=>400,
600=>400,720=>760)
bar(d)

pt="D:/Wasim/regio/rcm200/v4/rcm.ezg2"
agheat(pt)
cd(raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\ANALYSIS\batch_ki_dgm50")
rmdub()
du()
cd(raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman")
fdi()
du() # 148336.98 MB

cd(raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman/rcm")
rmdub()


stx=Raster("D:/remo/cordex/eobs/qq_v28.nc")


##for cor data
####
@vv "extract(st"
od = wa.lpro("D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt")
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
#st = Raster("D:/remo/cordex/eobs/v28/biascorr/pre_cor.nc";key="pre") #doesnt work in Rasters, due to dimension declaration
st = Raster("D:/remo/cordex/eobs/v28/pre/pre_cor.nc";key="pre")

# Extract values from raster
dou = collect(Rasters.extract(st, pnts;atol=.1)|>skipmissing)
dx = DataFrame(dou)
#GeoDataFrames.getgeometrycolumns(dx)
first(dx)
#k = dx[!,2][1]
k = points(dx.pre[1])
k = dx.pre[1]
first(k)
plot(k, title=dx.geometry[1])
plot(k)
typeof(k)
#make station df
k = dx.pre[1]
df = DataFrame([Rasters.lookup(k,1),
	Float64.(k.data)],[:date,:rr])

df.date .= parse.(DateTime, string.(df.date))

@df df plot(df.date,:pre)

#make station df loop
df = []
for i in 1:nrow(dx)
	k = dx.pre[i]
	tmp = DataFrame([Rasters.lookup(k,1),
	Float64.(k.data)],[:date,:pre])
	tmp.date .= parse.(DateTime, string.(tmp.date))
	push!(df, tmp)
end
df
ds = innerjoin(unique.(df, :date)..., on = :date, makeunique=true)

nn = map(x->replace(x, 	"\xfc" => "ue",    "\xdf" => "ss",
r",.*" => "",r"-.*" => ""),string.(od.name))



rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
cd("D:/remo/cordex/eobs/v28/")
wa.writewa("pre_stations.wa",ds)

first(ds)


##append headers.
crd = DataFrame(dx[!,1])
crd.name .= map(y->(y[1:2]*y[end-1]) ,nn)
rename!(crd,1=>:lon,2=>:lat)
crdf = permutedims(crd[!,[:name,:lon,:lat]])
#writedf("cords.wa",crdf)


fn="D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
# run(pipeline(`head -n 5 $fn`, 
#     stdout="station_crds.txt", 
#     stderr="errs.txt"))
#station_crds = readlines("station_crds.txt")
station_crds = readlines(fn)[1:5]
#cords_wa = readlines("cords.wa")

# Split each line into fields
station_crds = [split(line, ' ') for line in station_crds]
#cords_wa = [split(line, '\t') for line in cords_wa]

# # Replace fields in lines 3 and 4 in station_crds with fields from cords_wa
# station_crds[3][5:end] = cords_wa[3]
# station_crds[4][5:end] = cords_wa[4]
x = round.(Vector(crdf[2,:]),digits=2)
y = round.(Vector(crdf[3,:]),digits=2)
# Replace fields in lines 3 and 4 in station_crds with fields from cords_wa
station_crds[3][5:end] = string.(x)
station_crds[4][5:end] = string.(y)
# Join the fields back together
station_crds = [join(fields, '\t') for fields in station_crds]
# Write the modified contents back to "station_crds.txt"
open("tmp.txt", "w") do f
    for line in station_crds
        println(f, line)
    end
end

open("pre_stations.wa", "r") do f
    lines = readlines(f)
    open("tmp.txt", "a") do f2
        for line in lines
            println(f2, line)
        end
    end
end

npp("tmp.txt")


#######best solution.
# Read the contents of "pre_stations.wa"
pre_stations_lines = open("pre_stations.wa", "r") do f
    readlines(f)
end
# Write station_crds and all but the first line of pre_stations_lines to "pre_stations2.wa"
open("pre_stations2.wa", "w") do f
    for line in [station_crds; pre_stations_lines[2:end]]
        println(f, line)
    end
end

#now just append raw header to pre_stations.wa
station_crds = readlines(fn)[1:5]
station_crds = [split(line, ' ') for line in station_crds]
station_crds = [join(fields, '\t') for fields in station_crds]
open("pre_stations.wa", "w") do f
	for line in [station_crds; pre_stations_lines[2:end]]
		println(f, line)
	end
end
##need umlauts correction...
pwc()
cd /mnt/d/remo/cordex/eobs/v28/
#umlauts pre_stations.wa


pco = dfr("pre_stations.wa")
#pobs = dfr("D:/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/preci_1970.txt")
#pobs = dfr(fn)
###non-corrected data
fn2=raw"D:\remo\qm\corgrids\pre\pre_cor.wa"
pobs = dfr(fn2)
names(pobs)
mx = mall(selt(pco,"Wasserkuppe"),selt(pobs,"Wasserkuppe"))
qplot(mx)
@doc tline
tline(mx,:date)

selt(pobs,"Bischbrunn")|>dfm
selt(pobs,"Bischbrunn")|>dfp
dfp(pobs)
selt(pobs,"Wasserkuppe")|>baryrsum
selt(pco,"Wasserkuppe")|>baryrsum



function tline(df::DataFrame, date_col::Symbol, annotate_trendline::Bool=true)
    # Get the date column and column names for trendlines
    date_data = df[!, date_col]
    trendline_cols = setdiff(names(df), [string.(date_col)])

    p = plot()

    for y_col in trendline_cols
        y_data = df[!, y_col]

        # Perform linear regression to get the slope and intercept
        X = hcat(ones(length(date_data)), 1:length(date_data))
        y = y_data
        β = X \ y  # Linear regression

        # Extract the intercept and slope
        intercept, slope = β[1], β[2]

        # Generate the trendline using the linear equation
        trendline(x) = intercept + slope * x

        # Add the trendline to the plot
        plot!(p ,date_data, 
            trendline.(1:length(date_data)), 
            label=false,
            linewidth=2, linestyle=:dash)

        # Add annotation for the intercept and slope
        if annotate_trendline
            annotate!(mean(1:length(date_data)),mean(y),
                text("y = $(round(slope, digits=2))x + $(round(intercept, digits=2))", 8))
        end
    end

    return p
end
tline(mx,:date,annotate_trendline=true)

function tline(df::DataFrame, date_col::Symbol)
    # Get the date column and column names for trendlines
    date_data = df[!, date_col]
    trendline_cols = setdiff(names(df), [string.(date_col)])

    p = plot()

    # Initialize an empty string to store the trendline equations
    trendline_equations = ""

    for y_col in trendline_cols
        y_data = df[!, y_col]

        # Perform linear regression to get the slope and intercept
        X = hcat(ones(length(date_data)), 1:length(date_data))
        y = y_data
        β = X \ y  # Linear regression

        # Extract the intercept and slope
        intercept, slope = β[1], β[2]

        # Generate the trendline using the linear equation
        trendline(x) = intercept + slope * x
		# Add the trendline equation to the trendline_equations string
        #trendline_equations *= "$y_col: y = $(round(slope, digits=2))x + $(round(intercept, digits=2))\n"
        trendline_equations = "$y_col\ny = $(round(slope; sigdigits=2))x + $(round(intercept, sigdigits=1))"

        # plot
        plot!(p ,date_data, 
            trendline.(1:length(date_data)), 
            label=trendline_equations,
            linewidth=2, linestyle=:dash)

        
    end
    return p
end


tline(mx,:date)
ky = yrsum(pco)
tline(ky,:year)

dp = yrsum(pobs)

tline(selt(dp,4;dtcol=:year),:year)
tline(selt(dp,"Bischbrunn";dtcol=:year),:year)
tline(selt(ky,"Bischbrunn";dtcol=:year),:year)
fn
dwd = dfr(fn)|>yrsum

tline(selt(dwd,"Bischbrunn";dtcol=:year)|>dropmissing,:year)
title!("measured dwd data")

wa.tline(selt(dwd,6;dtcol=:year)|>dropmissing,:year)
selt(dwd,6;dtcol=:year)|>dfp


### from https://discourse.julialang.org/t/julia-3-times-slower-than-fortran-reading-integer-data-from-ascii-file/78516
n = 2_000_000
line = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19\n"
f = open("input.txt", "w")
write(f, "$n\n")
for i = 1:n
  write(f, line)
end
close(f)

@time open("input.txt", "r") do f
    n = parse(Int, readline(f))
    for i = 1:n
      readline(f)
    end
  end
#0.17s  

using BenchmarkTools
@btime open("input.txt", "r") do f
  n = parse(Int, readline(f))
  for i = 1:n
    readline(f)
  end
end

#152.350 ms (2002942 allocations: 122.47 MiB)
@doc CSV.File
@btime CSV.File("input.txt", delim = ' ', 
    header = false, skipto = 1) 
    
    #975.962 ms     #|> DataFrame 
@btime CSV.File("input.txt",header = 1)
#233.822 ms = 0.233822 s
@time dx=CSV.File("input.txt",header = 1)|>DataFrame;
#0.240306 seconds (2.00 M allocations: 152.811 MiB)
@time dfr(r"qges")
#0.684510 seconds 2nd: 0.005617 seconds
@edit rmlat()



###read from ncdu output
using JSON

function read_json_file(file_path::AbstractString)
    try
        # Read the content of the file
        content = read(file_path, String)
        
        # Parse the JSON data
        data = JSON.parse(content)
        
        @show data[3]
        # Process the data as needed
        #dicts = vcat(ds...)
        d = data[4]
        dfs = [DataFrame(x) for x in d]
        
        return dfs
    catch e
        println("Error reading or parsing JSON file: $e")
        return nothing  # or any other default value or behavior
    end
end

# Example usage
file_path = "D:/Wasim/regio/out/rc200/x3/cl/l.json"
ds = read_json_file(file_path)
ds[[1]]==ds[1]

dicts = vcat(ds...)
#dicts = mergewith(dicts...)
dfs = [DataFrame(x) for x in dicts]
ot = vcat(dfs[2:end]...)

#dicts = Dict.(["name", "asize", "dsize", "ino"] .=> 
#            getfield.(ds, 2))

#stack(ds)
#reduce(ds,vcat)
#reduce(ds,stack)
function convert_to_dict(file_content::AbstractString)
    try
        content = read(file_content, String)
        # Parse the JSON data
        data = JSON.parse(content)

        # Extract the relevant information and create dictionaries
        dicts = Dict.(["name", "asize", "dsize", "ino"] .=> 
            getfield.(data[end], data[2:end]))

        return dicts
    catch e
        println("Error parsing JSON data: $e")
        return nothing  # or any other default value or behavior
    end
end


content = read(file_path, String)
data = JSON.parse(content)
getfield.(data[end],3)
#di = convert_to_dict(file_path)

homes()
df = findlog(;lb=.6)
writedf(df,"qklog.txt")
npp(lat())

cd(df.nm[1]|>dirname)
@edit qbb()
hydro(r"qges")
cd("../../")
dfp(df.nm[2])
hydro(df.nm[6])

cd(df.nm[4]|>dirname)
wa.qba()
findlog(;lb=.6)
dfp("qgk")
glob("jl")
qg = dfr(r"qgk")
selt(qg,13)|>hydro
@doc routg

cd(df.nm[end]|>dirname)
qgk()

infile = raw"D:\temp\saale\control\smf180.v2.ctl"
infile = raw"D:\temp\saale\control\smf200.v7_spin1.ctl"

grep_with_context("routing_model", infile, 1)
##observed runoff! 
grep_with_context("[SpinUp", infile, 1)
grep_with_context("[groundwater_f", infile, 1)
@time facets(r"vapo")
@time facets(r"ts")

grep_with_context("routing_model", infile, 4)
npp(infile)
#saale.subset.txt

#a = readdf("D:/Wasim/Tanalys/DEM/Input_V2/meteo/saale.subset.txt")

begin 
    ofl = "route.txt"
    routeg(infile, ofl)
    sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
    sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/"
    specfile=joinpath(sfpt,sfn)
    obs = readdf(specfile)

    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
    df.name=map(x->replace(x,r"_>.*" => ""),df.name)
    sort!(df, :sim)
    sim = r"qges"|>glob|>first|>waread
    #map(x->rm(x),glob("qoutjl"))     #rm all 
    @info "
    taking prefix of sim and colnumber -> more robust merge with regex of obs!"
    outd = []
    for i in eachrow(df)
        println(i[1],"->",i[3])
        try
            dm =innerjoin(
                #sim[!,Cols(Regex(string(i[1]),"i"),end)],    
                sim[!,Cols("C"*string(i[1]),end)],
                obs[!,Cols(Regex(i[3],"i"),end)],    
                on=:date)
            onam = i[3]*"-qoutjl"
            wawrite(dm,onam)
            println("$onam saved!")
            DataFrames.metadata!(dm, "filename", onam, style=:note);
            push!(outd,dm)
            println(names(dm)," on $onam pushed!")
        catch
            onam = i[3]*"-qoutjl"
            @warn "merge is empty for $onam ! ..."
            continue
        end
    end
    vz = wa.cntcolv("outjl")
    vz = filter(xx->xx[2]==5,vz) #|>x->x[1][1]

end
npp(infile)

# kr = kge_rec()
ds = wa.kge_df3()
wa.nsx(sort(ds,2))
dfp(r"qbas")
wa.dfm(r"qbas";leg=false)
ggofjl_nosvg()
kgegrep()
kgewrite()
wa.waba()
wa.dfm(r"preci";leg=false,ann=false,fun=yrsum)
wa.dfm(r"^prec";leg=false,ann=false,fun=yrsum)
r"^prec+.*sum"|>facets
r = readras(r"^prec+.*sum")
r = r./3
@edit facets(r)


@vv "GeoStats"
##raster to Geotable +  correlation
using GeoStats
import CairoMakie as Mke

ras = Raster("D:/Wasim/regio/rcm200/v12/rcm.dep";missingval=-9999)
Rasters.replace_missing!(ras, 0)
plot(ras)
gri = CartesianGrid(size(ras))
val = ras.data

findmax(val)
findmax(ras)
Z = ras.data|>vec
table = (; Z)
typeof(table)
img = georef(table, gri)
# @doc vec
## gri2 = CartesianGrid(reverse(size(ras)))
gri2 = CartesianGrid((size(ras)))
img2 = georef((;Z=ras.data|>vec|>reverse), gri2)
img2.Z|>findmax
viewer(img2)
plot(ras)

using Random
samples = img |> Sample(10000, replace=false, rng=MersenneTwister(123))
samples |> viewer
#robust estimator:cressie
g = EmpiricalVariogram(samples, "Z", maxlag = 100.0, estimator = :cressie)
Mke.plot(g)
γ = fit(SphericalVariogram, g)
interp = samples |> InterpolateNeighbors(img.geometry, Kriging(γ))
interp |> viewer
img |> viewer
samples |> viewer

fn=raw"D:\Relief\usgs\lai\LAI-RCM200-MCD15A3H-061-results.csv"
df = CSV.read(fn,DataFrame)
println(names(df))

se = select(df, Cols(r"lai|Latitude|Longitude"i))
xmin,xmax = extrema(se.Longitude)
ymin,ymax = extrema(se.Latitude)
#lgrd = CartesianGrid(xmin,ymax,1)
#CartesianGrid((100, 100), (xmin,ymin), (1.0, 1.0))
using GeoStats
# Define the dimensions of the grid
dims = (x = (xmin, xmax), y = (ymin, ymax))
# Create the CartesianGrid
grid = CartesianGrid(dims...)
# Populate the grid with data
sdata = georef((lai = se.MCD15A3H_061_Lai_500m,), grid)
extrema(se.MCD15A3H_061_Lai_500m)
sdata = georef((;lai = se.MCD15A3H_061_Lai_500m), grid)
sdata|>viewer

density(se.MCD15A3H_061_Lai_500m)
#georef((a=rand(10), x=rand(10), y=rand(10)), (:x, :y))
tab = georef((lai=se.MCD15A3H_061_Lai_500m,
    x=se.Longitude,y=se.Latitude),(:x,:y))

tab|>viewer


#aggregate to mean by date
ag = select(df, Cols(r"MCD15A3H_061_Lai_500m|Latitude|Longitude|date"i))
map(typeof,eachcol(ag))
y = filter(x->!occursin(r"date|^L"i,x),names(ag))
ag[!, :year] = year.(ag[!,:Date])
agy = DataFrames.combine(groupby(ag, :year), y .=> mean .=> y)
dfp(agy)
rename!(agy,2=>:lai)
boxplot(agy.lai)

rename!(ag,4=>:lai,3=>:date)
#ag2 = unique(agy,:Latitude)
wa.vio(selt(ag,:lai))
wa.mbx(selt(ag,:lai)) #annotated boxplot with mean

using StatsPlots, ColorSchemes
df = selt(ag,:lai)
# Define the color gradient
colorgrad = ColorSchemes.rainbow
str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
ln = Symbol.(filter(x->!occursin(r"date|year|month"i,x),names(df)))
df.Month = month.(df.date)
mn = [ monthabbr(x) for x in (df.Month) ]

means = DataFrames.combine(groupby(df,:Month ), ln[1] => mean)
month_abbr = unique(mn)
colors = StatsPlots.colormap("Blues",12)
p = plot()  # create an empty plot
for (i, col) in enumerate(month_abbr)
    @df df boxplot!(str, col,
        fillalpha=0.75, 
        linewidth=0.25,
        seriescolor=colors[i],
        notch = true,
        whisker_width = :match,
        legend=false)
end
for i in eachrow(means)
    m = i[2]
    annotate!(i.Month - 0.5, m, #+ 1 
    text(round(m; digits=2), 6, :center, :top))
end
xticks!(0.5:11.5 , month_abbr)


samples = tab |> Sample(4000, replace=false, rng=MersenneTwister(123))
#samples |> viewer
#robust estimator:cressie
g = EmpiricalVariogram(samples, "lai", maxlag = 100.0, estimator = :cressie)
Mke.plot(g)
γ = fit(SphericalVariogram, g)
interp = samples |> InterpolateNeighbors(img.geometry, Kriging(γ))
interp |> viewer
img |> viewer
samples |> viewer

colorfunction.(means[!,2])
x = means[!,2]
x ./ maximum(x) |> Mke.plot

map(colorfunction, means[!,2])



v = means[!,2]
colors = colorfunction(v)
# Create the bar plot
bar(v, color=colors)
boxplot(v, color=colors)

wa.mbx(selt(ag,:lai)) 

function tbx(df::DataFrame)
    ti = try
        DataFrames.metadata(df)|>only|>last|>basename
    catch
    @warn "No basename in metadata!"
    ti = raw""
    end    
    str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    if (any(x->occursin("year",x),names(df)))
        df = df[!,Not("year")]
        s = Symbol.(filter(x->!occursin("date",x),names(df)))
        @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
        xticks!(0.5:11.5 , month_abbr)
        title!(ti)
    else    
    s = Symbol.(filter(x->!occursin("date",x),names(df)))
    fnt = Plots.font(
        family="sans-serif",
        pointsize=8, 
        valign=:bottom, 
        rotation= 0.0, 
        color=:black)    
    #anns = map(x->string.(df[!,ncol(df)]),monmean(df))
    #anns = monmean(df)[!,2]
    dx = monmean(df)
    anns = map(x->string.(round(x,digits=2)),dx[!,2])
    xans = map(x->Plots.text(x, fnt), anns)
    #y_offset = map(x->minimum(dx[!,x]),s) .+ map(x->mean(dx[!,x]),s)
    y_offset = dx[!,2]
    col_annotations = (Vector(dx[!, 1]) .- 0.5, y_offset, xans) #x y val            
    @df df StatsPlots.boxplot(str, cols(s), linewidth=0.1, 
        #colors=colorfunction(dx[:,2]),
        colors=:lightgrey,
        annotations=col_annotations,legend=false)

    xticks!(0.5:11.5, month_abbr)
    title!(ti)
    end
end

df = selt(ag,:lai)
tbx(df)
string.(monmean(df))

colors = colorfunction(dx[:,2])
s = Symbol.(filter(x->!occursin("date",x),names(df)))
# Plot each box separately with its own color
p = plot();  # create an empty plot
for i in 1:12
    bar(dx[:,2],
        linewidth=0.1,
        seriescolor=colors[i],
        #annotations=col_annotations,
        legend=false)
end
xticks!(0.5:11.5 , month_abbr)


bar(dx[:,2],color=colorfunction(dx[:,2]))


function tbx(df::DataFrame)
    ti = try
        DataFrames.metadata(df) |> only |> last |> basename
    catch
        @warn "No basename in metadata!"
        ti = raw""
    end

    str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
    month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    if any(x->occursin("year",x),names(df))
        df = df[!,Not("year")]
        s = Symbol.(filter(x->!occursin("date",x),names(df)))
        @df df StatsPlots.violin(str, cols(s), linewidth=0.1)
        xticks!(0.5:11.5, month_abbr)
        title!(ti)
    else
        s = Symbol.(filter(x->!occursin("date",x),names(df)))
        fnt = Plots.font(
            family="sans-serif",
            pointsize=8,
            valign=:bottom,
            rotation=0.0,
            color=:black
        )

        dx = wa.monsum(df)
        y_offset = dx[!,2]
        anns = map(x->string.(round(x,digits=2)), y_offset)
        xans = map(x->Plots.text(x, fnt), anns)
        #col_annotations = (Vector(dx[!, 1]) .- 0.5, y_offset, xans) # x y val
        #col_annotations = (Vector(dx[!, 1]) .- 0.5, y_offset, xans) # x y val
        col_annotations = (Vector(dx[!, 1]) .+ 0.5, y_offset, xans) # x y val

        # Create a colormap from the `monmean` values
        #k=wa.colorfunction(dx[!,2])
        #k = cf2(dx[!,2])
        #mat = reshape(k, (1,nrow(dx)))
        #colors = colormap.(dx[!,2])

        @df df StatsPlots.boxplot(
            str,
            cols(s),
            group = str, # group by month, needed for colormap
            linewidth=0.1,
            outliers = false,
            whisker_width = :match,
            #colors=map(x -> (x > 3 ? :green : :red), y_offset),
            #fillcolor = mat|>reverse,
            #fillcolor = Colors.gray.(1:12),
            #fillcolor = :grey,
            colors = Plots.colormap("RdBu",size(dx, 2),mid=0.2),
            annotations = col_annotations,
            legend=false,
        )

        xticks!(0.5:1.4:16, month_abbr)
        #xticks!(0.5:1.5:12.5, month_abbr)
        title!(ti)
    end
end

fn="D:/Wasim/Tanalys/DEM/Input_V2/meteo/saale_spec_ordered.09"
da = fread(fn)
dropmissing!(da)
tbx(da)
ka = dfr(r"qg")
tbx(ka)


df = selt(ag,:lai)
df = unique(df,:date)
tbx(df)

da = yrsum(da)
wa.tbx(da)

using Colors, ColorSchemes

"""
returns different colors for each month
"""
function cf2(v::Vector)
    min_val, max_val = minimum(v), maximum(v)  # calculate the min and max once, outside the loop
    range_val = max_val - min_val
    #cividis_colors = ["#00204C", "#003466", "#004680", "#005A9A", "#0071B2", "#0089CA", "#009FE3", "#2BB1EF", "#64C7F1", "#9EDBF5", "#D8EFF9", "#FFFFFF"]  # Get the colors from the cividis color scheme
    cividis_colors = [
        "#081d58",
        "#29429e",
        "#4e65cd",
        "#7989e5",
        "#a7b5ee",
        "#d3c2f8",
        "#f8e6f2",
        "#fdf0e8",
        "#f7f9e6",
        "#e0f6d5",
        "#b7f2c8",
    ]
    
    num_colors = length(cividis_colors)
    return [cividis_colors[Int(round(((x - min_val) / range_val) * (num_colors - 1)) + 1)] for x in v]
end

cf2(dx[:,2])


"""
returns different shades of grey for each month
"""
function cf2(v::Vector)
    min_val, max_val = minimum(v), maximum(v)  # calculate the min and max once, outside the loop
    range_val = max_val - min_val
    return [string("#", lpad(hex(Int(round(((x - min_val) / range_val) * 255))), 2, '0')) for x in v]
end

cf2(dx[:,2])


df_monthsum = monmean(df)
# Create box plot
rename!(df_monthsum, 2 => :mean)
boxplot(df_monthsum.mean, notch=true, annotate=true)

@df df_monthsum boxplot(group=:month,cols(:mean), notch=true, 
    median=true, label=false,notch_width=0.4, 
    whiskerprops=Dict(:linewidth=>2), notchprops=Dict(:linewidth=>2))
# annotate!(df_monthsum.month,df_monthsum.mean,df_monthsum.mean,
#     xticks=Dates.monthabbr.(df_monthsum.month) 
#     )



s = Symbol.(filter(x->!occursin("date",x),names(df)))
@df df StatsPlots.boxplot(
        cols(s),
    group = monmean(df).month)


using Plots, Dates
gr()
x   = Dates.Date(1960,1,1):Dates.Month(1):Dates.Date(2020,5,1)
plt = plot(x, randn(length(x)), alpha = 0.3, label = "");

recs = [Dates.Date(1960,04,1):Dates.Month(1):Dates.Date(1961,02,1),
        Dates.Date(1969,12,1):Dates.Month(1):Dates.Date(1970,11,1),
        Dates.Date(1973,11,1):Dates.Month(1):Dates.Date(1975,03,1),
        Dates.Date(1980,01,1):Dates.Month(1):Dates.Date(1980,07,1),
        Dates.Date(1981,07,1):Dates.Month(1):Dates.Date(1982,11,1),
        Dates.Date(1990,07,1):Dates.Month(1):Dates.Date(1991,03,1),
        Dates.Date(2001,03,1):Dates.Month(1):Dates.Date(2001,11,1),
        Dates.Date(2007,12,1):Dates.Month(1):Dates.Date(2009,06,1)]

for t in 1:8
    vspan!(plt, recs[t], label = "", color = :grey)
end

display(plt)



dates = Date(1960,04,1):Day(1):Date(2009,06,1)
ta = DataFrame(date=dates, val=rand(length(dates)))
recs = [Dates.Date(1960,04,1):Dates.Month(1):Dates.Date(1961,02,1),
        Dates.Date(1969,12,1):Dates.Month(1):Dates.Date(1970,11,1),
        Dates.Date(1973,11,1):Dates.Month(1):Dates.Date(1975,03,1),
        Dates.Date(1980,01,1):Dates.Month(1):Dates.Date(1980,07,1),
        Dates.Date(1981,07,1):Dates.Month(1):Dates.Date(1982,11,1),
        Dates.Date(1990,07,1):Dates.Month(1):Dates.Date(1991,03,1),
        Dates.Date(2001,03,1):Dates.Month(1):Dates.Date(2001,11,1),
        Dates.Date(2007,12,1):Dates.Month(1):Dates.Date(2009,06,1)]

p = plot(ta.date, ta.val, alpha = 0.3, label = "");
vspan!(p, recs, label = "", color = :grey)
for t in 1:8
    vspan!(p, recs[t], label = "", color = :grey)
end
display(p)



#D:\Wasim\regio\control\rcm200_x27-loc2.ctl
cd(raw"D:\Wasim\regio\out\rc200\x22")
kd=findlog()
sort!(kd,4,rev=true)
cdof(kd[1,end])
@edit ctl()
ctlf = wa.ctl2()
npp(ctlf)

split(ctlf,"\"")

"""
returns the position of each element in a vector
"""
function findalls(x::Vector)
    map(e -> e => findall(==(e), x), unique(x))
end
findalls(dx[:,2])
score = findalls(df.lai)
findmax(score)


function mbx(df::DataFrame)
    df.Month = month.(df.date)
    str = [ @sprintf("%02i", x) for x in (df.Month) ];
    
    month_abbr = ["Jan", "Feb", "Mar", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]

    ln = Symbol.(filter(x->!occursin(r"date|year|month"i,x),names(df)))
    # values = means[:, ncol(means)]
    # colors = colorfunction(values)

    # ln = Symbol.(filter(x->!occursin(r"date|year|month|mean"i,x),names(df)))

    p = @df df StatsPlots.boxplot(
        str , cols(ln),
        fillalpha=0.75, 
        linewidth=0.25,
        # seriescolor=:heat,
        # color = colors,
        notch = true,
        whisker_width = :match,
        legend=false)
    xticks!(0.5:11.5 , month_abbr)
    means = DataFrames.combine(groupby(df,:Month), ln[1] => mean)
    
    for i in eachrow(means)
        m = i[2]
        annotate!(i.Month - 0.5, m, #+ 1 
        text(round(m; digits=2), 6, :center, :top))
    end
    return p
end

setup()
d = raw"D:\Relief_DGMs\FABDEM\franken_selection"
rs = RasterStack(d)

@vv "erStack"
:N51E009_FABDEM_V1-0
rs|>first|>plot

r = rs|>last
r[X(Contains(9)), Y(Contains(49))]
k = grep(r"50|11",names(rs))
rs2 = rs[k]

m = mosaic(rs2)

# Combine with mosaic
#mosaic(f, regions...; missingval, atol)
#mos = mosaic(first, rs...) #breaks
plot(mos)

using GDAL
cd(d)
@doc GDAL.gdalbuildvrt("mosaic.vrt", "*.tif")

dem=Raster("fabdem.nc")
@vv "crop("
using GeoDataFrames
gd = GeoDataFrames.read(raw"D:\Wasim\main_basin.geojson")
s = "D:/Wasim/regio/rcm200/v12/catchment.shp"
gd = GeoDataFrames.read(s)
gd = wa.lpro(s)
# geom=gd.geometry|>first
# ArchGDAL.reproject(geom,
#     GeoFormatTypes.EPSG(25832),
#     GeoFormatTypes.EPSG(4326))
rc = crop(dem; to=gd.geometry|>first)
rplot(rc)
geom=gd.geometry|>first
plot!(geom,fillcolor=:transparent)
plot(dem)
plot!(geom,fillcolor=:transparent)
extrema(geom)

#pwsh
gdalbuildvrt mosaic.vrt *.tif
#gdal_translate -of GTiff -co "COMPRESS=LZW" mosaic.vrt fabdem.tif
gdal_translate -of NetCDF mosaic.vrt fabdem.nc
#gdalwarp -cutline "D:\Wasim\main_basin.geojson" -dstnodata -9999 mosaic.vrt output.nc

#on ubu18
gdalwarp -cutline "saale.shp" -dstnodata -9999 mosaic.vrt saale.nc

r=Raster("saale.nc")
plot(r)

D:\Relief_DGMs\FABDEM\franken_selection
gdalwarp -s_srs EPSG:25832 -t_srs EPSG:4326 "D:/Wasim/regio/rcm200/v12/catchment.shp" "saale.shp"
#gdalwarp -s_srs EPSG:25832 -t_srs EPSG:4326 "/mnt/d/Wasim/regio/rcm200/v12/catchment.shp" "saale.shp"

#-cutline "D:\Wasim\main_basin.geojson" -dstnodata -9999 mosaic.vrt output.nc

s = "D:/Wasim/regio/rcm200/v12/catchment.shp"
gd = wa.lpro(s)
@vv "geojson"
import GeoFormatTypes as GFT
import GeoDataFrames as GDF
outname = "saale.geojson"
geom_columns = [gd[!, :geometry]...]|>first
GDF.write(outname,gd; 
        layer_name="rcm", 
        crs=GFT.EPSG(4326), 
        driver="GeoJSON",
        geom_columns=geom_columns)
#        geom_columns=[:geometry])
using Shapefile
outname = "saale.shp"
Shapefile.write(outname,geom_columns; )
#ArchGDAL.write(outname,gd.geometry|>first) 
#Rasters.write(outname,gd.geometry|>first) 
ng = GDF.read(outname)
ng
plot(ng.geometry|>first,fillcolor=:transparent)
plot(ng.geometry,fillcolor=:transparent)
# geom=gd.geometry|>first
# ArchGDAL.reproject(geom,
#     GeoFormatTypes.EPSG(25832),
#     GeoFormatTypes.EPSG(4326))

using RCall
@rimport terra as tr
fn="D:/Wasim/regio/rcm200/v12/catchment.shp"
v = tr.vect(fn)
z = tr.project(v, "EPSG:4326")
tr.writeVector(z, "saale.shp",overwrite=true)
dem=tr.rast("mosaic.vrt")
scut = tr.crop(dem, z) #VERY FAST!
tr.res(scut)
tr.writeRaster(scut, "saale.nc",overwrite=true)
tr.writeRaster(scut, "saale.tif",names="saale",
    gdal=hcat("COMPRESS=DEFLATE", "TFW=YES"),
    overwrite=true,verbose=true)

r=Raster("saale.nc")
plot(r)
ng=GDF.read("saale.shp")
plot!(ng.geometry,fillcolor=:transparent)
missingval(r)

using Rasters, Downloads
using ArchGDAL
Yfilename = Downloads.download("https://lef-nancy.org/files/index.php/s/xGPBB8ksEiHNFsY/download")
Yraster = Raster(Yfilename)
soil_filename = Downloads.download("https://lef-nancy.org/files/index.php/s/pA8J6f5K3DxHb2T/download")
soil = Raster(soil_filename) # 1km grid
plot(soil)
missval       = missingval(soil)
missval_float = Float32(soil.missingval)
soil_cl12 = map(x-> ( (x == missval) ? missval_float : (x == 12 ? 1.0f0 : 0.0f0 ) ), soil)
soil_cl12_lr = resample(soil_cl12,to=Yraster,method=:average)
plot(soil_cl12_lr)

shp = Shapefile.Handle(fn)
cr = crop(soil_cl12_lr; to=shp.shapes[1:end])
plot(cr)

@doc wa.project
sc = crs(soil_cl12_lr)
cr = wa.project(soil; misval=missval_float, 
    #src="EPSG:6258", dst="EPSG:4326")
    src="EPSG:6258", dst="EPSG:25832")
plot(cr)
ng=GDF.read("saale.shp")
ng = GDF.read(fn)
#shp = Shapefile.Handle("saale.shp")
#crc = crop(cr; to=shp.shapes[1:end])
crc = crop(cr; to=ng.geometry)
crc = mask_trim(cr,ng.geometry)
crc = mask_trim(cr,shp.shapes[1:end])
plot(crc)
rplot(crc)
Rasters.write("soil.nc",soil)
pwc()

using GeoIO
using GeoStats
#geotable = georef
geotable = GeoIO.load("saale.shp")
#reduce(vcat,geotable) #Couldn't find dissolve
GeoIO.save("saale_catchment.geojson", geotable)

geotable = GeoIO.load("N49E011_FABDEM_V1-0.tif")
#ncg = GeoIO.load("saale.nc";flags="NetCDF")
GeoIO.formats() #NO netcdf.
ncg = GeoIO.load("saale.tif",flags="OF_VERBOSE_ERROR") #fails wit compression of tiff..
su = ncg[200:654,:]
using Makie,CairoMakie
vgjl("viz")
Triangle((0, 2), (1, 0), (1, 1))|>viz
viz(su)
su |> viewer
@time ncg |> viewer #20.355889 
#maby https://github.com/MakieOrg/GeoMakie.jl
@time ncg |> viewer #13s 2nd

@vv "bounds"
r
bnd = Rasters.bounds(r) |> collect
bnd_matrix = hcat(bnd...)
bnd_matrix = hcat(map(collect, bnd)...)
vcat(map(collect, bnd)...)|>cb

|>parent|>vec
Rasters.bounds(r)[2]

fr=raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\saale\elevation_data\merit_adapted_basins.geojson"
gd = GeoDataFrames.read(fr)
plot(gd.geometry,fillcolor=:transparent)
x = Raster(raw"C:\Users\chs72fw\.hydromt_data\artifact_data\v0.0.8\saale\elevation_data\merit_adapted\slope.tif")
plot(x) #1.15 GB
#xr = resample(x,Dict(tr=>))


# ##ubu18
# "/mnt/d/Relief_DGMs/FABDEM/wasim"
# cl="/mnt/d/Relief_DGMs/FABDEM/franken_selection/saale.shp"
# mos="/mnt/d/Relief_DGMs/FABDEM/franken_selection/mosaic.vrt"
# gdalwarp -s_srs EPSG:4326 -t_srs EPSG:25832 -r cubicspline -tr 500 500 -cutline $cl -dstnodata -9999 $mos dem2.tif
# gdalwarp -s_srs EPSG:4326 -t_srs EPSG:25832 -r cubicspline -tr 5000 5000 -cutline $cl -dstnodata -9999 $mos dem2.tif
# gdalwarp -s_srs EPSG:4326 -t_srs EPSG:25832 -r cubicspline -tr 200 200 -te 525700 621100 5546100 5603100 $mos dem2.tif
# gdalwarp -s_srs EPSG:4326 -t_srs EPSG:25832 -r cubicspline -tr 1000 1000 -te 525700 621100 5546100 5603100 $mos dem3.tif


# gdalwarp -s_srs EPSG:4326 -t_srs EPSG:25832 -r bilinear -tr 1000 1000 $mos dem2.tif
#     -te 9.38 50.07 10.69 50.58
# gdalwarp -s_srs EPSG:25832 -t_srs EPSG:4326 -te 9.38 50.07 10.69 50.58 -tr 1000 1000 ../ufra50m_fabdem.tif dem3.tif


"d:/Relief_DGMs/FABDEM/wasim"|>cd
ls()
x=readras("rcm.pur")
# p = points(x) #das sind id crds.
# first(p)
lookup(x, 1, 1)
wa.stats(x;missingval=-9999)

@doc rasterize

using GeoIO,GeoStats
#bd = GeoIO.load("rcm.asc",flags="OF_VERBOSE_ERROR")
r = ArchGDAL.readraster("rcm.asc")
lyr = 1
dx = ArchGDAL.getband(r,lyr)
x = ArchGDAL.getgeotransform(r)[2]
y = ArchGDAL.getgeotransform(r)[6]
@info "gdal cellsize info: x: $x y: $y"
ArchGDAL.getgeotransform(r)
umask = maximum(dx)
msk = 0
println("MIN:",minimum(r),"\nMAX:",maximum(r))
println("masking value range to $msk and $umask ...")
bitmat = (dx .> msk) .& (dx .<= umask)
# Filter the dx matrix using bitmat
dx_filtered = dx .* bitmat
#dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
dx_output = Matrix{Int8}(undef, size(dx, 1),size(dx,2))
#dx_output .= NaN
dx_output[bitmat] .= dx_filtered[bitmat]
dx = dx_output
xr = Raster(dx,dims=(X,Y))
plot(xr)

wa.agheat("rcm.dhk")

typeof(dx)(undef, size(dx, 1),size(dx,2))
dx_output = Matrix{Int8}(undef, size(dx, 1),size(dx,2))

add(t::Tuple{Int64, Int64}) = return (first(t) + last(t))
add(size(dx))
Rasters.cellsize(size(dx))
Rasters.cellsize(dx)
ncell(t::Tuple{Int64, Int64}) = return (first(t) * last(t))
ncell(size(dx))

wa.ncell(size(xr))
wa.ncell(xr)
lk="D:/Relief_DGMs/FABDEM/franken_selection/saale.tif"
cdof(lk)
using Rasters
using ArchGDAL
#wa.agheat(lk)
#wa.ncell(xr)
dx = Raster(lk)
dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
#dx_output .= NaN
agheat()
setup()

mw = Raster("D:/Relief_DGMs/FABDEM/wasim/merit/Saale_RiverWidth.tif",missingval=NaN)
mw = replace_missing(mw,missing)
plot(mw)
mean(mw|>skipmissing)
wa.stats(mw,missingval=missing)

#using Shapefile
fn="D:/Wasim/regio/rcm200/ezg.shp"
#shp = Shapefile.Handle(fn)
shp=wa.lpro(fn)
plot!(shp.geometry|>first,fillcolor=:transparent)

#mc = crop(mw; to=shp.geometry|>first)
mc = mask_trim(mw,shp.geometry)
describe(mw)
stats(mc)
stats(mw)

od=wa.stats("D:/Relief_DGMs/FABDEM/wasim/merit/Saale_RiverWidth.tif")
cd("D:/Relief_DGMs/FABDEM/wasim/merit/")
glob("S")
wa.stats("Saale_UpstreamArea.tif")
#mw = replace_missing(mw,missing)


#k = Raster("Saale_Elevation.tif",missingval=NaN)
k = Raster("Saale_Elevation.tif")
plot(k)
el = project(k;src="EPSG:4326",dst="EPSG:25832")
plot(el)
#replace!(el,NaN=>missing)
el = replace_missing(el,missing)
write("dem_utm.nc",el;force=true)
vgr("rcm.art1000")
vgctl("rcm.art1000")


r = Raster("D:/Relief_DGMs/FABDEM/wasim/merit/Saale_Elevation.tif")
bds = Rasters.bounds(r)
nr = Raster("D:/Relief_DGMs/FABDEM/wasim/merit50/saale50_Elevation.tif")

cr = crop(nr; to=r)
Rasters.bounds(cr)

bnd_matrix = hcat(map(collect, bds)...)
extrema(r)

using Shapefile
using GeoInterface
const GI = GeoInterface
@doc Shapefile.Polygon
GI.bbox(r)
pol = Shapefile.Polygon(;MBR=GI.bbox(cr))


fn=raw"C:\Users\chs72fw\Downloads\ee-chart.csv"
df = CSV.read(fn,DataFrame)
df = unique(df,1)
rename!(df,1=>:date)
a=df.date
[Date(d, dateformat"u d, y") for d ∈ a] 

using CSV, DataFrames, Dates

function eeread(filename)
    df = CSV.read(filename, DataFrame)
    df = unique(df, 1) # remove duplicate date rows
    rename!(df, 1 => :date)
    df.date .= [Date(d, dateformat"u d, y") for d in df.date]
    z = basename(filename)
    DataFrames.metadata!(df, "filename", z, style=:note)
    return df
end

# Call the function
df = eeread(raw"C:\Users\chs72fw\Downloads\ee-chart.csv")
df = selt(df,4)
dropmissing!(df)
dfm(df;fun=monmean,mode=:bar)
hydromon(df;col=2)
@edit hydro(df;col=2)
wa.hydromon(df;)



using DataFrames
using StatsPlots
using Dates
using Statistics

function monc(x::DataFrame)
    dmean = copy(x)
    s = map(Symbol, filter(x -> !occursin("date", x), names(dmean)))
    dmean[!, :month] = month.(df[!,:date]);
    select!(dmean, Not(:date))
    # Calculate monthly mean and confidence interval using combine and groupby
    dmean = DataFrames.combine(DataFrames.groupby(dmean, :month), 
        s .=> (dmean -> (mean(dmean), 1.96 * std(dmean) / sqrt(length(dmean)))) 
        .=> s)
    
    dmean.mn .= [first(only(x)) for x in eachrow(dx[!, Not(1)])]
    dmean.mx .= [last(only(x)) for x in eachrow(dx[!, Not(1)])]
    return dmean
end

@vv "um maximum"
dx = monc(df)

msx  = [minimum maximum]
cols = names(df)[Not(1)]
DataFrames.combine(df, cols .=> msx) 
mn = [first(only(x)) for x in eachrow(dx[!, Not(1)])]
mx = [last(only(x)) for x in eachrow(dx[!, Not(1)])]

using DataFrames
using StatsPlots
using Dates
using Statistics

function monc(x::DataFrame)
    dmean = copy(x)
    s = map(Symbol, filter(x -> !occursin("date", x), names(dmean)))
    dmean[!, :month] = month.(dmean[!,:date])
    select!(dmean, Not(:date))

    # Calculate monthly mean and confidence interval using combine and groupby
    dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
        s .=> (dmean -> (mean(dmean), 1.96 * std(dmean) / sqrt(length(dmean)))) 
        .=> s)

    # Create separate columns for mean, lower bound, and upper bound
    for col in s
        dmean[!, string(col)*"_mean"]  .= [first(only(x)) for x in eachrow(dmean[!, col])]
        dmean[!, string(col)*"_lower"] .= [first(only(x)) - last(only(x)) for x in eachrow(dmean[!, col])]
        dmean[!, string(col)*"_upper"] .= [first(only(x)) + last(only(x)) for x in eachrow(dmean[!, col])]
    end

    # # Drop the original columns representing the mean and standard deviation
    # select!(dmean, Not.(:month, s...))

    return dmean
end

dx = monc(df)


@df dx StatsPlots.plot(
            Vector(dx[!, 1]),
            cols(dx[!, Not(1:2)]|>propertynames),
            ribbon = true,  # Add ribbon for confidence interval
            #legend = leg,
            #title = ti,
            #seriestype = mode
        )

#mn = [ monthabbr(x) for x in unique(month.(df.date)) ]


#xticks!(st, mn)

theme(:dao)
function laiplot(df::DataFrame)
    ti = try
        DataFrames.metadata(df) |> only |> last |> basename
    catch
        @warn "No basename in metadata!"
        ti = raw""
    end
    dmean = copy(df)
    s = map(Symbol, filter(x -> !occursin("date", x), names(dmean)))
    dmean[!, :month] = month.(dmean[!,:date])
    select!(dmean, Not(:date))
    # Calculate monthly mean and confidence interval using combine and groupby
    dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
        s .=> (dmean -> (mean(dmean), 1.96 * std(dmean) / sqrt(length(dmean)))) 
        .=> s)

    # Create separate columns for mean, lower bound, and upper bound
    for col in s
        dmean[!, string(col)*"_mean"]  .= [first(only(x)) for x in eachrow(dmean[!, col])]
        dmean[!, string(col)*"_lower"] .= [first(only(x)) - last(only(x)) for x in eachrow(dmean[!, col])]
        dmean[!, string(col)*"_upper"] .= [first(only(x)) + last(only(x)) for x in eachrow(dmean[!, col])]
    end

    plt = @df dx plot(:month, :Lai_mean, 
        ribbon=(dx.Lai_lower, dx.Lai_upper), 
        labels = ti[1:4],   # "Lai_mean", 
        xlabel = "", #"Month", 
        ylabel = "",  #"Lai", 
        title = ti, #"Lai Mean with Ribbon",
        legend = :topleft,
        xticks = (1:12, 
        [ monthabbr(x) for x in unique(dx.month) ]),
        xrotation = 35
        )    
    return plt
end


laiplot(df)
df=dfr(r"qg")
laiplot(df)

wa.dfrib(selt(df, 6))
wa.dfrib(r"qges")
wa.dfrib(r"sb05";col=3)


wa.dfrib(r"wind")


using DataFrames
using StatsPlots
using Dates
using Statistics

function monc(df::DataFrame; confidence_level=0.95)
    # Make a copy of the input DataFrame
    dmean = copy(df)
    
    # Extract the columns of interest
    columns = filter(x -> !occursin("date", x), names(dmean))
    
    # Add a 'month' column based on the 'date' column
    dmean[!, :month] = month.(dmean[!, :date])
    
    # Drop the 'date' column
    select!(dmean, Not(:date))
    
    # Calculate monthly mean and confidence interval using combine and groupby
    dmean = DataFrames.combine(DataFrames.groupby(dmean, :month)) do group
        result = DataFrame(month = group.month)
        for col in columns
            values = group[!, col]
            mean_value = mean(values)
            confidence_interval = (mean_value, quantile(Normal(), 0.5 + confidence_level / 2) * std(values) / sqrt(length(values)))
            mv = string(col) * "_mean"
            lb = string(col) * "_lower"
            ub = string(col) * "_upper"
            result.mean .= mean_value
            result.lb .= confidence_interval[1] - mean_value
            result.ub .= confidence_interval[2] - mean_value
        end
        return result
    end

    return dmean
end


df

dm = monc(df)

plot(dm.month, dm.mean, 
    ribbon=(dm.lb, dm.ub), labels = "Lai_mean", 
    xlabel = "Month", ylabel = "Lai", title = "Lai Mean with Ribbon", legend = :topleft, xticks = (1:12, [ monthabbr(x) for x in unique(dm.month) ]), xrotation = 35)

#subset prec ########
using GeoDataFrames
const gdf = GeoDataFrames
fn = raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\pre_genRE\Rhine basin.kml"
cdof(fn)
rh = gdf.read(fn)
rh.geometry
#main catchment from HydroSHEDS
vgpy("HydroSHEDS")
vgjl("2070428090")

fn = raw"D:\Hydrologie\HydroATLAS_2019\basins\BasinATLAS_v10_ufra\BasinATLAS_ufra_lev07.shp"
ab=GeoDataFrames.read(fn)
saale = @rsubset ab :2==2070428090
saale.geometry|>plot

# fn = raw"D:\Hydrologie\HydroATLAS_2019\basins\BasinATLAS_v10_ufra\BasinATLAS_ufra_lev06.shp"
# ab=GeoDataFrames.read(fn)
# findindf(ab,"206042")
#main = @rsubset ab :2==2060420340
# main = findindf(ab,"207042")
# #main = @rsubset ab :2=
# main.geometry|>plot
# vgjl("annotations = ")
vgjl("centroi")

###plots with annotations.
fn = raw"D:\Hydrologie\HydroATLAS_2019\basins\BasinATLAS_v10_ufra\BasinATLAS_ufra_lev06.shp"
ab=GeoDataFrames.read(fn)
ab.geometry|>plot
for (i, pt) in enumerate(ab.geometry)
    pt = ArchGDAL.centroid(ab.geometry[i])
    x = ArchGDAL.getx(pt, 0)
    y = ArchGDAL.gety(pt, 0)
    name = ab.HYBAS_ID[i]
    annotate!(x, y, text(name, 8, :black, :bottom, :left))
end
Plots.plot!()


main = @rsubset ab :HYBAS_ID==2060429670
ofn = "D:/Relief_DGMs/FABDEM/wasim/lowres/main.geojson"
gdf.write(ofn,main)
ofn = "D:/Relief_DGMs/FABDEM/wasim/lowres/saale.geojson"
gdf.write(ofn,saale)

lpt = raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\pre_genRE"
ofn = joinpath(lpt,"saale.geojson")
gdf.write(ofn,saale)
ofn = joinpath(lpt,"main.geojson")
gdf.write(ofn,main)

###dl link.
#https://opendap.4tu.nl/thredds/catalog/data2/uuid/c875b385-ef6d-45a5-a6d3-d5fe5e3f525d/catalog.html?dataset=scanDatasets2/uuid/c875b385-ef6d-45a5-a6d3-d5fe5e3f525d/genRE_precipitation_hour_1.1.nc

#pre = Raster(joinpath(lpt,"genRE_precipitation_hour_1.1.nc");key=:pr,lazy=true)
#ds = NCDatasets.read(joinpath(lpt,"genRE_precipitation_hour_1.1.nc"))
#using RCall
using PyCall
lpt = raw"D:\remo\genRE"

using PyCall

@pyimport rasterio as rio
@pyimport rasterio.crs as rcrs
# Define the EPSG code
epsg_code = "EPSG:32632" #from cdo griddes
# Open the dataset
@pywith rio.open(joinpath(lpt, "genRE_precipitation_hour_1.1.nc"), "r+") as src begin
    # Set the CRS of the dataset
    src.set_crs(rcrs.CRS.from_string(epsg_code))
end



@pyimport xarray as xr
ds = xr.open_dataset(joinpath(lpt,"genRE_precipitation_hour_1.1.nc"))
# Resample the data to daily timesteps
ds_daily = ds.resample(time="D").sum("time")
# Save the resampled data to a new file
ds_daily.to_netcdf(joinpath(lpt, "genRE_daily.nc"))
#ds_daily = xr.open_dataset(joinpath(lpt,"genRE_daily.nc"))
@pyimport geopandas as gpd
using PyCall

@pyimport geopandas as gpd
@pyimport rasterio as rio
@pyimport rasterio.mask as rmask
@pyimport xarray as xr
@pyimport shapely.geometry as sgeom

# Open the dataset
ds_daily = xr.open_dataset(joinpath(lpt, "genRE_daily.nc"))

lpt = raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\pre_genRE"
ofn = joinpath(lpt,"saale.geojson")
# Load the GeoJSON file
gdf = gpd.read_file(ofn)
# Convert the GeoDataFrame to a list of shapes
shapes = [sgeom.shape(geometry) for geometry in gdf["geometry"]]
@pywith rio.open(joinpath(lpt, "genRE_precipitation_hour_1.1.nc")) as src begin
    # Print the CRS of the NetCDF file
    println(src.crs)
end


using PyCall

@pyimport rasterio as rio
@pyimport rasterio.warp as rwarp
@pyimport pyproj
@pyimport geopandas as gpd

# Open the raster dataset
@pywith rio.open(joinpath(lpt, "genRE_daily.nc")) as src begin
    # Get the CRS of the raster dataset
    src_crs = src.crs.to_string()
    # Open the GeoJSON file
    gdf = gpd.read_file(ofn)
    # Get the CRS of the GeoJSON file
    dst_crs = gdf.crs.to_string()
    # If the CRSs do not match, reproject the raster dataset
    if src_crs != dst_crs
        # Create a transformer
        transformer = pyproj.Transformer.from_crs(src_crs, dst_crs, 
        always_xy=true)
        # Reproject the raster dataset
        out_image, out_transform = rwarp.reproject(
            src.read(),
            src.transform,
            src_crs=src_crs,
            dst_crs=dst_crs,
            resampling=rwarp.Resampling.nearest,
            dst_transform=out_transform,
            dst_shape=src.shape,
            transformer=transformer
        )

        # Update the metadata to reflect the new CRS
        out_meta = src.meta.copy()
        out_meta.update(crs=dst_crs)

        # Save the reprojected raster dataset to a new file
        @pywith rio.open(joinpath(lpt, "genRE_daily_reprojected.nc"), "w", Dict{Symbol, Any}(map(Pair, keys(out_meta), values(out_meta)))) as dest begin
        dest.write(out_image)
        end
    end
    end
end

# Open the dataset with rasterio
@pywith rio.open(joinpath(lpt, "genRE_daily.nc")) as src begin
    out_image, out_transform = rmask.mask(src, shapes, crop=true)
    out_meta = src.meta
end

# Update the metadata to reflect the new shape
out_meta.update({"driver": "NetCDF",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})




using GeoDataFrames
using Rasters
import NCDatasets
cd("d:/remo/genRE")
lpt = raw"D:\remo\genRE"
cdof(fn)
ls()
#EPSG:32632"
pre = Raster("genRE_precipitation_hour_1.1.nc";key=:pr,
        lazy=true,crs=EPSG(32632),mappedcrs=EPSG(32632)) 

rh = GeoDataFrames.read(raw"D:\Wasim\main_basin.geojson")
rh.geometry
import ArchGDAL
#rh2 = lpro(rh.geometry; dst=EPSG(25832), src=EPSG(4326))
out = ArchGDAL.reproject(rh.geometry,
                    GeoFormatTypes.EPSG(4326),GeoFormatTypes.EPSG(32632))

@doc Where
#A[X(Where(x -> x > 15)), Y(Where(x -> x in (19, 21)))]
pre
p2 = pre[Ti(Where(x -> year(x) == 1999))]
ot = crop(p2; to=out|>first)
setup()
ot = mask_trim(p2,out|>first)
plot(ot[Ti=1])

ap = wa.project(ot[Ti=1];dst="EPSG:25832",src="EPSG:32632")
#contour(ap)
plot(ap)
@doc Near
#plot(ot[X(Near(225650),Y(Near(51294e6)))])
plot(ot[X=1,Y=4])

pre = Raster("genRE_daily.nc";key=:pr,
        lazy=true,crs=EPSG(32632),mappedcrs=EPSG(32632)) 
p2 = pre[Ti(Where(x -> year(x) == 2010))]
plot(p2[Ti=2])
bar(p2[X=34,Y=4])

ot = mask_trim(p2,out|>first)

using PyCall
@pyimport rasterio as rio
@pyimport rasterio.crs as rcrs
# Define the EPSG code
epsg_code = "EPSG:32632"
# Open the dataset
@pywith rio.open(joinpath(lpt, "genRE_daily.nc"), "r+") as src begin
    # Set the CRS of the dataset
    src.set_crs(rcrs.CRS.from_string(epsg_code))
end


z = out|>first
using GeoInterface
const GI = GeoInterface
GI.bbox(z)|>collect
#|>collect

ulx = 5.238845967940431e6
uly = 1.7250161769856259e6
lrx = 5.514078458420032e6
lry = 1.206493525543713e6
input_file = "genRE_daily.nc"
output_file = "genRE_crop.nc"
run(`gdal_translate -of NetCDF -projwin $ulx $uly $lrx $lry $input_file $output_file`)
lpt = raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\pre_genRE"
cd(lpt)
input_file = "genRE_precipitation_hour_1.1.nc"
output_file = "genRE_hour_crop.nc"
run(`gdal_translate -of NetCDF -projwin $ulx $uly $lrx $lry $input_file $output_file`)
#ERROR 1: netcdf error #-62 : NetCDF: One or more variable sizes violate format constraints 

#centroid(pre)

using RCall
tr=rimport("terra")
pr=tr.rast("genRE_daily.nc")
tr.crs("EPSG:32632")
#tr.summary(pr)
#pr.crs 
#rr=tr.project(pr,"EPSG:25832")

#p2[X=34,Y=4]
r = p2[Ti=16]
contourf(r)
plot!(out|>first,fillcolor=:transparent)
#ot = mask_trim(r,out|>first)

z = project(r;src="EPSG:32632",dst="EPSG:4326")
@doc lpro
rh = GeoDataFrames.read(raw"D:\Wasim\main_basin.geojson")
rh.geometry

plot(z)
plot!(rh.geometry,fillcolor=:transparent)
zc = crop(z; to=rh.geometry|>first) #this crops the raster to the extent of the polygon.
zc = wa.mask_trim(z,rh.geometry;pad=12) #so this really masks the raster.
plot(zc)
plot!(rh.geometry,fillcolor=:transparent)
@doc mask_trim

cd(raw"D:\remo\genRE")
pwd()
pre = Raster("genRE_daily.nc";key=:pr,
        lazy=true,crs=EPSG(32632),mappedcrs=EPSG(32632)) 
z = project(pre;src="EPSG:32632",dst="EPSG:4326") #takes looong
zc = crop(z; to=rh.geometry) #this crops the raster to the extent of the polygon.
plot(zc[Ti=99])
plot!(rh.geometry,fillcolor=:transparent)
Rasters.write(zc,"main_daily.nc")

plot(zc[Ti=90:99])
@vv "Where" #xs[Rasters.Where(xs .> msk)]
#: p2 = pre[Ti(Where(x -> year(x) == 1999))] 
su = zc[Ti(Where(x -> year(x) > 1999))]
dims(su)
typeof(su)|>cb
write("main_daily.tif",su;force=true) #errors
#rmlat()
plot(su[X=11,Y=1])
p2 = project(su,src="EPSG:4326",dst="EPSG:25832")
p2 = p2[Ti(Where(x -> year(x) < 2016))]
tst = p2[Ti(Where(x -> year(x) == 2014))]
plot(tst[Ti=1])
#plot!(rh.geometry,fillcolor=:transparent)
write("main_daily.nc",p2;force=true) #errors
fy = Raster(tst[Ti=1])

write("main_daily.asc",fy;force=true) 


x="D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
fl = CSV.read(x,DataFrame;limit=4)
xc = fl[2,5:end]|>collect
yc = fl[3,5:end]|>collect
pts = ArchGDAL.IGeometry[]
for i in 1:length(xc)
	pt = ArchGDAL.createpoint([xc[i],yc[i]])
	pt = ArchGDAL.reproject(pt,EPSG(25832),
		ProjString("+proj=longlat +datum=WGS84 +no_defs"))
	push!(pts,pt)
end
od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
# Convert observations to points
pnts = collect((ArchGDAL.getx(o,0), ArchGDAL.gety(o,0)) for o in (od.geometry) if !ismissing(o))  

dou = collect(Rasters.extract(su, pnts;atol=.1)|>skipmissing)
dx = DataFrame(dou)
GeoDataFrames.getgeometrycolumns(dx)
rename!(dx,2=>:pre)
#make station df loop
df = []
for i in 1:nrow(dx)
	k = dx.pre[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,:pre]))
end
df
ds = innerjoin(df..., on= :date, makeunique=true)
#nn = map(y->join(y[1:5],y[end],"_") ,string.(od.name))
#nn = map(y->hcat(y[1:5],y[end]) ,string.(od.name))
nn = map(y->(y[1:6]*"_"*y[end-1]) ,string.(od.name))
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
#map(x->replace(x,-Inf=>0.0),eachcol(ds))
for i in eachcol(ds)
    replace!(i,-Inf=>0.0)
end
wa.writewa("prec_ext.wa",ds)
op()
wa.hydro(ds;col=2)
dfm(ds;ann=false,fun=yrsum)
wa.dfrib(ds;col=3)

da = @rsubset ds year.(:date)==2016
dfp(da)

using PyCall
@pyimport rioxarray as rx
@pyimport geopandas as gpd
@pyimport xarray as xr
cd("d:/remo/genRE")
dataset_path = "genRE_precipitation_hour_1.1.nc"
###dataset_path = "genRE_daily.nc"
geojson_path = "main_basin_32632.geojson"
ds = xr.open_dataset(dataset_path)
ds = ds.rio.write_crs("EPSG:32632")
ds.crs
gdf = gpd.read_file(geojson_path)
#subset ds by year > 1999
ds = ds.sel(time=ds.time.dt.year > 1999)
#ds.lat
ds["pr"]
@doc ds.resample
ds = ds.sel(time=ds.time.dt.year < 2016).resample(time="D",skipna=true).sum()
ds.to_netcdf("genRE_daily_32632.nc") #2000-2016
#series = ds["pr"].resample(time="Y").sum().to_series()
@doc ds.rio.reproject("EPSG:4326")

@pyimport numpy as np
#ValueError('resampling must be one of: Resampling.nearest, Resampling.bilinear, Resampling.cubic, Resampling.cubic_spline, Resampling.lanczos, Resampling.average, Resampling.mode, Resampling.max, Resampling.min, Resampling.med, Resampling.q1, Resampling.q3, Resampling.sum, Resampling.rms')
#da = ds.rio.reproject("EPSG:4326", resampling="Resampling.cubic")
ds.close()
#ds = xr.open_dataset("genRE_daily_32632.nc")
# ds = ds.rio.reproject("EPSG:4326") #mem error.

using GeoDataFrames
using Rasters
import NCDatasets
cd("d:/remo/genRE")
lpt = raw"D:\remo\genRE"
pre = Raster("genRE_daily_32632.nc";key=:pr,
        lazy=true,crs=EPSG(32632),mappedcrs=EPSG(32632))
# df = ncdf(pre)
# wawrite(df,"mean_daily_32632.wa")
# baryrsum(df)
p2 = pre[Ti(Where(x -> year(x) == 2014))]
#p2 = wa.project(p2;src="EPSG:32632",dst="EPSG:4326")
src="EPSG:25832"
dst="EPSG:4326"
imethod="bilinear"
flags = Dict(
    "s_srs"=>src,
    "t_srs"=>dst,
    "r"=>imethod)
#replace_missing!(p2,0)
rs = Rasters.warp(p2,flags)
# using PyPlot
# pygui(true)
using Makie,CairoMakie
rh = GeoDataFrames.read(raw"D:\Wasim\main_basin.geojson")
rh.geometry

Rasters.rplot(rs[Ti=31])
#Rasters.write("pre2014.tif",rs)
@vv "to Geot"

using GeoStats
import CairoMakie as Mke
gri = CartesianGrid(size(rs))
val = rs.data
findmax(val)
findmax(rs)
table = (; Z = rs.data|>vec)
typeof(table)
img = georef(table, gri)
viewer(img)
plot(rs[Ti=99])
#plot!(rh.geometry,fillcolor=:transparent)


import ArchGDAL
setup()
#using GeoIO

pwd()
ls()
k="pre_daily.nc"
k="daily_4326.nc" #R crop
#rs = Raster(k,dims=(X,Y),key=:pre)
rs = Raster(k)
df = ncdf(rs)
@cmk
cmk.tsp(df)
cmk.tsbar(df)
Rasters.rplot(rs[Ti=101])
@doc Rasters.rplot(rs[Ti=101];plottype=Makie.Surface)
Rasters.rplot(rs[Ti=101];plottype=Makie.Surface)
Rasters.rplot(rs[Ti=101];plottype="Surface")
Rasters.rplot(rs[Ti=101];plottype=Makie.Contourf,colorrange=(1,30))

#Rasters.rplot(rs[Ti=101];plottype=Makie.Contourf,colorbar_padding=0.1,axistype="LScene")

@edit cmk.tsbar(df)
vec(Matrix(df[:, Cols(r"43")]|>dropmissing))|>Mke.scatter
@cmk
cmk.tspblack(df)
cmk.cloudplot2(df)

using NCDatasets
using Makie
# Open the NetCDF file
ds = Dataset("pr-obs.nc")
# Get the "pr" variable at t=101
pr = ds["pr"][:, :, 101,]

# Close the dataset
close(ds)
# Scale factor
scale_factor = 10
# Scale the z-values
pr_scaled = pr .* scale_factor
pr_scaled = replace(pr_scaled,missing=>-1)
# Create a surface plot
scene = Makie.surface(pr_scaled)
@doc Makie.surface
#prev = pr_scaled'
pr_scaled = pr .* 25
prev = replace(pr_scaled',missing=>0)
Makie.surface(prev;colorscale=Makie.pseudolog10)

#mkdir(raw"D:\Wasim\regio\out\vt\v2")

rs = Raster("pr-obs.nc")
r = "D:/remo/cordex/eobs/v28/fg_ens_mean_v28.0e_rcm.nc"|>Raster
r[Ti=1]|>plot
r[Ti=1]|>Makie.surface
import GeoInterface as GI
GI.bbox(r)
r[Ti=1]|>size

# Extent in the X and Y dimensions
extent_x = (9.049900000000001, 11.9499)
extent_y = (49.0499, 50.6499)
# Size of the raster
size_x, size_y = r[Ti=1]|>size
# Calculate the cell size
cell_size_x = (extent_x[2] - extent_x[1]) / size_x
cell_size_y = (extent_y[2] - extent_y[1]) / size_y
(cell_size_x, cell_size_y)

function getcs(r::Raster)
    # Extent in the X and Y dimensions
    extent_x = GI.bbox(r)[1]
    extent_y = GI.bbox(r)[2]
    # Size of the raster
    size_x, size_y = size(r)[1:2]
    # Calculate the cell size
    cell_size_x = (extent_x[2] - extent_x[1]) / size_x
    cell_size_y = (extent_y[2] - extent_y[1]) / size_y
    println("cell size x: $cell_size_x", " cell size y: $cell_size_y")
    return (cell_size_x, cell_size_y)
end

z=Raster("pr-obs.nc",key=:pr)
getcs(z)
using PyCall
@pyimport xarray as xr
fn=raw"D:\remo\genRE\genRE_daily_32632.nc"
ds = xr.open_dataset(fn)
o=ds["pr"].isel(time=101).squeeze() 
o.plot.surface() 

ds["pr"].mean(["x","y"]).groupby("time.year").sum("time").plot()
#looks good
df = ds["pr"].mean(["x","y"]).groupby("time.year").sum("time").to_dataframe()
#df2 = wa.pydf_to_julia(df) #errors.
df.iloc[:,1]
df.reset_index(inplace=true)
#col_names = string.(df.columns)  # Get the column names from the Python DataFrame
col_names = df.columns
col_arrays = [convert(Array, df[col]) for col in col_names]
cas = [convert(Vector{Float64}, z) for z in col_arrays]
julia_df = DataFrame(Symbol(col) => arr for (col, arr) in zip(col_names, cas))

od = pydf(df)
df[col_names[1]]
o = convert(Array, df[col_names[1]])
sdo = convert(Vector{Float64},o)
sdo = convert(Vector{Int64},o)
kd = pydf_to_julia(df)
adf = pydf(df)
select!(adf, Not(2))
using Makie
import CairoMakie as Mke
barplot(adf.year,adf.pr)
Mke.scatter(adf.year,adf.pr)

using Makie
import CairoMakie as Mke

# Create the bar plot
fig = Figure()
ax = Axis(fig[1, 1])
bars = barplot!(ax, adf.year, adf.pr)
# Add annotations
Mke.annotations!(ax, string.(round.(adf.pr, digits=2)), 
    #Point.(zip(adf.year .- 1, adf.pr));
    Point.(zip(adf.year, fill(500, length(adf.year))));
    #Point.(zip(adf.year, fill(mean(adf.pr), length(adf.year))));
    rotation = 45
    )
    
# Display the figure
fig

mean(adf.pr)

s = string.(ds.crs)
s = replace.(s,r"\s+"=>" ")
#s = split(s,"\n")[4:end]
s = split.(s,":")
#s = strip(s)
#DataFrame(s,:auto)
s = join(s,"\n")

k =raw"D:/remo/genRE/pr-dly.nc"



import JSON
fn="D:/Wasim/Tanalys/DEM/brend_fab/out/m4/qg.json"
d = map(DataFrame,open(fn, "r") do io;JSON.parse(io);end)
vs = filter(x->ncol(x)>1,d)
vs = vcat(vs...)
names(d)
#ns = filter(x->propertynames(x)==:column,d)
#ns = map(x->startswith(names(x),"col"),d)
ns = filter(x->ncol(x)==1,d)
names.(ns)
ns = filter(x->propertynames(x)|>first==:column,ns)
ns = vcat(ns...)
df = hcat(ns,vs)
f=Figure()
ax=Axis(f[1,1])
for z in  eachcol(vs[:,Not(:sum)])
       lines!(ax,vec(z))
       end
f


fn="D:/Wasim/Tanalys/DEM/brend_fab/out/m4/qg2.json"
d = map(DataFrame,open(fn, "r") do io;JSON.parse(io);end)
vs = filter(x->ncol(x)>3,d)
vs = vcat(vs...)
#ns = filter(x->propertynames(x)==:column,d)
#ns = map(x->startswith(names(x),"col"),d)
ns = filter(x->ncol(x)<3,d)
ns = vcat(ns...)
df = hcat(ns,vs)
f=Figure()
ax=Axis(f[1,1])
for z in eachcol(df[:,Not(Cols(1,2,end))])
    lines!(ax,vec(z),label=string(z))
    #CairoMakie.label!(ax,vec(z))
end
f

as = map(DataFrame,JSON.parsefile(fn))

fn = "D:/Wasim/Tanalys/DEM/brend_fab/out/m4/n.json"
open(fn, "r") do io;
    d=JSON.parse(io);
end
d = map(DataFrame,d)
vs = filter(x->ncol(x)>3,d)
vs = vcat(vs...)
ns = filter(x->ncol(x)<3,d)
ns = vcat(ns...)
df = hcat(ns,vs)
msn = names(df[:,Not(Cols(2,3,end))])
f=Figure();
ax=Axis(f[1,1],title="ndf");
#nm = names(df[:,Not(Cols(1,2,end))])
#for z in eachcol(df[:,Not(Cols(1,2,end))])
for z in eachcol(df[:,Not(Cols(2,3,end))])
    lines!(ax,vec(z),label = string(names(z)))
end
f[1, 2] = Legend(f, ax, "Leg title", framevisible = false)
f



#####################json parse. #########################
using JSON
using DataFrames
#invokation jqs tempfab.m4.2010 > t.json
function process_json(fn::String)
    # open(fn, "r") do io;
    #     d=JSON.parse(io);
    # end
    # d = broadcast(DataFrame,d)
    d = map(DataFrame,JSON.parsefile(fn))
    vs = filter(x -> ncol(x) > 3, d)
    df = vcat(vs...)
    #df = df[:, sort(names(df),rev=true)] #sort columns
    df = df[:, sort(names(df),rev=true)] #sort columns
    #ns = filter(x -> ncol(x) < 3, d)
    #ns = vcat(ns...)
    #df = hcat(ns, vs)
    return df
    # if j isa Dict
    #     d = DataFrame(j)
    #     return d
    # else
    #     d = map(DataFrame, j)
    #     vs = filter(x -> ncol(x) > 3, d)
    #     vs = vcat(vs...)
    #     ns = filter(x -> ncol(x) < 3, d)
    #     ns = vcat(ns...)
    #     df = hcat(ns, vs)
    #     return df
    # end
    
end
fn = "D:/Wasim/Tanalys/DEM/brend_fab/out/m4/n.json"
fn = "D:/Wasim/Tanalys/DEM/brend_fab/out/m4/t.json"
ks = process_json(fn)

fn = "D:/Wasim/sinn/out/v2/fl.json"
function jsread(fn::String)
    d = map(DataFrame,JSON.parsefile(fn))
    vs = filter(x -> ncol(x) > 3, d)
    df = vcat(vs...)
    return df[:, sort(names(df),rev=true)]
end

ks = jsread(fn)

total_minutes = div(1932,60) * 32
hours = div(total_minutes, 60)
minutes = total_minutes % 60 #calculates the remaining minutes.
hours, minutes


"D:/Wasim/regio/out/rc200/x22/re-0/"|>cd
pwd()  
fn = raw"D:\ClimateExplorer\precipitation\iknmi_radar_daily_Rhine.nc"
da = Raster(fn)
#df = wa.ncdf(da)
plot(da)

ti = Rasters.lookup(da,1)
df = DataFrame(date = ti, pr = da.data)
DataFrames.metadata!(df, "filename", fn, style=:note);
cmk.tsp(df)
cmk.tsbar(df)
plot(df.pr)

cd("D:/Wasim/Tanalys/DEM/brend_fab/out/c0")
qgk()
x = @gl "temper"
cmk.mkemon(x;col=1)
x = @gl "qbas"
cmk.mkemon(x;msk=true)
x = @gl "qs"
cmk.dfp(x)


"D:/Wasim/Tanalys/DEM/brend_fab/out/w3/penman/re/"|>cd
#xr to mke heat.
using PyCall
@pyimport xarray as xr
fn = @nco "wind"
ds = xr.open_dataset(fn)
#ds["wind"].isel(t=1).squeeze() |> plot
k = ds.keys()|>collect|>first
A = ds[k].isel(t=0).squeeze()[:values]
Mke.heatmap(A)
A_rev = reverse(A, dims=1)
A_trans = transpose(A_rev)
Mke.heatmap(A_trans)
# Reverse the transposed matrix along the first dimension (equivalent to terra::flip in R)
A_flip = reverse(A_trans, dims=1)
Mke.heatmap(A_flip)
r=Raster(fn)
r[Ti=1]|>Mke.heatmap

#thats it!
Mke.heatmap(reverse(A, dims=2))

"""
uses xarray to read netcdf and makie to plot.
using PyCall
@pyimport xarray as xr
"""
function mkheat(x::Union{String,Regex};msk=true,layer=0)
    if isa(x,Regex)
        #fn = filter(z->endswith(z,"nc"),glob(x))[1]
        #ne = join([x,"nc","\$"],"+.")
        fn = filter(z->endswith(z,"nc"),Grep.grep(x,readdir()))[1]        
    else
        fn = x
    end
    ds = xr.open_dataset(fn)
    k = ds.keys()|>collect|>first
    A = ds[k].isel(t=layer).squeeze()[:values]
    A = reverse(A, dims=2) #very important!
    if msk
        A = map(x -> x <= 0 ? NaN : x, A)
    end
    fig = Figure(
        #size=(800, 600), 
        fontsize=22);
    ti = replace(fn,".nc"=>"","_"=>" ")
    axs = Axis(fig[1,1],title=ti, aspect=1, 
        xlabel="x", 
        ylabel="y")
    #replace NaN values with missing
    B = replace(A, NaN => missing)
    # Find rows and columns where all values are missing
    xm = .!all(ismissing, B, dims=2)
    ym = .!all(ismissing, B, dims=1)
    # Remove rows and columns where all values are missing
    A = B[xm[:], ym[:]]

    p1 = heatmap!(A,colormap=(:turbo, 0.9))
        contour!(axs, A; color=:black) #, levels=lscale
        Colorbar(fig[1, 2], p1, 
            width=20, ticksize=20, 
            tickalign=1)
        xmax,ymax = size(A)
        limits!(axs, 1, xmax, 1, ymax)
    # hideydecorations!(axs, grid=true, ticks=false)
    # hidexdecorations!(axs, grid=true, ticks=false)
    return fig
end

mkheat(r"tem";msk=true)
mkheat(r"soil";layer=2,msk=false)
mkheat(r"^vap";msk=true)
mkheat(r"vap";msk=true)

#YES
ds = xr.open_dataset("tsoilfab__stack.2015.nc")
k = ds.keys()|>collect|>first
A = ds[k].isel(t=3).squeeze()[:values]
B = replace(A, NaN => missing)
# Find rows and columns where all values are missing
xm = .!all(ismissing, B, dims=2)
ym = .!all(ismissing, B, dims=1)
# Remove rows and columns where all values are missing
C = B[xm[:], ym[:]]
C = reverse(C, dims=2)
Mke.heatmap(C)
Mke.surface(C*50)
Mke.surface(C*100)

b=Raster(B,(X,Y))
#A = replace(A, NaN => missing)
# Convert the matrix to a DataFrame
df = DataFrame(B,:auto)
# Remove rows with missing data
df = dropmissing(df, disallowmissing=false)
# Convert the DataFrame back to a matrix
A = Matrix(df)
@pyimport numpy as np
# Get the array from the xarray Dataset
A = ds[k].isel(t=0).squeeze().values
# Create a new array that doesn't include the NaN values
Matrix(A[.!np.isnan(A)])
dx_output = Matrix{Float32}(undef, size(A, 1),size(A,2))
A[.!np.isnan(A)]
A = A[.!np.isnan(A)]


"""
uses Rasters to read netcdf and makie to plot.
"""
function mkrheat(x::Union{String,Regex};msk=true,layer=1)
    if isa(x,Regex)
        #fn = filter(z->endswith(z,"nc"),glob(x))[1]
        #ne = join([x,"nc","\$"],"+.")
        fn = filter(z->endswith(z,"nc"),Grep.grep(x,readdir()))[1]        
    else
        fn = x
    end
    ds = Raster(fn)
    #A = ds[Ti=layer].data
    #A = ds[:, :, 1].data
    #mx = missingmask(ds)
        
    if msk
        zm = (ds .> float(0.001))
        A = Rasters.mask(ds; with=zm)
    else
        A = ds
    end
    
    A = A[:, :, layer].data
    mrows = map(x -> !all(ismissing, x),eachrow(A))
    mcols = map(x -> !all(ismissing, x),eachcol(A))
    A = A[mrows, mcols]
    #A = reverse(A,dims=1) #very important!
    A = transpose(A)
    A = reverse(A, dims=2) #very important!
    
    
    fig = Figure(
        #size=(800, 600), 
        fontsize=22);
    ti = replace(fn,".nc"=>"","_"=>" ")
    axs = Axis(fig[1,1],title=ti, aspect=1, 
        xlabel="x", 
        ylabel="y")
    

    p1 = heatmap!(A,colormap=(:turbo, 0.9))
        contour!(axs, A; color=:black) #, levels=lscale
        Colorbar(fig[1, 2], p1, 
            width=20, ticksize=20, 
            tickalign=1)
        xmax,ymax = size(A)
        limits!(axs, 1, xmax, 1, ymax)
    # hideydecorations!(axs, grid=true, ticks=false)
    # hidexdecorations!(axs, grid=true, ticks=false)
    return fig
end

mkrheat(r"gws";msk=false,layer=1)
mkrheat(r"tem";msk=true)
mkheat(r"tem";msk=true)

mkrheat(r"soi";msk=false)
@doc cmk.mkrheat(r"soi";)

mkheat(r"gws";msk=false)
fn = r"gws"|>glob|>first
Raster(fn)|>heatmap
"D:/Wasim/regio/out/rc200/x22/re3" |>cd
cmk.mkrheat(r"soi";msk=false)
cmk.mkrheat(r"soi";msk=true,layer=2)

#YES
ds = xr.open_dataset("tsoilrcm__stack.2015.nc")
k = ds.keys()|>collect|>first
A = ds[k].isel(t=3).squeeze()[:values]

ds[k].isel(x=3,y=6).squeeze().values|>plot
##same as contourf
ds[k].mean("t").values|>v->reverse(v, dims=2)|>plot
#mean plots.
ds[k].mean("x").mean("y").mean()
ds[k].mean("x").mean("y").values|>lines

# Calculate the mean values
mean_values = ds[k].mean("x").mean("y").values
# Create a line plot
fig = Figure()
ax = Axis(fig[1, 1],
xlabel = rich(
    "Layers",
    subscript(" of unsatzone model", color = :slategray)
),  title="Mean Soil temperature values",
    ytickformat = "{:.2f} °C")
lines!(ax, mean_values)
fig

fig = Figure()
ax = Axis(fig[1, 1],yreversed=true,
ylabel = rich(
    "Depth",
    subscript(" #Layers of unsatzone model", color = :slategray)
),  title="Mean Soil temperature values",
    yscale = Makie.logit,
    xtickformat = "{:.1f} °C")
lines!(ax, mean_values,1:size(mean_values)[1])
fig


#transform to 0..1
mean_values = ds[k].mean("x").mean("y").values
mean_values = (mean_values .- minimum(mean_values)) ./ (maximum(mean_values) - minimum(mean_values))
mean_values = Float64.(mean_values)
#mean_values = mean_values .* 0.998 .+ 0.001

using CairoMakie
f = Figure()
for (i, scale) in enumerate([identity, log10, log2, log, sqrt])
    #, Makie.logit
    row, col = fldmod1(i, 3)
    Axis(f[row, col], yscale = scale, title = string(scale),
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(5))
    #lines!(range(0.01, 0.99, length = 200))
    lines!(mean_values .+ 0.0001) # Add a small constant to mean_values
end
f


##read ascii
cd("D:/Relief_DGMs/FABDEM/wasim/lr/")
k = Raster("rcm.kol",missingval=-9999)
k = read(Raster("rcm.aq1",missingval=-9999.0f0)) 
extrema(k)

k = read(Raster("rcm.aq1"))
Rasters.replace_missing!(k,-9999.0f0)
#k[k .== -9999.0f0] 
#da = replace.(k,-9999.0f0=>missing)

using ArchGDAL
# Open the dataset
dataset = ArchGDAL.read("rcm.aq1")
# Get the first band
band = ArchGDAL.getband(dataset, 1)
dx = permutedims(band, (2, 1))
A = dx.a
dx = reverse(A, dims=1)
B = replace(dx, -9999.0 => missing)
xm = .!all(ismissing, B, dims=2)
ym = .!all(ismissing, B, dims=1)
A = B[xm[:], ym[:]]
k = Raster(A,(X,Y))
extrema(k)


heatmap(k)
mean(k)

m = Raster(k.data;missingval=-9999)
r = Rasters.rebuild(k,missingval=minimum(k))
extrema(r)

function agread(x::String)
    dataset = ArchGDAL.read(x)
    # Get the first band
    band = ArchGDAL.getband(dataset, 1)
    dx = permutedims(band, (2, 1))
    A = dx.a
    dx = reverse(A, dims=1)
    B = replace(dx, -9999.0 => missing)
    xm = .!all(ismissing, B, dims=2)
    ym = .!all(ismissing, B, dims=1)
    A = B[xm[:], ym[:]]
    k = Raster(A,(X,Y))
    println(extrema(k))
    return k
end

k = cmk.agread("rcm.ezg")
heatmap(k)


###THIS IS RAW dat!
lk="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/"
x="Frnkische-24409003-Wolfsmuenster.was"

df = dfr(joinpath(lk,x))
using CairoMakie


rename!(df,1=>:Flows)
dropmissing!(df)
cmk.tsp(df)
using DataFrames
using Dates
# Assuming df is your existing DataFrame with a Date column named :date and a Flow column named :Flows
# Create a new column :year in df
df[!, :year] = year.(df[!, :date])
# Group by year
grp = groupby(df, :year)
# Initialize an empty DataFrame to store the results
result_df = DataFrame(year = Int[], MaxQ = Float64[], MinQ = Float64[], MeanQ = Float64[], Q90 = Float64[], Q10 = Float64[], Median = Float64[])
#result_df = []
for g in grp
    year_val = first(g.year)  # Extract the year value from the group
    
    # Calculate metrics and store them in a NamedTuple
    metrics = (
        MaxQ = maximum(g[:, :Flows]),
        MinQ = minimum(g[:, :Flows]),
        MeanQ = mean(g[:, :Flows]),
        Q90 = quantile(g[:, :Flows], 0.9),
        Q10 = quantile(g[:, :Flows], 0.1),
        Median = median(g[:, :Flows])
    )

    # Append a row to the result_df
    push!(result_df, (year=year_val, metrics...))
end

# Print the result_df
show(result_df, allcols=true)

using CairoMakie
# Create a new figure
f = Figure();
# Create an axis
ax = Axis(f[1, 1],
    ylabel = rich(
        "Flow",
        subscript(" #raw data", color = :slategray)
    ),  title="Saale",
        yscale = log10,
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(5),
        ytickformat = "{:.2f} [m³/s]")


# Create a scatter plot of the mean values
lines!(ax, result_df.year, result_df.MeanQ, color = :blue)
scatter!(ax, result_df.year, result_df.MaxQ, color = :red)
scatter!(ax, result_df.year, result_df.MinQ, color = :red)
# Add a band representing the 10th and 90th percentiles
band!(ax, result_df.year, result_df.Q10, result_df.Q90, color = (:blue, 0.2))
# Set the axis labels
ax.xlabel = "Year"
#ax.ylabel = "Flow"
# Show the figure
f


df = dfr(joinpath(lk,x))

dm = cmk.monmean(df)


using DataFrames
using Statistics

# Assuming result_df is your existing DataFrame with a :year column and other columns for flow metrics

# Create a new column :month in result_df
df[!, :month] = month.(df[!, :date])
# Group by year and month
grp_monthly = groupby(df, [:year, :month])
# Initialize a new DataFrame to store the monthly means



lk="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/"
x="Frnkische-24409003-Wolfsmuenster.was"
df = dfr(joinpath(lk,x))
rename!(df,1=>:Flows)

df = filter(x->x.date >= Date(1980,1,1),df)
dropmissing!(df)

# Create a new column :year in df
df[!, :year] = year.(df[!, :date])
# Calculate hydrological year starting in October
#df[!, :hydro_year] = ifelse.(month.(df[!, :date]) >= 10, year.(df[!, :date]), year.(df[!, :date]) - 1)
df[!, :month] = month.(df[!, :date])
# Group by hydrological year and month
grp_monthly = groupby(df, [:year, :month])
# Initialize an empty DataFrame to store the overall monthly metrics
result_monthly_df = DataFrame(hydro_year = Int[], month = Int[], MaxQ = Float64[], MinQ = Float64[], MeanQ = Float64[], Q90 = Float64[], Q10 = Float64[], Median = Float64[])

# Iterate over groups
for g in grp_monthly
    hydro_year_val = first(g.year)
    month_val = first(g.month)
    
    # Calculate metrics and store them in a NamedTuple
    metrics = (
        MaxQ = maximum(g[:, :Flows]),
        MinQ = minimum(g[:, :Flows]),
        MeanQ = mean(g[:, :Flows]),
        Q90 = quantile(g[:, :Flows], 0.9),
        Q10 = quantile(g[:, :Flows], 0.1),
        Median = median(g[:, :Flows])
    )

    # Append a row to the result_monthly_df
    push!(result_monthly_df, (hydro_year=hydro_year_val, month=month_val, metrics...))
end

# Print the overall monthly metrics DataFrame
show(result_monthly_df, allcols=true)

pwd()
writedf(result_monthly_df,"result_monthly_df.wa")



df = result_monthly_df
grp = groupby(df, :month)
months = ["Januar", "Februar", "März", "April",
    "Mai", "Juni", "Juli", "August", "September",
    "Oktober", "November", "Dezember"]

ti = "s"
col = :Median  
col = :MaxQ
f = Figure()
Axis(f[1, 1], title = string(col), #ti*" Basin:"*string(col),
    yticks = ((1:12) ./ 4,  reverse(months)))
for i in 12:-1:1
    values = select(grp[i], col)|>Matrix|>vec
    d = density!(values, offset = i / 4,
                color = :x, colormap = :thermal, 
                #colorrange = (-5, 5),
                strokewidth = 1, strokecolor = :black)
    #Apply an absolute translation to the Scene
    CairoMakie.translate!(d, 0, 0, -0.1i) 
end
return f




#grp = groupby(df, :month)

begin
    # Create a new figure
    f = Figure();
    # Create an axis
    ax = Axis(f[1, 1],
        xtickwidth = 20,
        xticksmirrored = true,
        xminorticks = IntervalsBetween(5),
        xminorticksvisible = true, 
        xminorgridvisible = true,
        xminortickwidth = 1,
        yticks = (1:12, reverse(months)) 
        )
    # Get the number of groups
    n = length(grp)
    # Create a boxplot for each month
    #for i in Iterators.reverse(enumerate(grp))
    for i in n:-1:1
        mg = grp[i]
        boxplot!(ax, mg.month, mg.MaxQ,  
        orientation= :horizontal,
        #show_notch = true, 
        #color = (:red,.3))
        color = (:grey,.3))
        boxplot!(ax,mg.month, mg.MeanQ, 
        orientation= :horizontal,
            show_notch = true,   
            color = :grey)
        #   lines!(ax, mg.month, mg.MinQ, color = :grey)
        #   band!(ax, mg.month,
        #     mg.Q10, 
        #     mg.Q90, color = (:blue, 0.2))

    end
    # Set the axis labels
    ax.xlabel = "Abfluss [m³/s]"
    #ax.ylabel = "Flow"
    # Show the figure
    f
end

begin
    # Create a new figure
    f = Figure();
    # Create an axis
    ax = Axis(f[1, 1],
        #xticks = (1:12, months) 
        )
    # Create a boxplot for each month
    #for i in eachrow(df)
          lines!(ax, df.hydro_year, df.MinQ, color = :grey)
          band!(ax, df.hydro_year,
            df.Q10, 
            df.Q90, color = (:blue, 0.2))
    #end
    # Set the axis labels
    ax.xlabel = "Abfluss [m³/s]"
    #ax.ylabel = "Flow"
    # Show the figure
    f
end


# Code für wlf mit annotations. theme_latexfonts()
begin
    f = Figure();
    # Create an axis
    ax = Axis(f[1, 1],
        xtickwidth = 20,
        xticksmirrored = true,
        xminorticks = IntervalsBetween(5),
        xminorticksvisible = true, 
        xminorgridvisible = true,
        xminortickwidth = 1,
        yticks = (1:12, reverse(months)) 
        )
    # Get the number of groups
    n = length(grp)
    # Create a boxplot for each month
    for i in n:-1:1
        mg = grp[i]
        boxplot!(ax, mg.month, mg.MaxQ,  
        orientation= :horizontal,
        color = (:grey,.3))
        boxplot!(ax,mg.month, mg.MeanQ, 
        orientation= :horizontal,
            show_notch = true,   
            color = :grey)
        # Calculate the median
        median_value = round(median(mg.MeanQ), digits = 2)
        #text(x, y; text, kwargs...)
        # weil mg.month[end] ein vector ist!
        #median_value - 45
        CairoMakie.text!(-35, mg.month[end];
        text = string(median_value), color = :black)
        maxmedian = round(median(mg.MaxQ), digits = 2)
        CairoMakie.text!(250, mg.month[end];
        text = string(maxmedian), color = :black)
    end
    # Set the axis labels
    ax.xlabel = "Abfluss [m³/s]"
    CairoMakie.text!(-39, 0.0001; text = "QMean", color = :black)
    CairoMakie.text!(250, 0.0001; text = "QMax ",  color = :black)
    
    # \nMedian"
    # \nMedian"
    #align=:center,
    #align=:center,
    
    # Show the figure
    f
end


boxplot(mg.month, mg.MaxQ)
median_value = median(mg.MaxQ)
CairoMakie.text!(mg.month, median_value;
text = string(median_value), color = :black)





using CairoMakie
using Statistics





using CairoMakie, DataFrames
# Transform the month column to a categorical variable
result_monthly_df.month = categorical(result_monthly_df.month, ordered=true)
# Reorder the levels to start from October
#levels!(result_monthly_df.month, [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9])
# Create a new figure
#result_monthly_df.month = Vector(result_monthly_df.month)

f = Figure()
# Create an axis
ax = Axis(f[1, 1])
# Create a scatter plot of the mean values
scatter!(ax, result_monthly_df.month, result_monthly_df.MeanQ, color = :blue)
# Add a band representing the 10th and 90th percentiles
band!(ax, result_monthly_df.month, 
    result_monthly_df.Q10, 
    result_monthly_df.Q90, color = (:blue, 0.2))
# Set the axis labels
ax.xlabel = "Month"
ax.ylabel = "Flow"
# Show the figure
f



# Create a new figure
f = Figure()
# Create an axis
ax = Axis(f[1, 1])
# Create a scatter plot of the mean values
scatter!(ax, result_monthly_df.month|>collect, result_monthly_df.MeanQ, color = :blue)
    
f

using CairoMakie, Random, Distributions
Random.seed!(13)
n = 3000
colors = Makie.categorical_colors(:spring, 8)[3:end]
CairoMakie.set_theme!(theme_latexfonts())
fig = Figure(resolution = (600, 400))
ax = Axis(fig[1,1]; palette = (; patchcolor = colors), 
    xticks = (1:7, ["cat 1", "A", "B", "C", "D", "E", "F"]), 
    yticks = ([-5], ["cat 2"]), yticklabelrotation = π/2)
boxplot!(ax, fill(-5,n), rand(Distributions.Normal(0, 0.5), n); orientation=:horizontal, 
    whiskerwidth = 1, width = 2, color = (:orange, 0.95), 
    whiskercolor = :red, mediancolor = :yellow, markersize = 8, 
    strokecolor = :black, strokewidth = 1, label = "horizontal")
boxplot!(ax, fill(1,n), rand(Distributions.Normal(1,  3), n); whiskerwidth = 1, 
    width = 0.5, color = :dodgerblue, whiskercolor = :dodgerblue, 
    mediancolor = :white, markersize = 5, strokecolor = :white, 
    strokewidth = 1, label = "vertical")
for i in 2:7
    boxplot!(ax, fill(i,n), rand(Distributions.Normal(rand(-2:5), 2*rand() + 0.3), n); 
        whiskerwidth = 1, width = 0.35)
end
axislegend(ax, position = :lt)
display(fig)



@vv "xtickformat"


using CairoMakie
odeSol(x, y) = Point2f(-x, 2y) # x'(t) = -x, y'(t) = 2y
fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", backgroundcolor = :black)
streamplot!(ax, odeSol, -2 .. 4, -2 .. 2, colormap = Reverse(:plasma),
    gridsize = (32, 32), arrow_size = 10)
fig

odeSol(123,123)
@doc Point2f




using PyCall

function drop_crs_xarray(input_file::String, output_file::String)
    # Import necessary Python modules
    xr = pyimport("xarray")

    # Load the NetCDF file
    ds = xr.open_dataset(input_file)

    # Check if the 'crs' variable exists as a data variable or coordinate variable
    if haskey(ds.variables, "crs") || haskey(ds.coords, "crs")
        # Drop the 'crs' variable
        ds = ds.drop_vars("crs", errors="ignore")
        println("Variable 'crs' dropped successfully.")
    else
        println("Variable 'crs' does not exist in the input file.")
    end

    # Save the modified dataset to a new NetCDF file
    ds.to_netcdf(output_file)
end

cd("D:/remo/genRE/")
drop_crs_xarray("./pre-gencrs.nc","./pre-gen.nc")
ds = xr.open_dataset("./pre-gencrs.nc")
haskey(ds.variables, "crs")
haskey(ds.coords, "crs")
ds = ds.transpose("x","y","t","pr","crs")

# Rename the dimensions
ds = ds.rename(Dict("longitude" => "x", "latitude" => "y", "time" => "t"))
# Add the new dimensions 'pr' and 'crs'
# Transpose the dimensions
#ds = ds.transpose("x","y","t","pr","crs")
output_file="./pre-gen.nc"
ds.to_netcdf(output_file)


r = Raster("D:/Wasim/regio/out/rc200/v2/temprcm.2012.nc")
z = r.data[:,:,1]
reverse!(z,dims=1)
z = z'
ex = extrema(z|>skipmissing)
lscale = ex[1]:10:ex[2] # Adjusted levels
fig = Figure(
            #size=(800, 600), 
            fontsize=22);
axs = Axis(fig[1,1], aspect=1, xlabel="x", ylabel="y")
      
xmax = findmax(size(z))[1] #gets the index of the max value
xmin = findmin(size(z))[1] #gets the index of the max value
streamplot!(axs, z)     
p1 = streamplot!(axs, 
z, ex[1] .. ex[2], ex[1] .. ex[2], 
colormap = Reverse(:plasma),
gridsize = (32, 32), arrow_size = 10)

vy =r.data[:,:,1]|>vec|>skipmissing
streamplot(vy)     
testField(x, y) = Point2f(-x, 2y) # x'(t) = -x, y'(t) = 2y

x,y = size(z)
x,y,z = r.data[:,1,1],r.data[1,:,1],r.data[1,1,:]

fig = Figure(fontsize = 22, fonts = (;regular="CMU Serif"));
ax = fig[1, 1] = Axis(fig, xlabel = L"x", ylabel = L"y")
streamplot!(ax,  colormap = Reverse(:plasma),
    gridsize = (32, 32), arrow_size = 10)                
fig


fn=raw"D:\Wasim\regio\out\rc200\x22\loc5\fspin\T_Lower_Boundary_Condition_000.nc"
Main.cmk.mkrheat(fn;msk=false)
cdof(fn)
Main.cmk.mkrheat(r"^temp")
Main.cmk.mkrheat(r"low"i)
mkrheat(r"sb05"i;msk=true,mskval=0.1)

fn="d:/Wasim/regio/out/r5/sb05rcm_1100.mit.nc"
fn="d:/Wasim/regio/out/r5/sb05rcm_1000.mit.nc"
cmk.mkrheat(fn;mskval=.1,umask=.44)
cmk.mkrheat(fn;msk=true,layer=1)
@doc cmk.mkrheat

agheat(fn)
fn="d:/Wasim/regio/out/r5/windrcm.2017.nc"
agcont(fn)
agcont2(fn)



@pyimport xarray as xr
cd("D:/remo/qm/corgrids/jlcor/")
ob = xr.open_dataset("rh-cor.nc")


using Pkg
Pkg.activate("D:/remo/qm/qm/")
Pkg.status()
function loadcdo()
    pt="C:/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end
loadcdo()


fn="D:/Wasim/sinn/out/qg_bestlin"
sim=dfr(fn)
obs=dfr("D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt")
"D:/Wasim/sinn/out/p0"|>cd
ctf=raw"C:\Users\chs72fw\Documents\EFRE_GIS\Hydrologie\wasim_pest\app\pest\sinn100.pest"
df = dfroute()
z = innerjoin(sim,obs[!,Cols(df.obs .-3,:date)],on=:date)
select!(z, :date, :)
wawrite(z,"m-qoutjl")
dfl(z)
wa.hydromon(z)
hyeval(z)

#mamba list --show-channel-urls --explicit | awk '{print $1,$2}' | sort -k2 -rn
####################json parsing#####

fn="D:/Wasim/regio/out/rc200/x36/cor6/q.json"
#d = map(DataFrame,open(fn, "r") do io;JSON.parse(io);end)
d = JSON.parsefile(fn)|>DataFrame
select!(d,:file,:)
d = JSON.parsefile(fn)|>DataFrame|>y->select(y,:file,:)
dst = rst.stats("D:/Wasim/regio/out/rc200/x36/cor6/qi__rcm_Layer_1.sum.nc")
d
fn="D:/Wasim/Tanalys/DEM/brend_fab/out/m3/t.json"
d = JSON.parsefile(fn)|>DataFrame #|>y->select(y,:file,:)


####################06.02.2024#####


xr  = pyimport("xarray")
sh="d:/remo/qm/prec/simh.nc"
simh=xr.open_dataset(sh)
sp="d:/remo/qm/corgrids/pre/pre-cor.nc"
simp=xr.open_dataset(sp)
tl=raw"D:\remo\genRE\pre-REcor_raw.nc"
obsh=xr.open_dataset(tl)
pyjl.doyplot(simh,simp,obsh;tosum=true)

k="pre"
grouped_obsh = dx[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()
grouped_simh = simh[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()
grouped_simp = simp[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()

plt.rc("font", family="serif", serif=["cmr10"])
vr = uppercase(first(k))
plt.figure()
plt.plot(grouped_simh, label="\$$vr _{sim,h}\$")
plt.plot(grouped_simp, label="\$$vr _{sim,p}\$")
plt.plot(grouped_obsh, label="\$$vr _{obs,h}\$")
plt.title("Historical modeled and predicted versus observed $(uppercasefirst(k))")
plt.xlim(0, 365)
plt.gca().grid(alpha=0.3)
plt.legend()
plt.show()

Main.pyjl.doyplot(simh,simp,obsh;tosum=true)
#ValueError("'longitude' not found in array dimensions ('time', 'y', 'x')")

sh="d:/remo/qm/prec/simh.nc"
simh=xr.open_dataset(sh)
sp="d:/remo/qm/corgrids/pre/pre-cor.nc"
simp=xr.open_dataset(sp)
tl="D:/remo/genRE/pr-dly.nc"
obsh=xr.open_dataset(tl).rename(Dict("pr"=>"pre"))
pyjl.doyplot(simh,simp,obsh;tosum=true)

sh="d:/remo/qm/prec/simh.nc"
simh=xr.open_dataset(sh)
sp = "D:/remo/qm/corgrids/jlcor/pre-cor.nc" ##<-time is not a dimension
sp = "D:/remo/qm/corgrids/jlcor/pre-REcor.nc"
simp=xr.open_dataset(sp)
tl="D:/remo/era_interim/25832_daymean_tp-mm_2000-2017.nc"
obsh=xr.open_dataset(tl).rename(Dict("tp"=>"pre"))
pyjl.doy(simh, simp, obsh; tosum=true)




z = simp
z = simp.assign_coords(time=z["time"])

grouped_obsh=obsh[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()
grouped_simh = simh[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()
grouped_simp = simp[k].mean("longitude").mean("latitude").groupby("time.dayofyear").sum()




pyjl.doy(simh, simp, obsh; tosum=true)
simp.variables|>collect
pyjl.doyplot(simh, simp, obsh; tosum=true)



#tl = "D:/remo/genRE/pr-dly.nc" #5min
tl="D:/remo/era_interim/projected_cropped/25832_daysum_tp-mm_2000-2017.nc"
fdoy(tl;tosum=true)
sp = "D:/remo/qm/corgrids/jlcor/pre-cor.nc" ##<-time is not a dimension
sp = "D:/remo/qm/corgrids/jlcor/rsds-cor.nc"
sp = "D:/remo/qm/corgrids/jlcor/pre-REcor.nc"
fdoy(sp)


z = xr.open_dataset("D:/remo/qm/corgrids/jlcor/rsds-cor.nc")
z.variables|>collect
z = z.assign_coords(time=z["time"])
#dimmns = Dict("longitude"=>"x","latitude"=>"y")
z = z.rename(dimmns)
dimmns = Dict("lon"=>"x","lat"=>"y")
z = z.rename(dimmns)

pyjl.fdoy(z)


z = xr.open_dataset("D:/remo/qm/corgrids/jlcor/rsds-cor.nc")

z.variables|>collect
k = z.keys()|>collect|>last
z = z.assign_coords(time=z["time"])
gry = z[k].mean("x").mean("y").groupby("time.dayofyear").mean()
vr = uppercase(first(k))
plt.figure()
plt.plot(gry, label="\$$vr _{sim,h}\$")
plt.title("$(uppercasefirst(k))")
plt.xlim(0, 365)
plt.gca().grid(alpha=0.3)
plt.legend()
plt.show()

#array method
z = xr.open_dataset("D:/remo/qm/corgrids/jlcor/pre-cor.nc")
grp = z.assign_coords(time=z["time"]).pre.mean("x").mean("y").groupby("time.dayofyear").mean()

#andersrum geht auch
grp = z.assign_coords(t=z.time).pre.mean("x").mean("y").groupby("t.dayofyear").mean()
plt.figure()
plt.plot(grp, label="\$array_{sim,h}\$")
plt.xlim(0, 365)
plt.gca().grid(alpha=0.3)
plt.legend()
plt.show()

z.assign_coords(t=z.time).pre.mean("x").mean("y").groupby("t.month").sum().plot()
z.assign_coords(t=z.time).pre.mean("x").mean("y").groupby("t.month").mean().plot().xlim(0,3)

ds = z.assign_coords(t=z.time)
#dsg = ds.groupby("time.month").sum().mean("month")
dsg = ds.groupby("time.month").mean()
dsg["pre"]
#ds_grouped.mean("x").mean("y").plot()


z.assign_coords(t=z.time)["pre"].isel(t=456).plot(mouseover=true)

out=z.assign_coords(t=z.time)
out = out.drop_vars("time")
out.to_netcdf("D:/remo/qm/corgrids/jlcor/pre-REcor2.nc")


rz = Raster("D:/remo/qm/corgrids/jlcor/pre-cor.nc")
rz = Raster("D:/remo/qm/corgrids/jlcor/pre-REcor2.nc")


#plot(rz[t=3])
findmax(rz[Dim{:t}=99])
rz.data[:,:,6]|>contourf
rz.data[:,:,99]|>contourf
rz.data[:,:,499]|>contourf

using Rasters
using Plots

# Load raster data
rz = Raster("D:/remo/qm/corgrids/jlcor/pre-REcor2.nc")

m = mean(rz, dims=3)
contourf(m)
#xmr = mean(skipmissing(rz), dims=(1, 2))
tm = lookup(rz,Ti)
@vv "CFTime"
using CFTime
# Compute mean along x and y axes
# Convert Raster to array
rz_array = convert(Array, rz)
# Skip missing values
rz_array_no_missing = skipmissing(rz_array)
# Calculate the mean along X and Y dimensions
rz_mean = mean(rz_array_no_missing, dims=(1, 2))

# Plot the mean values

# Group by day of year and calculate mean
grp = groupby(rz_mean, by=:dayofyear)
grp_mean = mean(grp)

# Plot
plot(grp_mean, label="array_sim_h")
xlims!(0, 365)
ylims!(minimum(grp_mean), maximum(grp_mean))  # Adjust ylims if needed
xlabel!("Day of Year")
ylabel!("Mean Value")
title!("Mean Value vs Day of Year")
grid(true)
legend(:topright)

@pyjl
lk="D:/remo/cordex/eobs/v28/rcm_eobs/rr_wa.nc"
lk="D:/remo/cordex/eobs/v28/rcm_eobs/qq_wa.nc"
Main.pyjl.fdoy(lk;tosum=false,txy=true)

pwd()
@edit qbb()

cd("D:/Wasim/regio/out/utm/utm_v6/station/")
@rgof


##run without compilation.
time julia --startup-file=no --color=yes --threads auto -q --compile=min --optimize=0 --compiled-modules=no "$jlpt/qsa.jl"


dx = dfroute(ofl=ofl)
rename!(dx,"sim"=>"Basin");

#qsa workflow
dx = rename(wa.dfroute(),"sim"=>"Basin")
dfs = getq();
dfs.Basin = parse.(Int64,dfs.Basin);
kd  = innerjoin(dfs, dx, on=:Basin);
pretty_table(kd,header=uppercasefirst.(names(kd));)

dfp("WLF")
dfp("qgko")

d = readall(r"WLF|qgko")
dm = mall(d)
selt(dm,r"7|W")|>baryrsum
selt(dm,r"7|W")|>dfl


cd("D:/Wasim/Tanalys/DEM/Input_V2/meteo/met0/")
v::Vector{String} = readdir()
files = v[(broadcast(x->occursin(r"txt",x),v))];
dds::Vector{DataFrame} = []
for file in files
    dtmp = nread(file)
    for nm in names(dtmp)
        newname = replace(nm, "\xc4" => "Ae", 
        "\xd6" => "Oe", "\xdc" => "Ue", "\xe4" => "ae",
        "\xf6" => "oe", "\xfc" => "ue", 
        "\xfc" => "ue", "/Main" => "", 
        #r"[_(].*" => "",
        r"[(].*" => "",
        "_" => " ",
        "\xdf" => "ss")
        rename!(dtmp, Dict(nm => newname))
    end
    push!(dds, dtmp)
end
hd(dds[1])
#dds[1]|>dfp
dst=("D:/Wasim/Tanalys/DEM/Input_V2/meteo/weekly/")

function toweek(pt::Union{String,DataFrame};fun=mean,agg=week)
    if pt isa String
        x = waread2(pt)
    else
        x = pt
    end
    
    df = copy(x)
    y = filter(x->!occursin("date",x), names(df))
    s = map(y -> Symbol(y),y)
    df[!, :date] .= agg.(df[!,:date]);
    df_agg = DataFrames.combine(groupby(df, :date), 
        y .=> fun .=> y);
    return(df_agg)
end

dds[1]|>toweek|>dfp

yearly_sum = combine(resample(dds[1], Year(1)), sum)


"""
running week and add confidence intervals
"""
function zzek(df::DataFrame;agg=week,fun=mean)
    # Make a copy of the input DataFrame
    dmean = copy(df)
    
    # Extract the columns of interest
    columns = filter(x -> !occursin("date", x), names(dmean))
    
    # Add a 'week' column based on the 'date' column
    dmean[!, :date] .= agg.(dmean[!,:date]);
    
    # Drop the 'date' column
    select!(dmean, Not(:date))

    dmean = DataFrames.combine(
        DataFrames.groupby(dmean, :date)) do group
        result = DataFrame(week = group.date)
        for col in columns
            values = group[!, col]
            result[!, col] .= fun(values)
        end
        return result
    end

    return dmean
end

function rweek(x::DataFrame; fun=mean)
    # Ensure the date column is of type Date
    df = copy(x)
    df.date = Date.(df.date)

    # Add a week column
    df.year = Dates.year.(df.date)
    df.week = Dates.week.(df.date)
    df.month = month.(df.date)
    
    # Group by week and calculate the mean
    aggdat = DataFrames.combine(groupby(
        df, [:year,:month, :week]), 
        names(df, Not([:date, :year, :week])) .=> fun;
        renamecols=false)

    # Recombine week to a date column
    #aggdat.date = Dates.firstdayofweek.(aggdat.week)
    aggdat.date = Date.(aggdat.year, aggdat.month)
    # Remove the week column
    select!(aggdat, Not([:week,:month, :year]))
    return aggdat
end

dmean = rweek(dds[1])
dfp(dmean)

cmk.dfp(selt(dmean,3))
cd(dst)
nm = replace(lowercase.(getnames(dds[1])), ".txt" => ".wk")
wawrite(dmean,nm)

getnames(dds)
for i in 2:length(dds)
    dmean = rweek(dds[i])
    nm = replace(lowercase.(getnames(dds[i])), ".txt" => ".wk")
    wawrite(dmean,nm)
end
#precisums:
dmean = rweek(dds[3])
nm = replace(lowercase.(getnames(dds[3])), ".txt" => ".wk")
wawrite(dmean,nm)
R"qt='c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe'"
24*7
#saale.raw0.txt
#inf='D:/Wasim/Tanalys/DEM/Input_V2/meteo/Pegel_1970_RAW.txt'

R"
inf='D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/ezg_sub.was'
onf='D:/Wasim/Tanalys/DEM/Input_V2/meteo/weekly/saale.wk'
system2(command = qt, args = paste(inf,onf,1,'all',24,2))
"
nread(lat())|>cmk.tsbar

mqt = @rget(qt)
#run((pipeline(`cmd.exe /c`,mqt)))
#qts=`c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe`
qts="c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe"
inf="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/ezg_sub.was"
onf="D:/Wasim/Tanalys/DEM/Input_V2/meteo/weekly/saale.wk"
#[<n1 headrows{def=5}>] [<n2 datacolumns{def=all}>] 
# [<timestep in h{def=1h}>]
# [<method 1(to mm/step) or 2(to m3/s) {def=1}>
#towsl(inf)|>cb
cp(inf,"qraw")
#nread(lat())|>cmk.tsbar
dfp(lat())

pwc()
inf="qraw"
run(`$qts $inf $onf 1 all 24 1`)

k = nread(lat())
# s=raw"C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/tst"
# oo = nread(s)
# oo|>cmk.tsbar
#ausprobieren
#and AL 1 ( modus = intern, 44.1 0.8 308 m3/s )


#fn=raw"L:\04-Phil\Geo1data\prj-Efre-Daten\GIS-Daten@Server K\Bodenschätzungskarten\Reichsbodenschätzung\Rhoen_BY_Teil\BR_Schichtinformation_utm.shp"
fn=raw"D:\Bodendaten\Bodenschätzungskarten\Reichsbodenschätzung\Rhoen_BY_Teil\BR_bodenschaetzung.shp"
fn=raw"D:\Bodendaten\Bodenschätzungskarten\Reichsbodenschätzung\Rhoen_BY_Teil\BR_Schichtinformation.shp"
import GeoDataFrames as gdf
a = gdf.read(fn)
b = unique(a,:FEUCHTE)
b = unique(a,:SBESCHR)
names(a)
view(b,2:3,3:ncol(b))
t = @rsubset a :OIDFESCH == "DEBYV07O8oWikzTz"
#plot(t.geometry)
cdof(fn)
#outname = "rbs.geojson"
outname = "schichtinfo.geojson"
gdf.write(outname,a; driver="GeoJSON")
op()
plot(a.geometry)

df=tsread("D:/Wasim/regio/out/rc200/eval-log.txt")
view(df,1,2)

rd = dropmissing(DataFrame(
    CSV.File("D:/Wasim/regio/out/rc200/x22/allroute",header=false,
    ntasks=1,
    delim=' ')),1)
    #comment="#",header=false)
sort!(rd,3;rev=true)


import GeoDataFrames as gdf
fn="D:/Wasim/regio/rcm200/v0/basins-v0.shp"
a = gdf.read(fn)
plot(a.geometry)

g="D:/Wasim/Tanalys/DEM/Input_V2/meteo/pegelstandorte_franken.txt"
g="D:/Wasim/Tanalys/DEM/Input_V2/meteo/Pegel_1970.txt"
peg = wa.stread(g)
plot!(peg.geometry,color=:black)

wa.stplot(g)

begin
    p = plot(a.geometry,color=:transparent)
    plot!(peg.geometry,color=:black)
    for (i, pt) in enumerate(peg.geometry)
        x = ArchGDAL.getx(peg.geometry[i], 0)
        y = ArchGDAL.gety(peg.geometry[i], 0)
        name = peg.name[i]
        annotate!(x, y, text(name, 8, :black, :bottom, :left))
    end
    plot!(p)
end
@vv ".crop("
# import Rasters as r
# pc = r.crop(peg.geometry;to=a.geometry)#fails

using RCall
@rimport terra as tr
ez=tr.vect("D:/Wasim/regio/rcm200/v0/basins-v0.shp")
R"""
wa.terra <- function(file, turn = F, crsval = 25832, encoding=encoding) {
  {
    if (!is.null(file))
    {
      st <- data.table::fread(
        input = file,
        skip = 0,
        header = T,
        nrows = 4,
        encoding = "UTF-8"
      )
    }
    else {
      cat("Error: please use filepath to read in")
    }
  }
  {
    if (turn == T)
    {
      cords <- data.frame(lon = as.numeric(t(st[3, 5:ncol(st)])),
                          lat = as.numeric(t(st[2, 5:ncol(st)])))
      ez = st[1, 5:ncol(st)]
      st = terra::vect(cords,crs=paste0("epsg:",crsval))
      st$rn = names(ez)#rownames(cords)
      st$ezg=as.numeric(ez)
      return(st)
    } else {
      cords <- data.frame(lon = as.numeric(t(st[2, 5:ncol(st)])),
                          lat = as.numeric(t(st[3, 5:ncol(st)])))
      ez = st[1, 5:ncol(st)]
      st = terra::vect(cords,crs=paste0("epsg:",crsval))
      st$rn = names(ez) #rownames(cords)
      st$ezg=as.numeric(ez)
      return(st)
    }
  }
}
"""
rpe = R"""wa.terra('D:/Wasim/Tanalys/DEM/Input_V2/meteo/Pegel_1970.txt')"""
tcr = tr.mask(rpe,ez)
cd("D:/Wasim/regio/rcm200/v0/")
tr.writeVector(tcr,"pegel_masked.shp",overwrite=true)

lat()
peg = gdf.read("pegel_masked.shp")
sort!(peg, :ezg;rev=true)
rename!(peg,:rn=>:name)
begin
    p = plot(a.geometry,color=:transparent)
    plot!(peg.geometry,color=:black)
    for (i, pt) in enumerate(peg.geometry)
        x = ArchGDAL.getx(peg.geometry[i], 0)
        y = ArchGDAL.gety(peg.geometry[i], 0)
        name = peg.name[i]
        annotate!(x, y, text(name, 8, :black, :bottom, :left))
    end
    plot!(p)
end

#lnk to shape
lnk = tr.as_polygons(tr.rast("rcm.lnk"))
tr.set_crs(lnk, "EPSG:25832" ) 
tr.writeVector(lnk,"riverlinks.shp",overwrite=true)

@edit vgctl("rcm200/v0")

R"""
lfp <- function(x) {
  require(leaflet)
  require(leaflet.opacity)
  require(leaflet.extras)
  require(leafem)
  
  if (inherits(x, "SpatRaster")) {
    nf = round(raster::ncell(x)/100, 0)
    pl = hcl.colors(n = nf, palette = "Dark 3")
    vl = pretty(as.numeric(stats::quantile(values(x),
                                           na.rm = TRUE,
                                           probs = seq(0, 1, length.out = length(pl)))))
    
    leaflet() %>%
      addTiles() %>%
      addProviderTiles(providers$OpenStreetMap.DE, group = c("OSM")) %>%  
      addProviderTiles("Stamen.TonerHybrid", group = c("Stamen")) %>% 
      addProviderTiles(providers$OpenTopoMap, group = c("Topography")) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = c("Esri.WorldImagery")) %>%
      addRasterImage(
        x = x,
        opacity = 0.8,
        colors = pl,
        layerId = names(x),
        group = names(x)
      ) %>%
      leaflet.opacity::addOpacitySlider(layerId = names(x)) %>%
      leafem::addMouseCoordinates() %>%
      addWMSTiles(
        "http://www.lfu.bayern.de/gdi/wms/geologie/dgk25?",
        layers = c("geoleinheit_dgk25", 'strukturln_dgk25'),
        options = WMSTileOptions(format = "image/png", transparent = TRUE),
        group = 'GK25'
      ) %>%
      addWMSLegend(
        uri = paste0("http://www.lfu.bayern.de/gdi/legende/geologie/dgk25/geoleinheit_dgk25.png")
      ) %>%
      addLayersControl(
        baseGroups = c('Topography', 'OSM', 'Stamen', 'Esri.WorldImagery'),
        overlayGroups = c('GK25', names(x)),
        options = layersControlOptions(collapsed = TRUE),
        position = "topleft"
      ) %>%
      addResetMapButton() %>%
      addSearchOSM() %>%
      addMeasure(
        position = "bottomleft",
        primaryLengthUnit = "meters",
        primaryAreaUnit = "sqmeters",
        activeColor = "#3D535D",
        completedColor = "#7D4479"
      )
  } else 
    if (inherits(x, "SpatVector")) {
    # Process the SpatialVectors data here
    dat = terra::as.data.frame(
      terra::project(x,"epsg:4326"),geom="XY")
    leaflet() %>%
      addTiles() %>%
      addProviderTiles(providers$OpenStreetMap.DE, group = c("OSM")) %>%  
      addProviderTiles("Stamen.TonerHybrid", group = c("Stamen")) %>% 
      addProviderTiles(providers$OpenTopoMap, group = c("Topography")) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = c("Esri.WorldImagery")) %>%
      addMarkers(lng=dat$x,lat=dat$y) %>% 
      addPopups(lng=dat$x,lat=dat$y,popup=dat$rn) %>% 
      #leafem::addFeatures() %>% 
      leafem::addMouseCoordinates() %>%
      addWMSTiles(
        "http://www.lfu.bayern.de/gdi/wms/geologie/dgk25?",
        layers = c("geoleinheit_dgk25", 'strukturln_dgk25'),
        options = WMSTileOptions(format = "image/png", transparent = TRUE),
        group = 'GK25'
      ) %>%
      addWMSLegend(
        uri = paste0("http://www.lfu.bayern.de/gdi/legende/geologie/dgk25/geoleinheit_dgk25.png")
      ) %>%
      addLayersControl(
        baseGroups = c('Topography', 'OSM', 'Stamen', 'Esri.WorldImagery'),
        overlayGroups = c('GK25', names(x)),
        options = layersControlOptions(collapsed = TRUE),
        position = "topleft"
      ) %>%
      addResetMapButton() %>%
      addSearchOSM() %>%
      addMeasure(
        position = "bottomleft",
        primaryLengthUnit = "meters",
        primaryAreaUnit = "sqmeters",
        activeColor = "#3D535D",
        completedColor = "#7D4479"
      )
    
    # Return the modified leaflet map
  } else {
    stop("Input must be either a SpatRaster or SpatialVectors object.")
  }
}
"""
@rput ez
@rput rpe
R"""lfp(rpe)"""
@rimport mapview as rmv
@rcall
pwd()
@doc rr.nctc("slp")
ra = rr.rastr(r"sl")
tr.summary(ra)
tr.plet(ra)

rmv.mapview(ra)

fn = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/Pegel_1970.txt"
#peg = tsread(fn)
peg = fread(fn)
obs = fread("D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt")

xm = mall([peg,obs])
selt(xm,r"Wolf")|>dfp

psp = fread("D:/Wasim/Tanalys/DEM/Input_V2/meteo/spec_discharge1970_2017.txt")
@vv "qtospend"
qts="c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe"
#[<n1 headrows{def=5}>] [<n2 datacolumns{def=all}>] 
# [<timestep in h{def=1h}>]
# [<method 1(to mm/step) or 2(to m3/s) {def=1}>
#towsl(inf)|>cb
inf = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/Pegel_1970.txt"
inf = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/peg_blank.wa"
onf = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/spec.wa"
# run(`$qts $inf $onf 5 all 24 1`)
# run(`$qts $inf $onf 5 16 24 1`)
cd("D:/Wasim/Tanalys/DEM/Input_V2/meteo/pegel")
#cp -v ../Pegel_1970.txt . 
#hier from sosubs but longer names.
#    #qts=$(cmd.exe /c "/mnt/c/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe");
#qts="/mnt/c/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe";

# splitcols () 
# { 
#     [[ -z "$1" ]] && { 
#         printf "Parameter is empty ...\nusage: <thiscmd> <time-series>\n";
#         return 1
#     };
#     b="$1";
#     END=$(perl -ae 'print ($#F) if $.==1' $b);
#     for ((i=4; i<=END; i++))
#     do
#         q=$(perl -ae 'print $F['$i'] if $.==1' $b);
#         perl -splae '$_=join q{ },@F[0..3,'$i'];s/\h+/\t/g;' $b > ${q:0:23}.wa;
#         GR='\033[0;32m';
#         printf "%s ${GR}${q:0:23}.wa written! \n${GR}";
#     done
# }

files = glob(".wa")
for f in files
    ofl=replace(f,".wa"=>".spec")
    run(`$qts $f $ofl 5 all 24 1`)
end

xm = mall(readall(glob("wolf")))
dfp(xm)
qplot(xm)
ftplin(xm)
dfm(xm;log=true)

pin = gdf.read("D:/Wasim/regio/rcm200/v0/pegel_masked.shp")
al = readdir()
want = pin.rn
#filter by readdir by want with suffix .wa
al = filter(x->occursin(r".spec",x),al)
se = [i[1:6] for i in al]
aa = [i[1:4] for i in want]
findalls(aa)
flist = []
for i in al
    for j in aa
        if occursin(j,i)
            println("$i found")
            push!(flist,i)
        end
    end
end
unique!(flist) # rm dubs from list
dfs = readall(flist)
dfs = mall(dfs) #147 cols
Main.pyjl.wawrite2(dfs,"newspec")
oo=tsread("newspec")
oo=dfread("newspec")
baryrsum(oo)
dfm(oo)

ns = @rsubset oo :date >= Date("2012-12-31")
baryrsum(ns)

##noch die header...
fn
hdr = CSV.read(fn,DataFrame;limit=5)
sk = hdr[!,1:4]
na = hdr[!,5:end]
#sort columns by name
na = na[!,sort(names(na))]
#subset by names of oo
nn = join([i[1:4] for i in names(oo)],"|")
na = select(na,Cols(Regex(nn))) #now only 144 cols
hdr = hcat(sk,na)
hdr = hdr[1:4,:]
@edit findindf(hdr, "-9999.0")
@subset hdr "-9999.0"
cols_with_neg9999 = [col for col in names(hdr) if any(hdr[4, col] .== -9999.0)]

dout = copy(dfs)
dout.YY = map(x ->year(x),dout.date);
dout.MM = map(x ->month(x),dout.date);
dout.DD = map(x ->day(x),dout.date);
dout[!, "HH"] .= 24;
dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))];
v = map(names,[dout,hdr])
map(length,v)
setdiff(v...)    #find differencies
for x in names(hdr)
    if occursin(r"_",x)    
    newname=replace(x,"_"=>"-")
    rename!(hdr,Dict(x=>newname))
    end
end
#subset by name
dout = select(dout,Regex(join(names(hdr),"|")))
hdr = select(hdr,Regex(join(names(dout),"|")))
dout = vcat(hdr,dout)
#47 bad cols.
cols_with_neg9999 = [col for col in names(dout) if any(dout[4, col] .== -9999.0)]

xo = select(dout,Not(Cols(cols_with_neg9999)))
misval=-9999.0
CSV.write("specdis_sel.txt",xo, 
    transform = (col, val) -> something(val, misval), delim="\t")  


#get areas: 
#gd = gdf.read("D:/Wasim/regio/rcm200/v0/basins-v0.shp")
sf="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
R"""
gd = wa.terra($sf)
ezv = gd$ezg
"""
@rget ezv
gd = wa.stread(sf)
gd.name
b = DataFrame(name=string.(gd.name),ez=ezv)
#arnstein testField
arn = dfr("Arnstein.wa")|>dropmissing
@rsubset b :name == "Arnstein"
arn.spec = (arn.Arnstein ./ (328.2)*1000)
obs
asp = selt(obs,r"Arnstein")

q = mall([arn,asp])
rename!(q,1=>"raw")
dfl(q[!,2:4]) #WHY? doublecheck mit qtospend



function calculate_runoff(files::Vector{String}, catchment_areas::Vector{Float64})
        # Überprüfen, ob die Anzahl der Dateien und die Länge des Einzugsgebietsvektors übereinstimmen
        if length(files) != length(catchment_areas)
            error("Die Anzahl der Dateien und die Länge des Einzugsgebietsvektors müssen übereinstimmen.")
        end
    
        # Initialisieren Sie den Vektor für die Abflussspende
        runoff_per_catchment = Vector{Float64}(undef, length(files))
    
        # Durchlaufen Sie jede Datei und berechnen Sie die Abflussspende
        for (i, file) in enumerate(files)
            # Lesen Sie die Abflussdaten aus der Datei
            discharge_data = readdlm(file, Float64)
    
            # Berechnen Sie die Abflussspende
            runoff_per_catchment[i] = sum(discharge_data) / catchment_areas[i]
        end
    
        return runoff_per_catchment
end
    
using RCall
vgr("Flowscreen")
#R""".libPaths(new="C:/Users/chs72fw/AppData/Local/R/win-library/4.2")"""
fc=rimport("FlowScreen")
lk="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/Frnkische-24409003-Wolfsmuenster.was"
wts = nread(lk)
wts = @rsubset wts :date >= Date("1950-01-01")
rename!(wts,1=>"Flow")
wts.Flow = parse.(Float64,wts.Flow)
wa.hydromon(wts)
wa.hydro(wts)
wa.dfm(wts;fun=monmean)
wa.mbx(wts) #<-eher da hier rein.

lk=
cd("D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_15Min_Teil1/")
glob("olfs")
cdu()
rglob("olfs")
fn=last(rglob("olfs"))
#fn="D:\\Wasim\\Pegeldaten\\2018_08_08_Bodenwasserhausaltsmodellierung_Franken\\qm_hourly\\qm1h.txt"
fn="D:\\Wasim\\Pegeldaten\\2018_08_08_Bodenwasserhausaltsmodellierung_Franken\\qm_hourly\\qm_saale_h.txt"
@edit nread(fn)
#d2 = nread(fn;skipto=20) ##VERY big.
d2 = dfr(fn)
names(d2)|>println
#selt(d2,r"saa")
sh = selt(d2,4)
rename!(sh,1=>"wlf")
#ERROR: ArgumentError: cannot parse String31("") as Float64
#sh.wlf = replace(sh.wlf,""=>missing)
#sh.wlf = replace.(sh.wlf, "" => missing)
#sh.wlf = parse.(Float64,sh.wlf)
#sh.wlf = tryparse.(Float64,sh.wlf)
#sh.wlf = replace.(x->isnothing(x)=>missing,sh.wlf)
#dropmissing!(sh)
wa.hydro(sh)
dropmissing!(sh)
wa.mbx(sh)
title!("WLF from qm_saale_h.txt")
savefig("wlf_qm_saale_h.png")

ser = fc.create_ts(wts.Flow)
@rput wts
#@doc rcopy(wts,"ts")

R"""
"""

#qtosp $a saale_spec4 5 all 24 2

#qtosp $a saale_spec2 5 all 24 1
#qtosp $a saale_spec3 
###check if spec was right:
x="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/qm_hourly/saale_spec2"
x="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/qm_hourly/saale_spec3"
x="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/qm_hourly/saale_spec4"
#qtosp $a.new saale_spec5 5 all 24 1
x="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/qm_hourly/saale_spec5"
spn=dfr(x)
baryrsum(spn)

##spn stimmt schon, es ist stündliche Gebietsspende

baryrsum(spo)
xk="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
spo = dfr(xk)
spo=selt(spo,r"Wolf")
spn=selt(spn,r"Wolf")
xm=mall(spo,spn)
qplot(xm)










?FlowScreen::bf_oneparam()
rebf <- bf_oneparam(wts$Flow, k=0.9)
plot(wts$Date, wts$Flow, xlab="", ylab="Q (m3/s)", type="l")
points(wts$Date, rebf, type="l", col="blue")
cat("mean baseflow",mean(rebf),"Q (m3/s)")
summary(rebf)
?bf.stats
wts[,hyear:=year(wts$Date)]
res <-bf.stats(wts, by = "hyear")
View(bf.stats)
res2 <- screen.metric(res[,2], "m3/s",text = names(res)[2])
#z=read.flows(lk)
res2 <- screen.metric(res$BFVol, "m3/s",text =
                        "BFVol - Baseflow volume for the summary period, in km3")


screen.metric(res$MinBF, "m3/s",text =
                "Min daily baseflow for the summary period, in m3/s")


                # Erstellen des DataFrames mit den gegebenen Daten
df = DataFrame(
                    Variable = ["Niederschlag", "Versickerung", "Transpiration", "Interzeption", "Bodenverdunstung", "Tatsächliche Evapotranspiration"],
                    min = [770, 200, 270, 150, 30, 480],
                    max = [800, 300, 300, 280, 50, 680]
                )
                
                # Erstellen des Boxplots
@df df boxplot(:Variable, [:min, :max], title="Niederschlagsdaten", ylabel="mm", xlabel="Variable", legend=:topright)                

# df = DataFrame(
#     Variable = repeat(["Niederschlag", "Versickerung", "Transpiration", "Interzeption", "Bodenverdunstung", "Tatsächliche Evapotranspiration"], inner=2),
#     Wert = [770, 800, 200, 300, 270, 300, 150, 280, 30, 50, 480, 680],
#     Typ = repeat(["min", "max"], 6)
# )

df = DataFrame(
                    Variable = ["Niederschlag", "Versickerung", "Transpiration", "Interzeption", "Bodenverdunstung", "Tatsächliche Evapotranspiration"],
                    min = [770, 200, 270, 150, 30, 480],
                    max = [800, 300, 300, 280, 50, 680]
                )
df.avg .= ((df.min .+ df.max) ./ 2)

sizes = df.avg
labels = df.Variable
pie(df.avg,labels = labels, legend=false, aspect_ratio=1)
# Add labels with percentages using annotate!
angles = cumsum(sizes) .- sizes ./ 2
for (label, size, angle) in zip(labels, sizes, angles)
    percentage = string.(round(size / sum(sizes) * 100, digits=1))
    x = cosd(angle)
    y = sind(angle)
    annotate!(x, y, 
    text("$label\n$percentage %",8,
    :center;
    family="Computer Modern"))        
end

title!("Pie Chart with Labels and Percentages")

percentage=[]
for size in df.avg
    push!(percentage,
    (df.Variable[size]*(string(
        round(size / sum(df.avg) * 100, 
        digits=2)))*" %"))
end

#df.prc = percentage
df.prc = df.avg / sum(df.avg) * 100

##annotated piechart.
using Makie
begin
    data   = df.avg    
    labels = hcat.(df.Variable .* "\n" , string.(round.(df.prc, digits=2)) .* " %")
    labels = join.(labels," ")

    wcolors = Makie.wong_colors() #[:yellow, :orange, :red, :blue, :purple, :green]
    unique_values = unique(data)
        # Map unique values to colors
    value_to_color = Dict(unique_values[i] => 
        wcolors[i % length(wcolors) + 1] for i in 1:length(unique_values))
        # Map values in a to colors
    scolors = [value_to_color[value] for value in data]
    
    r0, r1 = 1, 4

    fig, ax, plt = Makie.pie(data, color=scolors, 
        radius=r1, inner_radius=r0,
            strokecolor=:white, strokewidth=5, 
            axis=(autolimitaspect = 1, )
    )

    rot(θ) = 90<θ<270 ? θ*π/180 - π : θ*π/180
    radius(θ, r) = 90<θ<270 ? 1.25*r : 1.05*r   # adjust function of label sizes

    θ = (cumsum(data) - data/2) .* 360/sum(data)
    scθ = sincosd.(θ)

    for (li, θi, sci) in zip(labels, θ, scθ)
        ri = radius(θi, r1)
        Makie.text!(ax, li, 
        position=ri.*(sci[2], sci[1]), 
        rotation=rot(θi), fontsize=10, space=:data)
    end

    #Makie.ylims!(ax, -5.5, 5.5)
    hidedecorations!(ax)            # hides ticks, grid and lables
    hidespines!(ax)    
    fig
end


using Makie
begin
    #data = sort(data, rev=false) # Sort data in descending order
    data = df.avg
    data = sort(data, rev=false) # Sort data in descending order
    prc = sort(df.prc, rev=false) 
    labels = df.Variable[sortperm(data, rev=false)] # Sort labels according to sorted data
    #labels = labels[sortperm(data, rev=false)] # Sort labels according to sorted data
       
    unique_values = unique(data)
    # Map unique values to colors
    value_to_color = Dict(unique_values[i] => 
        wcolors[i % length(Makie.wong_colors())] for i in 1:length(unique_values))
        # Map values in a to colors
    scolors = [value_to_color[value] for value in data]
    

    nc = length(scolors)
    nδ = 1/nc
    
    r0, r1 = 1, 4
    f, ax, plt = Makie.pie(data,
                    #offset = 0.5, #turns plot
                    color = scolors,
                    radius = r1,
                    inner_radius = r0,
                    strokecolor = :white,
                    strokewidth = 2.5,
                    axis = (autolimitaspect = 1, )
                    );
    cbar = Makie.Colorbar(f[1,2], colormap = cgrad(scolors, categorical = true))
    cbar.ticks = (range(0+nδ/2, 1-nδ/2, nc), labels)
    hidexdecorations!(ax)
    hideydecorations!(ax)
    θ = (cumsum(data) - data/2) .* 360/sum(data)
    rot(θ) = 90<θ<270 ? θ*π/180 - π : θ*π/180
    # adjust function of label sizes
    radius(θ, r) = 90<θ<270 ? .75*r : .50*r    #1.25*r : 1.05*r   
    scθ = sincosd.(θ)
    for (li, θi, sci) in zip(string.(round.(prc,digits=1)), θ, scθ)
        ri = radius(θi, r1)
        Makie.text!(ax, li, 
        position=ri.*(sci[2], sci[1]), 
        rotation=rot(θi), fontsize=14, space=:data)
    end
    f    
end



#gkplot
fn="D:/Wasim/Tanalys/DEM/brend_fab/x0/br.gk1"
fn="D:/Wasim/regio/rcm200/v11/rcm.gkm"
agheat(fn)

fn="D:/Wasim/regio/control/rcm-c5_loc7.ctl"
sd=wa.fsoil(fn)
sd=wa.read_soildata_4(fn)
dx = split.(grep("thick",sd))
dx = [i[3:5] for i in dx]
dx = DataFrame(dx,:auto)
#parse.(Float64,dx)
dfloat!(dx)
ou = colsums(dx)
cuo = cumsum.(eachcol(dx[!,Not(Cols(r"date|month|year"))]))
lines(ou)
lines([first(i) for i in cuo])
lines!([last(i) for i in cuo])
title!("Soil thickness")

scatter([first(i) for i in cuo],
[last(i) for i in cuo])

wp=nread("D:/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/wern_2012")
wp = @rsubset wp :date >= Date("2013-01-01")
hydromon(wp)
wa.hydromon(wp,col=2)
wa.hydromon(wp,col=3)
wa.hydro(wp,col=3)
wa.hydro(wp,col=2)
cd("D:/Wasim/regio/out/lowres/c5/loc8/")
ds=pyjl.allkge()
sort!(ds,1)

dfs=readall(r"qoutjl")
dyr  = map(byear,dfs)
plot_grouped_metrics(dyr)
plot_grouped_metrics(dyr;col="kge",threshold=.4)
pxm(dyr)


#sting partial functions
filter(in("D"),pwd())
occursin("D",pwd())
in("Abracadabra")('a')
vv = ['a'+i:'z' for i in 0:4]
map(filter(in("Abracadabra")), vv)


cd("D:/Wasim/regio/out/lowres/c5/loc9")
outd = readall(r"qoutjl$")
plts = [ wa.hyeval(m;fun=mean) for m in outd];
describe(plts)
#all in one layout with shared axis labels:
plot(plts..., layout=(4,3), size=(1200,1600) .* 1.0, 
    sharex=true, sharey=true, xlabel="", ylabel="")


dyr=byear.(outd)
wa.plot_grouped_metrics(dyr)
wa.pxm(dyr;threshold=.6)#,yaxis=:log
yaxis!(:log)

raw"D:\Wasim\regio\out\lowres\v4\f0"|>cd
outd = readall(r"qoutjl$")
plts = [ wa.hyeval(m;fun=mean) for m in outd];
describe(plts)
#all in one layout with shared axis labels:
plot(plts..., layout=(4,3), size=(1200,1600) .* 1.0, 
    sharex=true, sharey=true, xlabel="", ylabel="")
Plots.savefig("mon_score_v4.svg")

dyr=byear.(outd)
wa.plot_grouped_metrics(dyr)
wa.pxm(dyr;threshold=.0)


cd("D:/Wasim/regio/out/rc200/x22/cl4/")
outd = readall(r"qoutjl$")
plts = [ wa.hyeval(m;fun=mean) for m in outd];
describe(plts) #9
#all in one layout with shared axis labels:
plot(plts..., layout=(3,3), size=(1200,1600) .* 1.0, 
    sharex=true, sharey=true, xlabel="", ylabel="")
Plots.savefig("mon_score_x22-cl4.svg")
    

###wern###########

araw="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/Wern-24382304-Arnstein.was"
ad = dfr(araw)
dfm(ad)
hydromon(ad)
ads = @rsubset ad year.(:date) > 1990
dfm(ads;mode=:bar,
fun=yrsum,ann=false)

dfp(ads)
fn="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/Wern-24381503-Ettleben.was"
edf = waread2(fn)|>dropmissing
dfp(edf)
fn="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_TagMit/wasim/Wern-24385007-Sachsenheim.was"
sadf = waread2(fn)|>dropmissing
dfp(sadf)

wm = mall([ad,edf,sadf])
dfp(wm)
hydromon(wm)
cdof(fn)
#another renamer.
rename!(wm, Dict(i => name for (i, name) in zip(1:3, ["Arnstein","Ettleben","Sachsenheim"])))
wm[!,[2,1,3,4]]     #reorder
wawrite(wm[!,[2,1,3,4]],"wern.wa")
@vv "qtospend"
qts="c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe"
inf="wern.wa"
#take headers from:
#D:\Wasim\Tanalys\DEM\Input_V2\meteo\ts-0\wern_2012
npp(inf)
onf="wern.spec"
#[<n1 headrows{def=5}>] [<n2 datacolumns{def=all}>] 
# [<timestep in h{def=1h}>]
# [<method 1(to mm/step) or 2(to m3/s) {def=1}>
#towsl(inf)|>cb
cp(inf,"wern-raw.wa")
run(`$qts $inf $onf 5 all 24 1`)
npp(lat()) ##check an cleanup
nread(lat())|>cmk.tsbar
dsp = nread(lat())
dsp = @rsubset dsp year.(:date) > 2012
hydro(dsp)
hydro(dsp,col=2)
hydro(dsp,col=3)
#dfp(lat())
pwc()

#compaire with specdis_kmu
xk="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
obs = nread(xk)
eto = selt(obs,r"Ettleben")
tm = mall(eto,selt(dsp,r"Ettleben"))
tm = reorder_df(tm)

qplot(tm)
dfl(tm)
cmk.cloudplot2(tm)
#wa.corrbar(select(tm,1),select(tm,3))
mbx(tm)
dfm(tm)

dst="D:/Wasim/Tanalys/DEM/Input_V2/meteo/wern.spec"
cp("wern.spec",dst)
#p024 wern.spec

using RCall
ad = rimport("AnomalyDetection")
@doc ad.AnomalyDetectionTs
pt=raw"D:\Wasim\regio\rcm\pestpp\v1\dsub.wa"
df = dfr(pt)
df=df[,c(1,3)]

as = selt(df,r"Arnstein")
# ad.AnomalyDetectionTs(as, max_anoms=0.02,longterm = true,
#                    verbose = true,
#                     direction="both", plot=true)
as = select(as, :date, :)
@rput as
#require(lubridate)
#as$date = ymd(as$date)
R"""
require(AnomalyDetection)
#as$date = as.POSIXct(strptime(as$date,format="%Y-%m-%d"))
as$date = as.POSIXct(as$date,format="%Y-%m-%d")
AnomalyDetectionTs(as,max_anoms=0.02,longterm = F,
verbose = T,
 direction='both', plot=TRUE)
"""


#CSV.read("D:/Wasim/regio/rcm/sld",DataFrame; delim=" ",header=false)

using DataFrames

function uct(file)
    counts = Dict()
    open(file) do f
        for line in eachline(f)
            words = split(line, "\t", keepempty=false)
            for word in words
                if word != "-9999"
                    counts[word] = get(counts, word, 0) + 1
                end
            end
        end
    end
    df = DataFrame(value = tryparse.(Float64, #cause of header:nrows etc.
        collect(keys(counts))), 
        count = 
        collect(values(counts))) #|>dropmissing
        #tryparse.(Int64,)
    df = df[.!isnothing.(df.value), :]
    #transform!(df, :value => (x -> coalesce.(x, missing)) => :value) #add columns
    sort!(df, :count, rev=true)  # sort by count in descending order
    return df
end

fn = "D:/Wasim/regio/rcm/rcm.wit"
af = uct(fn)
#transform!(af, :value => (x -> coalesce.(x, nothing)) => :value)
#plot(af.count,af.value, seriestype=:bar, legend=false, 
#    xlabel="Soil type", ylabel="Count", title="Soil type distribution")
#r = Raster(fn;missingval=-9999.0)
#r = readras2(fn)
r = Raster(fn)
r.data
#lookup(r,2)
z = r.data
missval = -9999
reverse!(z,dims=1)
reverse!(z,dims=2)
replace!(z, missval=>missing)

@edit rst.stats(fn;missingval=-9999.0)

fn="D:/remo/qm/corgrids/jlcor/pre-REcor.nc"
#fn="D:/remo/qm/corgrids/jlcor/pre-REcor2.nc" #<--corrput.
#fn="D:/remo/qm/corgrids/jlcor/pre-cor.nc" #<--corrput.?
fn="D:/remo/qm/corgrids/pre/pre-cor_extremes.nc"
pr=Raster(fn,key=:pre)
@vv "lookup"
@vv "Near("
@vv "CFTime"
name(pr)
lookup(pr,:Ti)
plot(pr)

tm = lookup(pr,Ti)
tm = parse.(DateTime,string.(tm))

#tm = map(x->tryparse(DateTime,string(x)),tm)
pr[1,3,2]
pr[:,:,2]|>contour
pr[:,:,end-2]|>contour
pr[:,:,end-2]|>cmk.mkrheat
pr[:,:,end-2]|>plot

#assing tm to Ti
# pr[:Ti] = tm error
#new = Raster(pr.data, dims=(lookup(pr,X), lookup(pr,Y), tm))
ln, lt = X(lookup(pr,X)), Y(lookup(pr,Y))
#ti = Ti(parent(tm))
#ti = Ti(DateTime(1997):Month(1):DateTime(2098))
ras = Raster(pr.data, dims=(ln,lt,tm)); #<-semicolon.Display error
#ras[X(Rasters.Near(10.0)),Y(Rasters.Near(49.80))]|>plot
#pr[X(Rasters.Near(10.0)),Y(Rasters.Near(49.80))]|>plot
s = pr[X(Rasters.Near(10.0)),Y(Rasters.Near(49.80))]
using Makie
#Rasters.rplot(s)
Rasters.rplot(pr)
Rasters.rplot(pr[X(Near(10.0))])

sdx = rst.ncdf(pr)
vio(sdx)
aa = monsum(sdx)
wa.tbx(sdx)
vibx(sdx)
ovio(sdx)

df = rst.ncdf(pr)
hydro(df)
dfm(df,fun=yrsum,ann=false,mode=:bar,log=false)
hydromon(df)
mbx(df)
vio(df)
findmax(df.pre)
using RCall
ad = rimport("AnomalyDetection")
#@doc ad.AnomalyDetectionTs
df = reorder_df(df,true)
@vv "weekl"
# ad.AnomalyDetectionTs(df,  #<-DateTime needed.
#     max_anoms=0.02,
#     longterm = true,
#     verbose = true, direction="both", 
#     plot=true)

ot = ad.AnomalyDetectionVec(select(df,2), 
    # only_last =,
    max_anoms=0.02, 
    direction="pos", 
    threshold="p95",    #"p99",
    period=365,
    plot=true)

#######task: find occurence of rcm rastercells in ug############
import GeoDataFrames as gdf
pr[Ti=1]
pt = raw"D:/Wasim/regio/rcm200/v13/cmtv13.shp"
#cdof(pt)
ug = gdf.read(pt)
#ras[Ti=At(DateTime(2001))]
#plot!(ug.geometry, color=:transparent, legend=false)
#plot(g.geometry, fillcolor=false)
ug = wa.reverse_coords(ug) ##err
@rimport terra as tr
g = tr.vect(pt)
g = tr.project(g, "epsg:4326")
tr.writeVector(g, "D:/Wasim/regio/rcm200/v13/v13_4326.json")
#gd = tr.crds(g)
#@rget gd
#gp = gdf.DataFrame(nm=g.rcm,geometry=pol)
pt = "D:/Wasim/regio/rcm200/v13/v13_4326.json"
gp = gdf.read(pt)

fn="D:/remo/qm/corgrids/pre/pre-cor_extremes.nc"
pr = Raster(fn,key=:pre)
#r = pr[Ti=Near(DateTime(2024,1,3))]
r = pr[Ti=Near(DateTimeNoLeap(2024,1,3))]
using Plots
plot(r)
plot!(gp.geometry, fillcolor=false)
#plot(gp.geometry, fillcolor=false)
r = pr[mean(Ti)]
#rx = mask_trim(r, gp;pad=1)
rx = mask_trim(pr, gp;pad=1)
#rx = Rasters.mask(r; with=gp,boundary=:touches)
#mean(A, dims=X) # Ti if time were available would also be possible
rmean = mean(rx, dims=Ti)
plot(rmean,title="Mittlerer Niederschlag REMO",xlabel="",ylabel="")
plot!(gp.geometry, fillcolor=false)
#Rasters.
write("D:/remo/qm/corgrids/pre/p_mean_crop.nc",rmean)

##subset by time with [..]
ot = view(rx, Ti(DateTimeNoLeap(1980,1,1,12)..DateTimeNoLeap(2098,1,1,12)))
ot = sum(ot, dims=Ti)
yrs = 2098-1980 #<-years
#ot[Where()]
plot(ot ./ yrs,title="Jährlicher Niederschlag REMO",xlabel="",ylabel="";c=:matter)

##get outline only.
ex = reduce(gdf.union,gp.geometry)
coords = GeoInterface.coordinates(ex)
extrema.(coords)
flat_coords = vcat(coords...)
min_x = minimum(x -> x[1], flat_coords)
max_x = maximum(x -> x[1], flat_coords)
min_y = minimum(x -> x[2], flat_coords)
max_y = maximum(x -> x[2], flat_coords)
#ox = extrema(flat_coords) #overall_extrema
#cu = Rasters.trim(Rasters.crop(pr; to=ex, touches=true);pad=2)
#cu = Rasters.trim(Rasters.mask(pr; with=ex, touches=true);pad=2)
#cu = sum(pr[X(ox[1][1]..ox[2][1]),Y(ox[1][2]..ox[2][2])],dims=Ti)
#cu = sum(pr[X(min_x..max_x),Y(min_y..max_y)],dims=Ti)
yrs = unique(year.(parse.(DateTime,string.(lookup(pr,Ti)))))
cu = sum(pr[X(9.2..10.8),Y(50.0..50.8)],dims=Ti)
cy = cu ./ length(yrs)
contourf(cy,title="Jährliche Niederschlagssummen",
    xlabel="",ylabel="";c=:matter)
plot!(ex, fillcolor=false)
#writeit
write("D:/remo/qm/corgrids/pre/p_yrmean_crop.nc",cy)
towsl(dirname("D:/remo/qm/corgrids/pre/p_yrmean_crop.nc"))|>cb


mean(cu, dims=Ti) |> plot
plot!(ex, fillcolor=false)

ot = view(cu, Ti(DateTimeNoLeap(1980,1,1,12)..DateTimeNoLeap(2098,1,1,12)))
ot = sum(ot, dims=Ti)
yrs = 2098-1980 #<-years
contourf(ot ./ yrs,title="Jährlicher Niederschlag REMO",xlabel="",ylabel="";c=:matter)
plot!(ex, fillcolor=false)
plot!(gp.geometry, fillcolor=false)


#sc = [DateTimeNoLeap(1980,1,1,12),DateTimeNoLeap(2098,1,1,12)]
#make a range of dates
using CFTime
sc = DateTimeNoLeap(1980,1,1,12):Day(1):DateTimeNoLeap(2098,1,1,12)
ag = Rasters.aggregate(Rasters.Center(), rx, sc)
#rx[X=3,Y=3]|>plot
#extend(rx[X=3,Y=3])
#ts = pr[X(10.0..10.5), Ti(DateTimeNoLeap(2001, 5)..DateTimeNoLeap(2001, 6))]
#plot(ts)
rmean[X(Near(10.0)),Ti(1)] |> plot
pr[X(Near(10.0)),Y(Near(50.33))] |> plot

tm = lookup(pr,Ti)
tm = parse.(DateTime,string.(tm))

#######task: find occurence of rastercells in ug############
gp = gdf.read("D:/Wasim/regio/rcm200/v13/cmtv13.shp")
px = first(gp.geometry)
src = EPSG(25832)
dst = EPSG(4326)
rpp = ArchGDAL.reproject(px, src, dst)
coords = GeoInterface.coordinates(rpp)
ptc = coords[1]
ptc_tuples = [tuple(pt...) for pt in ptc]
prev = [(Y,X) for (X,Y) in ptc_tuples]
#reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
reversed_polygon = ArchGDAL.createpolygon(prev)

agd = gdf.DataFrame(nm=gp.rcm[1],geometry=reversed_polygon)
plot!(agd.geometry, fillcolor=false)

pol = ArchGDAL.IGeometry[]
for i in 1:size(g,1)
    polygon = g.geometry[i]  # assuming the first geometry is a polygon
    src = EPSG(25832)
    dst = EPSG(4326)
    reprojected_polygon = ArchGDAL.reproject(
        polygon, src, dst) #;order=:compliant
    #z = []
    #for j in reprojected_polygon
        coords = GeoInterface.coordinates(j)
        ptc = coords[1] #muss sein
        ptc_tuples = [tuple(pt...) for pt in ptc]
        prev = [(Y,X) for (X,Y) in ptc_tuples]
        #reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
        reversed_polygon = ArchGDAL.createpolygon(prev)
        push!(z,reversed_polygon)
    end
    push!(pol,z)    
    #push!(pol,reversed_polygon)
end


s = "D:/Wasim/regio/out/rc200/n1/loc2/tst.json"
k = jsread(s)


s="D:/Wasim/regio/out/rc200/n1/loc0/qgesrcm.n1.2015"
cd(raw"D:\Wasim\regio\out\rc200\pest\p12")
s = @gl "qgk"
sa=fread(s)
#sa = @rsubset sa :date >= Date("2015-01-01")
dfp(sa)

s="D:/Wasim/Tanalys/DEM/Input_V2/meteo/pre_ce.txt"
sa=wa.waread(s,station=true)
stpx,stdat=wa.waread(s,station=true,pts_and_data=true,proj=true)
stpx.geometry

@edit waread(s)

"D:/temp/TESTOUT"|>cd
@ncrm

od = selt(stdat,r"kissing"i)|>dropmissing
cdof(s)
pwd()
wawrite(od,"prec_ce_bad_kissingen.wa")
#|>Matrix|>cb

x=range(0,stop=1,length=100)
y=range(1,stop=2,length=100)
[x;y] #<-vertical stacked
[x y] #<-horizontal stacked i.e. matrix
[x y]' #<-transposed
[x;y]'|>reverse


#fliessgeschwindigkeit ws. zu niedrig.
# formel dafür: 
# v = Q/A
# Q = A*v
Q=440 #m3/s
A=1700 #m2
v=Q/A  #unit: m/s 0.26
df = fread("D:/Wasim/regio/out/rc200/pest/p13/qgesrcm.p13.2016")
dfp(selt(df,r"22"))
plot!(df.date,fill(v,length(df.date)),
    label="v=Q/A",xlabel="Date",ylabel="m/s")

vday = v*3600 #m/day

"D:/Wasim/regio/out/rc200/cor1/"|>cd
wa.cmplot()
#"D:/Wasim/Goldbach/out/sglo_v45/post/perc_stack.nc"|>Raster|>plot
rs = Raster("D:/Wasim/Goldbach/out/sglo_v45/post/perc_stack.nc", key=:perc_stack)
rs[Ti=1]|>plot
rs[:,:,3]|>plot
rst.stats(rs)
ad = ncdf(rs)

dx = rs
tmp = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
#m = (X=(lookup(rs,Dim{:easting})),Y=(lookup(rs,Dim{:northing})),Ti=lookup(rs,Ti))

rn = Raster(tmp,
    dims = (X=(lookup(rs,Dim{:easting})),Y=(lookup(rs,Dim{:northing})),Ti=lookup(rs,Ti)))


    using Rasters

    # Extract dimensions
    X_dim = rs.dims[1]
    Y_dim = rs.dims[2]
    Ti_dim = rs.dims[3]
    
    # Create a 3D array with random numbers
    values = rand(X_dim, Y_dim, Ti_dim)
    
    # Create a raster with these values
    raster = Raster(values, dims)

import DimensionalData as ddd
z = ddd.rebuild(rs,values)

#Rasters._maybename(rs)
#setproperty!(rs, :dims, (:X, :Y, :Ti))
rs_renamed = Raster(
    data = copy(rs.data),
    dims = (
        X = rs.dims[1].val,
        Y = rs.dims[2].val,
        Ti = rs.dims[3].val
    )
)

#replace(rs, missing => 0.0)

Rasters.rebuild(rs,refdims=(X=(lookup(rs,Dim{:easting})),Y=(lookup(rs,Dim{:northing})),Ti=lookup(rs,Ti)))
#Rasters.rebuild(rs,dims=(:X,:Y,:Ti))

rcs = rebuild_dimensions(rs)


fn="D:/Wasim/regio/rcm200/pest/v13/ctls/kr1.tsv"
df = CSV.read(fn,DataFrame,delim=':',header=false)
@time_imports setup()
scatter(df.Column4)
sort!(df,4)

x="D:/Wasim/regio/out/pestout/v13/p13/WaSiM_Diag.xml"
@edit wa.extract_duration_from_xml(x)
# run took 2 hrs and 24 min...


using PackageCompiler
using Pkg
pt,nm=(pwd(),"tsjl-win.so")
pts = joinpath(pt,nm)
pkgs = Pkg.installed()
pkgs = Symbol.(keys(pkgs))
create_sysimage(pkgs,sysimage_path=pts)


setup()
fn="D:/remo/cordex/eobs/v28/proj/bak/sfcWind-obs.winlist.bak"
#wa.dubplot(fn)
df = tsread(fn;header=5)
#Date.(reduce(hcat,select(df,1:3)|>eachcol),DateFormat("yyyy mm dd"))
#Date.(map(row -> string(row...), eachrow(select(df, 1:4))), DateFormat("yyyymmddHH"))
df.date = Date.(
    map(row -> string(join(row,"-")), 
    eachrow(select(df, 1:3))),
    DateFormat("yyyy-mm-dd"));
#select!(df, Not(1:4))
#transform!(df, :date => (x -> coalesce.(x, missing)) => :date)
#transform!(df, :proprow => :new)
transform!(df, :Liste => (x -> length.(x)) => :cnt)
selt(df,:cnt)|>wa.dubplot
selt(df,:cnt)|>baryrsum


fn = raw"D:\Wasim\regio\control\new_soil_upper1m"
df=wa.fsoil(fn)
fl = raw"D:\Wasim\regio\control\rcm200_x22-f18-spin-loc.ctl"
fl = raw"D:\Wasim\Goldbach\revision\control\fab90-v2_loc.ctl"
@doc lutab
lu=wa.lutab(controlfile=fl)

laid = wa.luvars(lu)
@df laid boxplot(:class,:LAI,label="LAI",
    fillcolor=:transparent,xrotation=45,bottom_margin = 15mm)

wa.luplot(lu)

findmax(laid.LAI)

"D:/Wasim/Goldbach/revision/fab150/pestpp/v2/output/pestout/qg-files/"|>cd
files = rglob("qgko")

cd("D:/Wasim/Goldbach/out/sglo_v45/")
df = dfr(r"qout")
@doc kge()
kge2(getproperty(df,names(df)[1]),getproperty(df,names(df)[2]))
@pj
pyjl.allkge()


x="D:/Wasim/sinn/control/sinn100._e1_EOBS.ctl"
de = fsoil(x)
thk= extract_layers(de,"thickness")
de.sums|>histogram


fn=raw"D:\Wasim\regio\control\rcm-c9_win_spin5.ctl"
rd=fsoil(fn)
rd.sums|>histogram
rd.sums|>findmin
rd[findmin(rd.sums)[2],[:Name,:key]]

rs = @rsubset rd :sums >= 1


#art1000 -> raw"D:\Wasim\regio\out\rc200\x22\eobs\loc3"
fn = raw"D:\Wasim\regio\control\rcm200_x22-eobs-loc3.ctl"
#rd = wa.read_soildata_4(fn)
v = replace.(grep(r"MultipleHorizons|ksat\s+",readlines(fn)),
    r"\t\t" => " ",r"\t" => " ",r"\s+" => " ")

# Split the rd.Name
split_names = split.(rd.Name)
# Create the DataFrame
df = DataFrame(Dict("n$i" => [split_names[j][i] for j in 1:length(split_names)] for i in 1:3))
ks = hcat(df,extract_layers(rd,"ksat"))
#[i[3] for i in v]
@rsubset ks :n1=="Sl3"



fl=raw"L:\04-Phil\Geo1data\prj-Efre-Daten\Datenmanagement_Modelle\Maps_of_indicators_of_soil_hydraulic_properties_for_Europe\FS_water_UTM_AOI\ths_fao_octop_utm_aoi.tif"
r=Raster(fl)
plot(r)
@cmk
r1 = r[Band=1]
cmk.mkrheat(r1;turn=false)
#Rasters.rplot(r1;padding = (15mm, 15mm, 15mm, 60mm))
Rasters.rplot(r1;title="",Color)

ksr = raw"L:\04-Phil\Geo1data\prj-Efre-Daten\Datenmanagement_Modelle\Maps_of_indicators_of_soil_hydraulic_properties_for_Europe\FS_water_UTM_AOI\ks_fao_utm_aoi.tif"
ksr = Raster(ksr)
#mkecont(r::Raster;wasim=true,missval=0)
#cmk.mkecont(ksr[Band=2];wasim=false,missval=-1.7e308)
cmk.mkrheat(ksr[Band=2];wasim=false,turn=false,missval=-1.7e308)

raw"D:\Wasim\regio\out\rc200\pest\p13"|>cd
dfs = map(dfr,glob(r"qoutjl$")) |>  λ -> map(skipyr, λ) |>  λ ->map(byear, λ)
pxm(dfs)
plot_grouped_metrics(dfs)
dfp(r"gwst")
isroute()

Q=pyjl.allkge()
Q.fn = replace.(Q.fn,r"-qoutjl" => "","_"=>" ")
@vv "tex"

dm = (pwd())|>splitpath|>last
latex_str = sprint(io -> pretty_table(io, Q, backend = Val(:latex),header=names(Q)))
write("KGE_eval-table-$dm.tex", latex_str)
pwc()
dfp("D:/Wasim/regio/out/lowres/c5/loc5/ALL")

df=waread2("D:/Wasim/regio/out/rc200/cl7/tout")
heat(df)
heat(df[!,Cols(r"Bad|Wolf|date")])
se = df[!,Cols(r"Sal|Bad|Wolf|date")]
heat(se)

# case.pdc.csv A summary of prior-data conflict information 
# case.N.pcs.csv A summary of parameter changes
# by group compared to the initial parameter ensemble.
cd("D:/Wasim/regio/out/lowres/pest_864m/main-v2/out")
df = CSV.read("rcm.pdc.csv",DataFrame)
describe(df)
plot(df.distance,df.obs_mean)
se = df[!,Cols(r"dist|obs_mean")]
qplot(se)
qqplot(se[!,1],se[!,2], qqline = :fit)
qqplot(Cauchy,se[!,2])
qqnorm(df[!,2], qqline = :R)

wa.qpl(se)
#title!("Prior-data conflict")
ds = findlog(;lb=.5)
sort!(ds,3,rev=true)
qg = dfr(r"qges")
nr = selt(dfr(ds[1,:nm]),r"C4")
dfp(selt(qg,r"C4"))
wa.dfp!(nr;label="non-routed")
#readobs
obs = dfr("D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu_12")
wlf = selt(obs,r"Wolf")
wa.dfp!(wlf;label="OBS")

xm = mall([selt(qg,r"C4"),nr,wlf])
rename!(xm,1=>:nonrouted,2=>:routed,3=>:observed)
wa.dfp(xm)
heat(xm;c=:transparent,cbar=false)
#qqplot!(qqline = :fit)
#heat(xm;c=:greys)
#heat(xm;c=:white)
wa.heat(xm)
wa.heat(selt(obs,r"mittel|wolf|bad|ober"i))
@edit rst.ezplot()
rst.ezplot(raw"D:/Wasim/regio/rcm200/v15/")
rst.fzplot(raw"D:/Wasim/regio/rcm200/v15/")

cd(raw"D:\Wasim\regio\out\rc200\z0")

v = rglob("qges")
dfs = map(dfr,v)
baryrsum(mall(dfs))

using PyCall
@pyimport datatable as dt
using Conda
Conda.add("xlrd")
@pyimport xlrd
pt=raw"C:\Users\chs72fw\Desktop\OneDrive - Universität Würzburg\divers\annot_63wasim00english.xlsx/edit"
pt=raw"C:\Users\chs72fw\Desktop\OneDrive - Universität Würzburg\divers\annot_63wasim00english.xlsx"
A = dt.fread(pt)
Conda.rm("xlrd")
pwd()

"qge"|>rglob
s = "D:\\Wasim\\regio\\out\\rc200\\z0\\spin2\\qgesrcm.z0.2015"
#awk \'NR==1{print $NF,FS,FILENAME,"\\n","year",FS,"mean",FS,"sum"};$1~/^[0-9]{4}$/&&seen[$1]++{A[$1]++;a[$1]+=$NF}END{for(x in a)print x,a[x]/A[x],a[x]}\'
#dt.fread(;cmd = """awk \'NR==1{print \$NF,FS,FILENAME,"\\n","year",FS,"mean",FS,"sum"};\$1~/^[0-9]{4}\$/&&seen[\$1]++{A[\$1]++;a[\$1]+=\$NF}END{for(x in a)print x,a[x]/A[x],a[x]}\' $s """)
#dt.fread(;cmd = """cat qgesrcm.z0.2015""")

#weil es ja cmd ist. type statt cat.
#fread(cmd=""" type qgesrcm.z0.2015 """,fill=True)
A = dt.fread(cmd=""" type $s """,fill=true).to_pandas()
#or gawk
#fread(cmd=""" awk "/$0~2015/;{print}" qgesrcm.z0.2015 """,fill=True)
#DT=fread(cmd=""" type qgesrcm.z0.2015 | findstr "2015" """)
dt.fread(cmd=""" awk "/\$0~2015/;{print}" $s """,fill=true).to_pandas()
dt.fread(cmd=""" type $s | findstr "2015"  """,fill=true)
pt=raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\wasim\wern\wrn_v3.zip"
A = dt.fread(pt,fill=true).to_pandas()
#A = dt.iread(pt)

#py"A.to_dict('list')"
#D = convert(DataFrame,A.to_pandas())


"D:/Wasim/regio/out/pestout/m_eval/"|>cd
df = dfroute()
kernelplot(df)
klog(df)
dfp("wurz")
dfp("gwst")
dfp("sb05")
facets("sb05")
@cmk
cmk.mkrheat(r"sb0"|>readras)
cmk.tsp(r"sb0"|>dfr)
cmk.tsp(r"qoutjl"|>dfr)
cmk.tsp2(r"Wolf"|>dfr)

fl = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/pre_ce.txt"
ce = waread2(fl)
st,ce = wa.waread(fl;station=true,pts_and_data=true,proj=true)
fn="D:/Wasim/regio/rcm200/v13/cmtv13_4326.json"
import GeoDataFrames as gdf
p = gdf.read(fn)
import GeoInterface as gi
#gi.extent(p.geometry)
et = reduce(gdf.union,p.geometry)
et .+ 0.1
?buffer(et,0.1)
import ArchGDAL as ag
ptsub = ag.buffer(et,0.1) #.5
ptsub = wa.reverse_coords(ptsub)
@vv "within"
#point in polygon
v = [ gdf.within(i, ptsub) for i in st.geometry ] 
stsub = st[[ gdf.within(i, ptsub) for i in st.geometry ],:]
#subset df colnames 
#v = join([string(stsub.name),:date],"|")
v = vcat(stsub.name,:date)
v = string.(v)
v = replace.(v,r"-" => "_","__"=> "_")
#first.(split.(v,"-"))
#get only first 4 characters
v = [i[1:4] for i in v]
v = join(v,"|")
v = Regex(v)
cesub = select(ce,Cols(v))
ofile = replace(fl,".txt" => "_sub.txt")
#wawrite(cesub,ofile)
cesub = @rsubset cesub year.(:date) >= 1980
ofile = replace(fl,".txt" => "_smallsub.txt")
wawrite(cesub,ofile)
@edit baryrsum(cesub)
#dropmissing(cesub)
wa.baryrsum(cesub;leg=false,grid=false,ticks=false)
#tline(cesub)
wa.baryrsum(selt(cesub,r"Wuer"i);leg=false)

cmk.tsbar(selt(cesub,r"Wuer"i))
dfm(selt(cesub,r"Wuer"i);fun=yrsum,ann=false,mode=:bar)
rmm(selt(cesub,r"Wuerzbu"i);)|>dfp

#get a mean of cesub
#transform!(cesub, rowmeans -> :submean)
#numeric_cols = names(cesub, Real)
#cesub.submean .= mean.(eachrow(cesub[:, Not(:date) .∈ numeric_cols]))
#To replace missing values in a DataFrame with 0 in Julia, you can use the `coalesce.` function. Here's how you can do it:
#- `coalesce.(cesub, 0)` replaces all missing values in `cesub` with 0. The `.` after `coalesce` is used to broadcast the function over the entire DataFrame.
cesub .= coalesce.(cesub, 0) 
# dat = tryparse.(Float64,string.(cesub[:, Not(:date)]))
# dat.mean = mean.(eachrow(dat))
cesub.submean .= mean.(eachrow(cesub[:, Not(:date)]))
println(names(cesub))
cesub = @rsubset cesub year.(:date) >= 1970
baryrsum(selt(cesub,r"submean"))
cmk.tsbar(selt(cesub,r"submean"))

selt(cesub,r"submean")|>tline
df2 = ce
df2 .= coalesce.(df2, 0)
df2.pmean .= mean.(eachrow(df2[:, Not(:date)]))
selt(df2,r"pmean")|>tline
#stimmt nicht, weil 0 mit in den trend berechnet wird.
#besser jede station einzeln.
st,ce = wa.waread(fl;station=true,pts_and_data=true,proj=true)
ptsub = et
ptsub = wa.reverse_coords(ptsub)
@vv "within"
#point in polygon
stsub = st[[ gdf.within(i, ptsub) for i in st.geometry ],:]
v = vcat(stsub.name,"date")
v = replace.(v,"." => "_","-" => "_")
v = replace.(v,"__"=>"_")
nsu = select(ce,Cols(v))
nsu = @rsubset nsu year.(:date)>=1980
#p1 = Plots.plot();
p1 = tline(select(nsu,Cols([1,ncol(nsu)]))|>dropmissing);
for i in 2:size(nsu)[2]-1
    tmp = select(nsu,Cols([i,ncol(nsu)]))|>dropmissing
    tline!(tmp)
end
Plots.title!("Niederschlag Trends")

wa.baryrsum(nsu;leg=false)
ofile = replace(fl,".txt" => "_smallsub.txt")
wawrite(nsu,ofile)
##better reverse:
npt = ag.createpoint.(
ag.gety.(stsub.geometry,0),
ag.getx.(stsub.geometry,0))
stsub.geometry .= npt
#write header
plot(stsub.geometry)
plot!(et,fillcolor=:transparent)
sufile = replace(fl, ".txt" => "_hdrs")
writedf(sufile,stsub)
#wa.revcrds(stsub)
hdf = DataFrame(
    YY=fill("YY",4),
    MM=fill("MM",4),
    DD=fill("DD",4),
    HH=fill("HH",4)
)
#dat = stsub.name
#reduce(vcat,vcat)
#höhe X,Y,name
dat = DataFrame([1:nrow(stsub),stsub.xc,stsub.yc,stsub.name],:auto)
dat = permutedims(dat)
rename!(dat,stsub.name)
outdf = hcat(hdf,dat)
sufile = replace(fl, ".txt" => "crds_hdrs")
writedf(outdf,sufile)

#now merge via npp
#ofile = replace(fl,".txt" => "_smallsub.txt")
@edit stread(ofile)
k=stp(ofile;proj=true)
plot(k.geometry)

fn=ofile
fl = CSV.read(fn,DataFrame;limit=4)
xc = fl[2,5:end]|>collect
yc = fl[3,5:end]|>collect

####
@vv "zrxp"
vgr("zrxp")
fn="D:/Wasim/Pegeldaten/2018_08_08_Bodenwasserhausaltsmodellierung_Franken/Q_15Min_Teil3/24409003__Wolfsmünster__Q__Prod__2018-08-07__mit_Freigabe.zrxp"
df = CSV.read(fn,DataFrame;skipto=9,header=1)
nh=read_delim(file = lnk,delim = " ", 
skip = 8, 
col_types = 
cols( X1 = col_date(format = "%Y%m%d %H%M %S"), 
X2 = col_double(), X3 = col_skip() ), 
col_names = F, na = "-777", trim_ws = T, 
progress = T, escape_backslash = F) 

using CSV
#dateformat"yyyymmdd HHMM SS"
nh = CSV.read(fn, DataFrame,
               delim = " ", 
               skipto = 9,
               limit = 100,
               header = false, 
               #types = [Date, Float64, String],
               #types = [Date, Float64],
               types = Dict(1=>DateTime,2=>Float64),
               dateformat = dateformat"yyyymmddHHMMSS",
               select = [1,2],
               missingstring = "-777",
               silencewarnings = true)
#[Date(d, dateformat"u d, y") for d in df.date]
# hdr = CSV.read(fn, DataFrame,delim="|",limit=7
#     ,select = [1,2,3])
#headers:
replace.(split.(readlines(fn)[2],"|")[[3,5]],
   "\xfc"=>"ue","\xe4"=>"ae","\xdf"=>"ss")

replace.(readlines(fn)[1:8],
"\xfc"=>"ue","\xe4"=>"ae","\xdf"=>"ss",
"\xf6"=>"oe","\xc4"=>"Ae","\xd6"=>"Oe",
"\xb3"=>"³")

"""
reads .zrxp Pegeldaten
"""
function qread(fn::String;drop=false)
    # hdrs=replace.(readlines(fn)[1:8],
    #     "\xfc"=>"ue","\xe4"=>"ae","\xdf"=>"ss",
    #     "\xf6"=>"oe","\xc4"=>"Ae","\xd6"=>"Oe",
    #     "\xb3"=>"³")
    # println(hdrs)
    df = CSV.read(fn, DataFrame,
               delim = " ", 
               skipto = 9,
               limit = 2_000_000,
               header = false, 
               types = Dict(1=>DateTime,2=>Float64,3=>String),
               dateformat = dateformat"yyyymmddHHMMSS",
               select = [1,2,3],
               missingstring = "-777",
               silencewarnings = true)
    hd=replace.(split.(readlines(fn)[2],"|")[[3,5]],
               "\xfc"=>"ue","\xe4"=>"ae","\xdf"=>"ss")
    rename!(df,1=>"date",2=>hd[1],3=>hd[2])
    if drop
        select!(df, Not(3))
    end
    return df
end

df = qread(fn;drop=true)
@time setup()
dfp(selt(df,r"Wolf"))
df=selt(df,r"Wolf")
cdof(fn)
baryrsum(df)
# s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
# @df df Plots.plot(:date,cols(s))
# plot(df.date,df.SNAMEWolfsmuenster)

"D:/Wasim/Tanalys/DEM/Input_V2/meteo/"|>cd
function tswrite(df::DataFrame,file::AbstractString;HH::Int64=0)
    dout = copy(df)
    if in("year",names(dout))
        @warn "yearcol found!"
        CSV.write(file, dout, 
        transform = (col, val) -> something(val, missing), delim="\t")  
        return
    end
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout.HH = map(x ->hour(x),dout.date)
    dout.mm = map(x ->Minute(x).value,dout.date)
    
    dout = dout[!,Cols([:YY,:MM,:DD,:HH,:mm],Not(Cols(r"date")))]
    CSV.write(file, dout, 
    transform = (col, val) -> something(val, missing), delim="\t")  
    nothing
end
tswrite(df,"wlf_15min_raw.wa")
npplat()
qts="c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe"
inf = "wlf_15min_raw.wa"
onf = "wlf_dly_spec.wa"
run(`$qts`)
run(`$qts $inf $onf 5 all 24 1`)
npplat() #nope, wsl cut
pwc()
a=wlf_15min_raw.wa
cut -f1-4,6- $a > wlf_15
npplat()

inf = "wlf_15"
onf = "wlf_dly_spec.wa"
run(`$qts $inf $onf 5 all 24 1`)
npplat()

first(df)
A=2125.700
q=11.6
q/A
q*60*60*24/A/1000


onf2="wlf_m3.wa"
run(`$qts $inf $onf2 5 all 1 1`)
npplat()

da = nread(onf2)
max_value, max_index = findmax(da.Wolfsmuenster)
@view da[max_index, :]
#	440 m3/s

A=2125.700
Q=0.745215
Q*A
60*60*24*Q/(A*1000)
@view df[findmax(df.SNAMEWolfsmuenster)[2], :]

Q=440/A

sp = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu_12"
obs = waread2(sp)
obs = (selt(obs,r"Wolf"))
@view obs[findmax(obs.Wolfsmuenster)[2], :]

s = raw"D:\Wasim\Tanalys\DEM\Input_V2\meteo\wlf_dly_spec.wa"
df = wa.waread(s)
@view df[findmax(df.Wolfsmuenster)[2], :]


"""
Calculates the runoff rate (MQ) given the mean discharge and catchment area.

Parameters:
mean_discharge (Float64): The mean discharge in m³/s.
catchment_area (Float64): The catchment area in km².

Returns:
Float64: The runoff rate (MQ) in l/s/km².
"""
function calculate_runoff_rate(mean_discharge, catchment_area)
    # Convert catchment area from km² to m²
    catchment_area_m2 = catchment_area * 1_000_000
    # Calculate the runoff rate (MQ) in m³/s/m²
    runoff_rate_m3sm2 = mean_discharge / catchment_area_m2
    # Convert the runoff rate to l/s/km²
    runoff_rate_lskm2 = runoff_rate_m3sm2 * 1_000_000_000
    return runoff_rate_lskm2
end
calculate_runoff_rate(440, 2125.7)

mean_discharge = 440  # m³/s
mean_discharge = 16.3  # m³/s
catchment_area = 2125.7  # km²
runoff_rate = calculate_runoff_rate(mean_discharge, catchment_area)
println("Runoff rate (MQ): $runoff_rate l/s/km²")


"""
Calculates the mean discharge given the runoff rate and catchment area.

Parameters:
runoff_rate (Float64): The runoff rate (MQ) in l/s/km².
catchment_area (Float64): The catchment area in km².

Returns:
Float64: The mean discharge in m³/s.
"""
function calculate_mean_discharge(runoff_rate, catchment_area)
    # Convert runoff rate from l/s/km² to m³/s/m²
    runoff_rate_m3sm2 = runoff_rate / 1_000_000_000
    # Convert catchment area from km² to m²
    catchment_area_m2 = catchment_area * 1_000_000
    # Calculate the mean discharge in m³/s
    mean_discharge = runoff_rate_m3sm2 * catchment_area_m2
    return mean_discharge
end

#runoff_rate = 7.6  # l/s/km²
runoff_rate =  17.8852 # l/s/km²
catchment_area = 2125.7  # km²
mean_discharge = calculate_mean_discharge(runoff_rate, catchment_area)
println("Mean discharge: $mean_discharge m³/s")


##aggregate to daily TS; fun = mean
using Dates

df_daily = copy(df)
# Convert the date column to Date type (removing the time part)
df_daily.date = Date.(df_daily.date)
#df_daily = DataFrames.combine(groupby(df, :date), :2 => mean)
@view df_daily[findmax(df_daily.SNAMEWolfsmuenster_mean)[2], :]
@rsubset df_daily :date >= DateTime(2003, 1, 3) && :date <= Date(2003, 1, 4)
su = @rsubset df_daily :date >= DateTime(2003, 1, 3) && :date <= Date(2003, 1, 6)
su|>dfp

Dates.Day.(df.date)[6].value

df.days = Dates.Day.(df.date)
df.mon = Dates.Month.(df.date)
df.yr = Dates.Year.(df.date)
#df.dt = Dates.Date.([df.yr,df.mon,df.days])
#map(x->value(x),[df.yr,df.mon,df.days])
dty = [i.value for i in df.yr]
dtm = [i.value for i in df.mon]
dtd = [i.value for i in df.days]
df.dt = Dates.Date.(dty,dtm,dtd)

df_daily = DataFrames.combine(groupby(df, :dt), :2 => mean)
@view df_daily[findmax(df_daily.SNAMEWolfsmuenster_mean)[2], :]
su = @rsubset df :date >= DateTime(2003, 1, 3) && :date <= Date(2003, 1, 4)
mean(tovec(su,2))

rename!(df_daily, :SNAMEWolfsmuenster_mean => :Wolfsmuenster,:dt=>:date)
pwd()
wawrite(df_daily,"wlf_daily.wa")

#transform!(df_daily, :Wolfsmuenster => (x -> calculate_runoff_rate(x,2125.700)) => :MQ)

"""
Calculates the runoff rate in mm/day given a list of discharge values and the catchment area.

Parameters:
discharge_values (Vector{Float64}): A list of discharge values in m³/s.
catchment_area (Float64): The catchment area in km².

Returns:
Float64: The runoff rate in mm/day.
"""
function calculate_runoff_rate_mmday(discharge_values, catchment_area)

    # Calculate the mean discharge
    mean_discharge = sum(discharge_values .* 1000) / length(discharge_values)

    # Convert catchment area from km² to m²
    catchment_area_m2 = catchment_area * 1_000_000

    # Calculate the runoff rate in mm/day
    runoff_rate_mmday = (mean_discharge * 86400) / catchment_area_m2
    #runoff_rate_mmday = (mean_discharge * 5760) / catchment_area_m2
    #to mm/day
    runoff_rate_mmday 

    return runoff_rate_mmday
end

calculate_runoff_rate_mmday(dfx[1].wlf,2125.700)
#dfx = DataFrames.combine(groupby(df, :dt), :2 => (x -> calculate_runoff_rate_mmday(x,2125.700)))
#@view dfx[findmax(dfx.SNAMEWolfsmuenster_function)[2], :]

rename!(df, :SNAMEWolfsmuenster => :wlf)
dfx = groupby(df, :dt)
dfx[1]
#map(x->findmax(x.wlf),dfx)
q=[]
for i in dfx
    tmp=DataFrame(wlf=calculate_runoff_rate_mmday(i.wlf,2125.700),date=i.dt[1])
    push!(q,tmp)
end
q = vcat(q...)
@view q[findmax(q.wlf)[2], :]
ofl="wlf_jlspec.wa"
wawrite(q,ofl)

#transform!(df_daily, :Wolfsmuenster => (x -> calculate_runoff_rate_mmday(x,2125.700)) => :TS)
#@view df_daily[findmax(df_daily.MQ)[2], :]
#plot(df_daily.date,df_daily.TS)
#@view df_daily[findmax(df_daily.MQ)[2], :]
#wawrite(select(df_daily,[1,3]),"wlf_spec.wa") #s.o.
#rm("wlf_spec.wa") 

#jetzt das Selbe mit QtoSpend
qts="c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasim_helptools_time/qtospend.exe"
inf = "wlf_daily.wa"
onf = "wlf_dly_wasim.spec"
run(`$qts`)
#usage: qtospend <discharge_m3/s> <spec_disch._mm/dt> 
#[<n1 headrows{def=5}>] [<n2 datacolumns{def=all}>] 
#[<timestep in h{def=1h}>] [<method 1(to mm/step) or 2(to m3/s) {def=1}>
run(`$qts $inf $onf 5 all 24 1`)
npplat()

waspec = waread("wlf_dly_wasim.spec")
jlspec = waread("wlf_jlspec.wa")
@view jlspec[findmax(tovec(jlspec,1))[2], :]
@view waspec[findmax(waspec.Wolfsmuenster)[2], :]

fu = mall([jlspec,waspec])
fu3 = @rsubset fu year.(:date) == 2003
dfp(fu3;yaxis=:log)

#das wäre falsch
#transform!(fu3 , :Wolfsmuenster => (x -> calculate_runoff_rate(x,2125.700)) => :MQ)
#@view fu3[findmax(fu3.MQ)[2], :]

#mm/day to m3/s
x = sum(fu3.Wolfsmuenster .* 1000) / length(fu3.Wolfsmuenster)
calculate_runoff_rate(x,2125.700)

sp = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu_12"
obs = waread2(sp)
obs = (selt(obs,r"Wolf"))
@view obs[findmax(obs.Wolfsmuenster)[2], :]

obs2 = waread("specdis_kmu.txt")
obs2 = (selt(obs2,r"Wolf"))
@view obs2[findmax(obs2.Wolfsmuenster)[2], :]
13.25/(60*60*24)*2125.7*1000

##!! zurückrechnen auf m3/s!!!!
transform!(obs2, :Wolfsmuenster => 
    (x -> (x/(60*60*24)*2125.7*1000)) => :q)
dfp(obs2;yaxis=:log)

dfm(selt(df_daily,"MQ"))
dfm(df_daily)

xm = mall(df_daily,obs2)
qplot(select(xm,Cols(r"Wolfsmuenster|q")))

select(xm,Cols(r"Wolfsmuenster$|q|date"))|>dfm

#obs2===obs
#obs2.Wolfsmuenster .- obs.Wolfsmuenster|>sum

obs = waread2(sp)
sb = selt(obs,"Steinbach")
pr = selt(obs,r"parte"i)
@view sb[findmax(sb.Steinbach)[2], :]

dfm(obs)
npp(sp)

@doc wa.waread
pts = wa.waread(sp;station=true)
sort!(pts,:ez)
#
A = 21519
A = @rsubset pts :name=="Steinbach"
A = only(A.ez)
transform!(sb, :Steinbach => 
    (x -> (x/(60*60*24)*A*1000)) => :q)

@view sb[findmax(sb.q)[2], :]
@view sb[findmin(sb.q)[2], :]

idx = argmax(sb.q)
@view sb[idx, :]

@view sort(sb, :q)[end, :]

"D:/Wasim/Tanalys/DEM/brend_fab/out/w3/penman/re/"|>cd
ctlf = ctl2()
ofl="route.txt"
routeg(ctlf,ofl)
readlines(ofl)
dfroute(;ofl)


###sinn
sp = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu_12"
obs = waread2(sp)
si = selt(obs,r"Gemuen")
@view si[findmax(tovec(si,1))[2], :]

stn = wa.waread(sp;station=true,proj=true,pts_and_data=false,src=EPSG(25832),dst=EPSG(4326))
#names(si)[1]
@rsubset stn occursin(r"Gem",:name)

A = 620.3
##!! zurückrechnen auf m3/s!!!!
transform!(si, names(si)[1] => 
    (x -> (x/(60*60*24)*A*1000)) => :q)

@view si[findmax(tovec(si,1))[2], :]

s=raw"D:\Wasim\regio\control\abflusstafel.txt"
replace.(grep(r"^M|Sinn - Haupt",readlines(s)),r"\s+" => " ")
#x = CSV.read(s,DataFrame,header=false,transpose=true)
@doc wa.readbetween
#x = wa.readbetween(s,r"^[Sinn - H]",r"HQ")
# x = wa.readbetween(s,r"kst",r"Main$")
# x = replace.(x,r"\s+" => " ")
# x = split.(x)
replace.(readlines(s)[17:23],r"\s+" => " ")

@doc readdlm(s)
CSV.read(s,DataFrame,header=41,skipto=42)

raw"D:\Wasim\sinn\out\y1"|>cd
@pyjl
ds = pyjl.allkge()
dfs = map(dfr,glob(r"qoutjl$")) |>  λ -> map(skipyr, λ) |>  λ ->map(byear, λ)
pxm(dfs)
plot_grouped_metrics(dfs)
plot_grouped_metrics(dfs,col=:kge)
plot_grouped_metrics(dfs,col=:nse)
pxm(dfs,threshold=0.10)

#vcat(map(x->x[findmax(x.kge)[2], :],dfs))
reduce(vcat, map(x -> DataFrame(x[findmax(x.kge)[2], :]), dfs))

using RCall
@rimport terra as tr



x=raw"C:\Users\chs72fw\Desktop\OneDrive - Universität Würzburg\divers\Param_Magdali.txt"
#Evaluation Statistic  
df = CSV.read(x,DataFrame;skipto=271,footerskip=165,delim=' ',header=270,stripwhitespace=true,ignoreemptyrows=true)
bar(df.PBIAS)
df.logNSE|>plot
df.NSE|>plot!

fn=raw"C:\Users\chs72fw\Desktop\OneDrive - Universität Würzburg\divers\Lfu_ezgs.txt"
#D = DataFrame(grep("änkische Saale",readlines(fn)))

lines = readlines(fn)
filtered_lines = filter(x -> occursin("änkische Saale", x), lines)
M = split.(filtered_lines, " ", limit=3)
D = DataFrame(Line = filtered_lines)

fn=raw"C:\Users\chs72fw\Desktop\OneDrive - Universität Würzburg\divers\Bodendaten_Proben.csv"
#da = gdf.read(fn)
da = CSV.read(fn, DataFrame)|>dropmissing
#using GeoDataFrames, GeoInterface
# Convert RECHTSWERT and HOCHWERT to Float64
#da.RECHTSWERT = parse.(Float64, da.RECHTSWERT)
#da.HOCHWERT = parse.(Float64, da.HOCHWERT)
# Create a new column that contains Point objects
#da[!,:geometry] = [ArchGDAL.createpoint(da.RECHTSWERT[i], da.HOCHWERT[i]) for i in 1:nrow(da)]

da[!,:geometry] = [ArchGDAL.createpoint(
    da.HOCHWERT[i],
    da.RECHTSWERT[i] 
    ) for i in 1:nrow(da)]
#gk4 <-CRS("+init=epsg:31464") #-->GK4, (DWD)
#inplace
[ArchGDAL.reproject(pt,EPSG(31464),EPSG(4326)) for pt in da.geometry]
#[ArchGDAL.reproject(pt,EPSG(31468),EPSG(4326)) for pt in da.geometry]
#[ArchGDAL.reproject(pt,EPSG(32632),EPSG(4326)) for pt in da.geometry]
#@edit rst.project(fn)epsg:31468

da.lat .= [ArchGDAL.getx(pt,0) for pt in da.geometry]
da.lon .= [ArchGDAL.gety(pt,0) for pt in da.geometry]

dx = select(da,Not(:geometry))
dn = dirname(fn)
writedf(dx,joinpath(dn,"BodenProben.tsv"))
db = select(da,Not([1,2]))
db = select(db, Not(:geometry), :)
#gdf.write(joinpath(dn,"BodenProben.geojson"),db; driver="GeoJSON")
#gdf.write(joinpath(dn,"BodenProben.shp"),db)
#s.C:/Users/Public/Documents/Python_Scripts/rfile/ezgchecks.R 
# z = read.table(z,header=TRUE,sep="\t")
# p = vect(z,geom=c("lon", "lat"),crs="+proj=longlat +datum=WGS84")
# pd = crop(p,el) %>% terra::as.data.frame(.,geom="XY")
fn = "C:/Users/chs72fw/Desktop/OneDrive - Universität Würzburg/divers/BodenProben_ez.csv"
#df = CSV.read(fn,DataFrame)
@pyimport pandas as pd
df = pd.read_csv(fn,skip_blank_lines=true, encoding="cp1252")
#make uniqe by lat,lon
df = df.drop_duplicates(subset=["x","y"]) #240
using StringEncodings
df = open(fn, enc"CP1252") do io
    CSV.File(io;decimal=',') |> DataFrame
end

using DataFramesMeta
#filter(x -> occursin(r"Schluff"i, x), df.BODENFORM)
#subset(df, :BODENFORM, :y => ByRow(!))
@subset(df, occursin.(r"Schluff", :BODENFORM))

DataFrames.subset(df,:BODENFORM => ByRow(==("Sch")))
DataFrames.subset(df, :BODENFORM => ByRow(x -> occursin(r"Sch", x)))
DataFrames.subset(df, 
    3 => ByRow(x -> occursin(r"Gem", x)) ,
    5 => ByRow(z -> !occursin(r"gley"i, z)))

subset(df, 
    3 => ByRow(x -> occursin(r"Grab", x)) ,
    5 => ByRow(z -> !occursin(r"gley"i, z)))

@edit findindf(df,"acker")
#StatsBase.countmap(df.G_HOEHE, df.HOCHWERT)
#StatsBase.countmap(df.HOCHWERT,df.G_HOEHE)
pwd()

using RCall
libstr="C:/Users/chs72fw/AppData/Local/R/win-library/4.2"
@rput libstr
#RCall.@R_str(".libPaths($libstr)")
R""".libPaths(new=libstr)"""
@rimport terra as tr
##2015 lucal landcover#########
fn="D:/Relief/landcover/Land_Cover_DE_2015.tif"
r =Raster("D:/Relief/landcover/lc_ufra.tif")
crs(r)
fl="LULC_DE_2014_nbg_200m.tif"
fl="DFD_LULC_DE_2014_v1_clip_AOI_Franken_Domain.tif"
fl="LULC_DE_2014_lanczos_200m.tif"
fldr="C:/Users/chs72fw/Documents/EFRE_GIS/Landcover/20181016_Deutschland_LC_clip_for_Ullmann"
r2=Raster(joinpath(fldr,fl))
cd(fldr)
sf("nbg")
r1=tr.rast(fn)
#tr.setMinMax(r1);
r2=tr.rast(first(sf("nbg")))
r1p=tr.project(r1,r2)
r1=tr.crop(r1,r1p)
#r3=tr.rast(tr.ext(r2),res=200,crs=tr.crs(r2))
#r3=tr.project(r3,tr.crs("epsg:25832"))
r3=tr.rast(tr.ext(r3),res=200,crs=tr.crs(r3))
r3 = tr.project(r1,r3,method="near")
ls()
#rm("D:/Relief/landcover/Land_Cover_2015_200m_25832.tif")
tr.writeRaster(r3,"D:/Relief/landcover/Land_Cover_2015_200m_25832.tif")
#rout=tr.project(r3,tr.crs("epsg:25832"),method="near")
#r3=tr.setMinMax(r3)

df=fread("D:/Wasim/regio/out/rc200/y7/l1/tout")
@edit skipyr(df)

ftp(df)
cdof(df)
mkdir(raw"D:\Wasim\regio\out\rc200\allqb")
cd(raw"D:\Wasim\regio\out\rc200\allqb")
@doc rglob()
@edit rglob("2x")
v = wa.rglob("qbas","D:/Wasim/regio/out/rc200/x22")
v = filter(x->!occursin(r"wkly",x),v)
#l = map(wcl,v)
dfs = []
for i in v    
    try 
        push!(dfs,fread(i))
    catch
        println("$i errord")
        continue
    end
end
dfs

xm=map(wa.yrsum,dfs)
xm
wa.filterplot("12",xm)

df = dfr("D:/Wasim/regio/out/rc200/qbfile/allqb")
df = skipyr(df)
baryrsum(df;leg=false)

a = Dates.Minute(8) * 200
#convert to hours
#a = div(a.value, 60) # integer division
#convert to hours
a = round(a.value / 60) # rounding to nearest integer
a = Dates.Hour(a)

a = Dates.Second(60) * 1600
#convert to hours
a =  Dates.Hour(round(a.value / 60 / 60)) # rounding to nearest integer


"D:/Wasim/regio/out/rc200/y7/l3"|>cd
dfs = map(dfr,glob(r"qoutjl$")) |>  λ -> map(skipyr, λ) |>  λ ->map(byear, λ)
pxm(dfs)
plot_grouped_metrics(dfs)
dfp(r"gwst")
isroute()

@pyjl
Q=pyjl.allkge()
Q.fn = replace.(Q.fn,r"-qoutjl" => "","_"=>" ")
dm = (pwd())|>splitpath|>last
latex_str = sprint(io -> pretty_table(io, Q, backend = Val(:latex),header=names(Q)))
write("KGE_eval-table-$dm.tex", latex_str)

fzplot("D:/Wasim/regio/rcm200/y7/x0")
savefig("D:/Wasim/regio/rcm200/y7/x0/fzplot.png")
oplat()

pyjl.pyplot_df(r"gwst")

#10272.8 minutes in hours
a = Dates.Minute(10273)
a = Dates.Hour(round(a.value / 60)) # rounding to nearest integer
#a = Dates.Hour(div(a.value, 60))   # integer division
Dates.Day(div(a.value, 24))         # integer division

div(1220.6,60)
round(1220.6 / 60 / 24;digits=2)
round(10273 / 60 / 24;digits=2)

#sf rcm.rec|xargs -I% grep -iwH 'took'  %
@cmk
using .cmk
fn="D:/Wasim/regio/out/rc200/y7/g3/uprsrcm_1200.sum.nc"
@doc mkrheat(fn;ubound=0.0)
mkrheat(fn;umask=0.1,mskval=-10.001,layer=1)
mkrheat(fn;umask=-4.6145383e-14,mskval=-100,layer=1)
mkrheat(fn;umask=-1,mskval=-100,layer=1,missval=0)

tsp(df)
df2 = dropmissing(df,2)
tsbar(selt(df,2))
tsbar(selt(df,r"schlim"i))
cloudplot(df)

"D:/Wasim/sinn/pestout/p1/"|>cd
dfs = map(dfr,glob(r"qoutjl$")) |>  λ -> map(skipyr, λ) |>  λ ->map(byear, λ)
pxm(dfs)
plot_grouped_metrics(dfs;col=:kge)
#@edit mkrheat(r"wind";mskval=.001,layer=1,umask=9.1,missval=0)
Main.cmk.mkrheat(r"wind";mskval=.001,layer=1,umask=9.1,missval=0)
mkrheat(r"vap";)
mkrheat(r"tsoi"i;layer=2,umask=20.1,mskval=-1000.001)
fn="D:/Wasim/regio/out/rc200/y22/rf-f16/precrcm_1300.sum.nc"
cmk.mkrheat(fn;layer=1,mskval=20.001)
cdof(fn)
cmk.mkrheat(r"perc";)
sf(r"perc";)
@doc regand
cmk.mkrheat(r"perc+.*14")
cmk.mkrheat(r"perc+.*16")
cmk.mkrheat(r"perc+.*17+.*sum")
cmk.mkrheat(regand("perc","17"))

"D:/Wasim/sinn/pestout/p1/"|>cd
cmk.mkrheat(r"perc+.*17+.*sum")
cmk.mkrheat(r"perc+.*16+.*sum")
op()
cd("D:/Wasim/regio/out/rc200/y12/f4/")
cmk.mkrheat(r"perc")
#cmk.tsplot((r"perc")|>dfr)
dfp("perc")
cd("D:/Wasim/regio/out/rc200/y12/f6")
cmk.mkrheat(r"perc")
cmk.mkrheat(r"pre+.*14")
cmk.mkrheat(r"gw")

#R = Raster("D:/Wasim/regio/out/rc200/y12/f6/sbst.nc")
x="D:/Wasim/regio/out/rc200/y7/l4/precrcm.2013.nc"
cdof(x)
@cmk
using .cmk
mkrheat(x)
v = sf("perc")
mkrheat(v[3])
@edit mkrheat(r"gwst";maskval=-100.001,layer=1)
mkrheat(r"gwst";mskval=-100.001,layer=1)
mkrheat(r"gwstrcm_";mskval=-10,umask=-1,hide=true)

x="D:/Wasim/regio/out/rc200/y12/f6/precrcm_1300.sum.nc"
mkrheat(x)

cd(("E:/qq/rf2/"))
glob("pre")
#r = Raster("pre_80_ens_0002.nc",key=:pre)
#ERROR: NetCDF error: Variable 'time_bnds' not found in file pre_80_ens_0002.nc 
#ncap2 -O -s 'defdim("bnds",2);time_bnds=make_bounds(time,$bnds,"time_bnds");' $a p_04.nc
#xrt p_04.nc pre

x="E:/qq/rf2/p_04.nc"
r = Raster(x,key=:pre)
mkrheat(r,layer=1345)
df = rst.ncdf(r)
using DataFramesMeta
#dfut = @rsubset df :date >= "2050-01-03T12:00:00"
dfut = @rsubset df year.(:date) >= 2050
cmk.tsbar(dfut)
baryrsum(dfut)



@pyjl
xr = pyimport("xarray")
ds = xr.open_dataset("pre_80_ens_0002.nc")
ds.keys()
#add time_bnds
time_bnds = ds["time"].copy()
ds["time_bnds"] = xr.DataArray(time_bnds, 
    dims=["bnds"], 
    coords=Dict("bnds"=>[0]))   
#ERROR: KeyError: key "time_bnds" not found
#ds.expand_dims(time_bnds.values)
ds = ds.close()


####works too
using NCDatasets
# Open your NetCDF file
ds = Dataset("pre_80_ens_0006.nc", "a")
# Check if the TIME variable exists
if haskey(ds, "time")
    # Get the time variable
    time = ds["time"]
    # Create a new dimension called 'bnds'
    defDim(ds, "bnds", 2)

    # Create a new variable called 'time_bnds'
    time_bnds = defVar(ds, "time_bnds", time, ("time", "x", "y"))

    # Calculate the bounds for each time value
    # # This is just an example, you might need to adjust this for your specific case
    # for (i, t) in enumerate(time)
    #     time_bnds[i, 1] = t - 0.5
    #     time_bnds[i, 2] = t + 0.5
    # end
else
    println("The time variable does not exist in the NetCDF file.")
end
# Close the NetCDF file
close(ds)

r=Raster("E:/qq/rf2/pre_80_ens_0005.nc",:pre)
#r=Dataset("E:/qq/rf2/pre_80_ens_0005.nc")
close(r)


k="pre_80_ens_0007.nc" 
#rst.add_time_bounds(k)
ds = Dataset(k, "a")
#renameVar(ds, "bnds", "time_bnds")
renameDim(ds, "time_bnds","bnds")
close(ds)

Main.rst.add_time_bounds(k)
#r=Raster("E:/qq/rf2/pre_80_ens_0007.nc",:pre)

k="pre_80_ens_0002.nc" 
ds = Dataset(k, "r")
ds.variables
#renameVar(ds, "bnds", "time_bnds")
close(ds)
cp(k,"tmp.nc",force=true)
ds = Dataset("tmp.nc","a")
time=ds["time"]
#NCDatasets.defVar(ds, "time_bnds", time, ("time"))
#NCDatasets.defVar(ds, "time_bnds", time, ("lon", "lat", "time"))
#NCDatasets.defVar(ds, "time_bnds", time, 1)
NCDatasets.defVar(ds, "time_bnds", time, ("time", "x", "y"))
close(ds)


using NCDatasets

# Open the old dataset
ds_old = NCDataset("tmp.nc", "r")
# Create a new dataset
ds_new = NCDataset("new_file.nc", "c")

# Define the new dimensions in the desired order
defDim(ds_new, "time", length(ds_old["time"]))
defDim(ds_new, "lon", length(ds_old["lon"]))
defDim(ds_new, "lat", length(ds_old["lat"]))

pre=ds_old["pre"]

# Define a new variable with the reordered dimensions
#var_new = defVar(ds_new, "pre", pre, ("lon", "lat", "time")) #errors
# Copy the data from the old variable to the new one
#var_new[:,:,:] = permutedims(ds_old["pre"][:,:,:,], [3, 1, 2])

# Close the datasets
close(ds_old)
close(ds_new)

"D:/Wasim/regio/rcm200/v4/"|>cd
@time_imports setup()
fzplot()
Main.rst.ezplot(;cbar=false)
savefig("D:/Wasim/regio/rcm200/v4/ezplot-jl.png")
oplat()

cd("E:/qq/")
sf("tas-qq2")

r=Raster("E:\\qq\\jlcor\\tas-qq2.nc",key=:tas)
plot(r[Ti=6])
##exkurs Rcall json###############
@rimport terra as tr
ve = tr.vect("D:/Wasim/regio/rcm200/v4/ezg.shp")
ve = tr.project(ve, "+proj=longlat +datum=WGS84 +no_defs")
tr.writeVector(ve,"D:/Wasim/regio/rcm200/v4/ezg_4326.json")
############
import GeoDataFrames as gdf
ge = gdf.read("D:/Wasim/regio/rcm200/v4/ezg_4326.json")
ge = gdf.read("D:/Wasim/regio/rcm200/v4/ezg.shp")
rcmoutline = reduce(gdf.union,ge.geometry)
plot!(ge.geometry,color=:transparent)

using NCDatasets
fn=raw"D:\remo\proj_main\25832_daymean_rh_2000.nc"
ds = Dataset(fn, "r")
close(ds)

fn="D:/remo/cordex/eobs/proj/hu_ens_mean_utm.nc"

rh=Raster(fn,key=:hu)
plot(rh[Ti=666])
plot!(ge.geometry,color=:transparent)
mrh = mean(rh,dims=Ti)
plot(mrh[Ti=1])
plot!(ge.geometry,color=:transparent)
#plot(rmh[X(Near(9.9)),Y(Near(50.0))],margins=5mm)
xn = mean(extrema(lookup(rh,:X)))
yn = mean(extrema(lookup(rh,:Y)))
plot(rh[X(Near(xn)),Y(Near(yn))],margins=5mm)


#subset by time while reading...
rh=Raster(fn,key=:hu;mappedcrs=EPSG(25832),crs=EPSG(25832))[Ti(Where(x->year.(x)>=2000))]
xn = mean(extrema(lookup(rh,:X)))
yn = mean(extrema(lookup(rh,:Y)))
plot(rh[X(Near(xn)),Y(Near(yn))],margins=5mm)


hu=23
√(100-hu)
#lt eobs site
daily averaged relative humidity HU
The dataset for relative humidity is based on the in-situ data holdings of the ECA&D dataset similar as the other E-OBS datasets. The gridding method used for this element is the same as the one used for temperature, precipitation and sea level pressure. To remove some of the skewness in the data, the relative humidity values (and the background field used in the gridding method) were transformed by 
prior to fitting. This also ensures that all interpolated values, when converted back to the unit of %, are equal or smaller than 100%
mean(√, [1, 2, 3])
mean(√, [9,9,9])
#!!! note If itr contains NaN or missing values, the result is also NaN or missing
mrh = mean(rh,dims=Ti) #so this is not what we want
rh = mask(rh,with=ge.geometry)
plot(rh[Ti=100])
#mrh = mean(rh,dims=Ti)
mrh = mean(rh,dims=(X,Y))
#plot(mrh[X=1,Y=1])
#ncdf(mrh)
#lookup(mrh,:X)
@cmk
using .cmk
mkrheat(mrh)

mrh = project(mrh;dst="EPSG:4326")
plot(mrh[Ti=1])

#Highly recommended, lakemodel groundwater model

(1,       2,       59 ) .* (0.1,    0.1,     0.3333)


########das ist mist.
fn="D:/remo/cordex/eobs/v28/rcm_eobs/rh_wa.nc"
rh=Raster(fn,key=:rh)
dims(rh)[1]
lookup(rh,rh[1,:,:].refdims)
rh
contourf(rh[t=3])
#;mappedcrs=EPSG(25832),crs=EPSG(25832)
rne = project(rh[t=3];dst="EPSG:25832")


fn="D:/remo/cordex/eobs/v28/proj/rh-obs.nc"
#subset by time while reading...
rh=Raster(fn,key=:rh;mappedcrs=EPSG(25832),crs=EPSG(25832))[Ti(Where(x->year.(x)>=2000))]
xn = mean(extrema(lookup(rh,:X)))
yn = mean(extrema(lookup(rh,:Y)))
plot(rh[X(Near(xn)),Y(Near(yn))],margins=5mm)

plot(rh[X(Near(xn)),Y(Near(yn)),
    Ti(Where(x->year.(x)==2000))],
#    yaxis=:log,
    margins=5mm,title="EOBS 28e")

plot(rh[Ti(Near(DateTime( today() )))],margins=5mm)
plot!(ge.geometry,color=:transparent)

rhm = mean(rh[Ti(Where(x->year.(x)==2022))],dims=(Ti))
contourf(rhm,margins=5mm)
#plot(rhm,margins=5mm)
plot!(ge.geometry,color=:transparent)

##mean over all dimensions
nme = mean(skipmissing(reshape(rh[Ti(Where(x->year.(x)==2022))], :)))
x = [1, 2, missing, 4, 5]
mean(skipmissing(x))
mean(x)

kras = rh[Ti(Where(x->year.(x)==2017))]
kras = replace_missing(kras,0)
kras = mean(kras,dims=(Ti))
contourf(kras,margins=5mm)
plot(kras,margins=5mm)
plot!(ge.geometry,
    color=:transparent,
    linecolor=:lightgrey,linesize=1.9,linestyle=:dashdot)

    #moved to rst.rmeanplot
function meanplot(ras,yr=2000;contour=true,kwargs...)
    kras = ras[Ti(Where(x->year.(x)==yr))]
    kras = replace_missing(kras,0)
    kras = mean(kras,dims=(Ti))
    if contour
        contourf(kras,
            margins=5mm;
            colorbar_title = "mean values of the year $yr",
            #colorbar_title_padding = 5,
            kwargs...)
    else
        plot(kras,margins=5mm;
        colorbar_title = "mean values of the year $yr",
        kwargs...)
    end
    #plot(kras,margins=5mm)
    plot!(ge.geometry,
        color=:transparent,
        linecolor=:lightgrey,linesize=1.9,linestyle=:dashdot)
end

meanplot(rh,2003;title="EOBS 28e",contour=false)
meanplot(rh,2003;title="EOBS 28e",contour=true)
rst.rmeanplot(rh,2022,ge;contour=false)
rst.rmeanplot(rh,2022;contour=true)


# yrs = rh.dims|>last|>size|>first 
# yrs = yrs ÷ 365
# cy = sum(Rasters.trim(Rasters.crop(rh;to=rcmoutline, 
#     touches=true);pad=2),dims=Ti) ./ yrs
# #preobs = sum(obsras,dims=tdim) ./ yrs
# plot(cy,margins=5mm)





@pyjl
readdir("d:/remo/qm/rh/")
##humidity "doyplot"
xr  = pyimport("xarray")
sh="d:/remo/qm/rh/simh.nc"
simh=xr.open_dataset(sh)
#sh="d:/remo/qm/rsds/simh.nc"
sh="d:/remo/qm/rh/simp.nc"
simp=xr.open_dataset(sh)
readdir("D:/remo/cordex/eobs/v28/rh/")
tl="D:/remo/cordex/eobs/v28/rh/rh-obs.nc"
obsh=xr.open_dataset(tl)
@edit 
pyjl.doyplot(simh,simp,obsh;mytitle="EOBS")

#D:/Wasim/regio/out/rc200/x4/loc3/tsoilrcm__stack.2015.nc
"D:/remo/cordex/eobs/v28/proj/rh-obs.nc"



fn="D:/Wasim/regio/out/rc200/x4/loc4/tsoilrcm__stack.2015.nc"
r = Raster(fn)
plot(r[t=2])
dims(r)
plot(r[t=27])
cdof(fn)
sb = readras(r"sb05")
plot(sb[t=1])
qgk()
nsegrep()
#=
Multi line 
comments 
=#

#warscher
#Water_balance_estimation_in_high_Alpine_terrai.pdf
#Comparison of the water storage (Smod) derived from results of hydrological model runoff before ANN correction. Monthly sums
#June–October 2002–2011. 
