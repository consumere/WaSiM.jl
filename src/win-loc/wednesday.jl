
function monsum(x::String)
    df = readdf(x)
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    df[!, :month] = month.(df[!,:date]);
    df_monthsum = combine(groupby(df, :month), y .=> sum .=> y);
    return(df_monthsum)
end

function monsum(x::DataFrame)
    df = x
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    df[!, :month] = month.(df[!,:date]);
    df_monthsum = combine(groupby(df, :month), y .=> sum .=> y);
    return(df_monthsum)
end

function monmean(x::String)
    df = readdf(x)
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    df[!, :month] = month.(df[!,:date]);
    df_monthsum = combine(groupby(df, :month), y .=> mean .=> y);
    return(df_monthsum)
end

function monmean(x::DataFrame)
    df = x
    y = filter(x->!occursin("date",x),names(df))
    s = map(y -> Symbol(y),y)
    df[!, :month] = month.(df[!,:date]);
    df_monthsum = combine(groupby(df, :month), y .=> mean .=> y);
    return(df_monthsum)
end

function barp(x::DataFrame)
    "with DataFrame input"
        df = x
        ti=DataFrames.metadata(df)|>only|>last|>basename
        if any(x->occursin("year",x),names(df))
            ln = Symbol.(filter(x->!occursin("year",x),names(df)))
            @df df Plots.plot(:year,
                cols(ln),
                legend = :topright, 
                title=ti,
                seriestype=:bar) #color=:lightrainbow
        elseif any(x->occursin("month",x),names(df))
            ln = Symbol.(filter(x->!occursin("month",x),names(df)))
            @df df Plots.plot(:month,
                cols(ln),
                legend = :topright, 
                title=ti,
                seriestype=:bar)
        elseif (
            any(x->occursin("month",x),names(df)) & 
            any(x->occursin("year",x),names(df))            
            )
            ln = (filter(x->!occursin("month",x),names(df)))
            ln = Symbol.(filter(x->!occursin("year",x),ln))
            @df df Plots.plot(:month,
                cols(ln),
                legend = :topright, 
                title=ti,
                seriestype=:bar)
        else
            dfp(df)        
        end
end

            
            #s = "apple banana"
            # occursin("apple", s) & occursin("banana", s)


"/mnt/d/Wasim/regio/out/lowres/c6/c5-stats"|>cd
df=readdf(r"clo")
df = yrsum(df)



regand(names(df),"month","year")



sfvar sb*nc
julia --startup-file=no --color=no --threads auto -q --optimize=0 -e "
using NCDatasets;
nc=NCDataset(ARGS[1]);
println(nc)
" $ff

#printstyled(nc,color=:blue)

julia --startup-file=no --color=yes --threads auto -q --optimize=0 -e "
using NCDatasets;
nc=NCDataset(ARGS[1]);
println(nc)" $ff


jlo --startup-file=no --color=yes --threads auto -q --optimize=0 -e "
using NCDatasets;
nc=NCDataset(ARGS[1]);
println(nc)" $ff




nc=NCDataset(ARGS[1])
printstyled(nc,color=:blue)

# sfvar sb*nc
# julia  --startup-file=no --color=yes --threads auto -q --optimize=0 -e "



pt="/mnt/c/Users/chs72fw/Documents/EFRE_GIS/Landcover/20181016_Deutschland_LC_clip_for_Ullmann/landuse_midpoints.geojson"
using GeoJSON
jsonbytes = read(pt);
fc = GeoJSON.read(jsonbytes)
df = DataFrame(fc)
describe(df)
first(fc)

fn="/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/w1/w1_nclogjl.txt"
df = CSV.read(fn,DataFrame,delim="\t")

#Plots.plot(df[:,end-1],df[:,end])
#x,y=df[:,1],df[:,2]

x,y=df[:,end-1],df[:,end]
s = propertynames(df)[end-1:end]
StatsPlots.plot( 
qqplot(x,y, qqline = :fit), 
qqplot(Cauchy,y), 
qqnorm(y, qqline = :R),title=String.(s)) 



#
fd()
cd("/mnt/d/Wasim/regio/out/lowres/c9")
map(parent, dims(z, (X, Y)))

z=readallras(r"gwst")
zs=stack(z)
zr = Raster(zs;dims=(X, Y, Z, Ti),missingval=0)
typeof(zr)
typeof(zs)

gr()
Plots.plot(zr[Ti=(2)])
zr[Ti=(1)]|>describe
zr[Ti=(9)]|>Plots.plot


###########how to stack
#zs=readallras(r"gwn")|>stack
zs=readallras(r"uprs")|>stack
zw = dropdims(zs; dims=3)

zstack = Raster(zw;dims=(Y, X, Ti),missingval=0)
#zstack = Raster(zw;dims=(X,Y,Ti),missingval=0)

zstack[Ti=(9)]|>Plots.plot
zstack[Ti=(2)]|>Plots.plot
zstack[Ti=(3)]|>Plots.contourf
# zstack[Ti=(3)]|>transpose|>Plots.plot
# A=transpose(zstack)
###########how to shp
using Shapefile
shp="/mnt/d/Relief/DLM250/dlm250_ufrasubset/DLM250_Fliessgew_Ufra.shp"
riv = Shapefile.Handle(shp)
riv.shapes[1:end]|>Plots.plot
riv.crs
tr = readras(r"albe")
cut = mask(tr; with=riv.shapes[1:end])
Plots.plot(cut)

rivsub=mask(riv; with=tr)
m=cut.dims|>DataFrame #|>Matrix
#m[1:2,:]|>Matrix
xm = (cut.dims[1],cut.dims[2],cut.dims[3])
typeof(xm)
map(parent,xm)
#tempr = Raster(xm;dims=(X,Y,Ti))
#tempr = Raster(map(parent,xm))
#Plots.plot(tempr)
Plots.plot(xm)

########four panels plots######
# lk https://github.com/dankelley/R2julia/blob/main/R2j.Rmd
using Plots
n = 10
x = 1:n
l = @layout[a b; c d]
a = plot(x, rand(n))
b = plot(x, 2*x)
c = plot(x, x.*(x.-n))
d = plot(x, exp.(-x))
plot(a, b, c, d, layout=l)



cd("/mnt/d/Wasim/regio/out/rc200/v0")
using Plots
bw=readdf("sobw.v0")
dfyrs(bw)
names(bw)

d=readdf("sobw.v0")
names(d)

dyr = bw

println("calculating yearly water balance of WaSiM special output data..\n
bwvr = rain + snow + uprs - perc - qb - qd - qi - etr_ - ei_ - etrs_\n")

#posbal = findfirst(r"^(Perc)", names(dyr)) - 1
posbal = findfirst(x -> occursin(r"^(Perc)", x), names(dyr))
propertynames(dyr)
dyr[!,Cols(r"^(prec)|^(Snow)|^(Capi)|year")]

pos = select(dyr,names(dyr)[posbal+1:end])
neg = select(dyr,Cols(names(dyr)[2:posbal],"year"))


#transform!(pos, [:year], :sum .=> ByRow(cumsum) => :cum_sum)
#transform!(neg, [:year], :sum .=> ByRow(cumsum) => :cum_sum)

df = DataFrame(A = [1, 2, 3], B = [4, 5, 6], C = [7, 8, 9])

df = pos
names(df)|>print
y = filter(x->!occursin(r"date|year",x),names(df))
df_csum = combine(groupby(df, :year), y .=> sum .=> y);
# Compute cumulative sum for each row of columns A and B
df = df_csum[!,Not(:year)]
# calculate the sum of each row
df_sum = DataFrame(sum_each_row = [sum(eachrow(df)[i]) for i in 1:size(df, 1)])
# add the new column to the original dataframe
hcat(df, df_sum)

# #cumsum(df)
# x=[]
# for row in eachrow(df)
#     x = cumsum(row,Cols(Not(:year)))
# end
# x
df = pos
y = filter(x->!occursin(r"date|year",x),names(df))
df_csum = combine(groupby(df, :year), y .=> sum .=> y)
# Compute cumulative sum for each row of columns A and B
df = df_csum[!,Not(:year)]
# calculate the sum of each row
psum = DataFrame(
    possums = [sum(eachrow(df)[i]) for i in 1:size(df, 1)],
    year=df_csum[!,:year]
)

df = neg
y = filter(x->!occursin(r"date|year",x),names(df))
df_csum = combine(groupby(df, :year), y .=> sum .=> y);
# Compute cumulative sum for each row of columns A and B
df = df_csum[!,Not(:year)]
# calculate the sum of each row
nsum = DataFrame(
    negsums = [sum(eachrow(df)[i]) for i in 1:size(df, 1)],
    year=df_csum[!,:year]
)


bw = innerjoin(psum, nsum, on=:year)
bw.bw = bw[!,:possums] .- bw[!,:negsums]
select!(bw, [:year, :bw])
ti=DataFrames.metadata(df)|>only|>last|>basename
@df bw Plots.plot(:year,:bw,
legend = :topright, 
seriestype=:bar,
xticks = bw.year,
xlabel = "Year",
ylabel = "bwvr[mm]",
title = ti,
fillcolor = ifelse.(bw.bw .> 0, "cornflowerblue", "coral2"))


###how to subset best...
x="/mnt/d/Wasim/regio/out/rc200/v0/sobw.yearv0"
df = CSV.read(x,DataFrame,delim="\t",comment="-")
df[!,filter(x->occursin(Regex("Sum","i"),x),names(df))]
df[!,Cols("SUM.x")]
df[!,filter(x->startswith(x,"S"),names(df))]
df[!,("SUM.x")]
pos
df[!,Cols(r"S")]        ##!!!
df[!,Not(Cols(r"S"))]   #negative

"v12"|>cd
dout = readdf(r"q")
dout.YY = map(x ->year(x),dout.date)
dout.MM = map(x ->month(x),dout.date)
dout.DD = map(x ->day(x),dout.date)
dout[!, "HH"] .= 24
#df = select!(df,Symbol.(filter(x->!occursin("date",x), names(df))))
#dout = select(df, Not(:date))
dout = dout[!,Cols([:YY,:MM,:HH,:DD],Not(Cols(r"date")))]
#cls = propertynames(df)|>sort|>reverse
#df = df[!,cls[2:end]] 
CSV.write(file, dout, transform = (col, val) -> something(val, missing),delim="\t")  
# df = readdf(r"ev")
# writewa("tst",df)

#dout[!,Not(Cols(r"date"))]


lk="/mnt/d/Wasim/Goldbach/stateini/evaporation_modis_500/2022/ctlcum"
nd=readdf(lk)

data = nd[!,Not(:date)]
rename!(data,1=>"x",2=>"y")
using LinearFitXYerrors, DataFrames
linearfitxy(data.x, data.y, isplot=false, ratio=:auto)


setup()
s="/mnt/d/Wasim/regio/rcm200/fine_soil_table.tsv"
df = CSV.read(s,DataFrame,
    types = [Int,Int,String],
    delim="\t")

#s = Symbol.(filter(x->!occursin("date",x),names(df)))
#o = DataFrames.metadata(df)|>collect
ti = "AndrewsPlot of "*basename(s)
cls = df[!,Not(Cols(r"bf"))]|>propertynames
@df df andrewsplot(:bf_id, cols(cls), legend = :topleft,title=ti)

bl = copy(df)

findall(bl.horizons[:] .>= 3)
findall(bl.horizons .> 3)
bl[bl.horizons .> 3,:]
#DataFrames.subset(df, :horizons => ByRow(isodd) .& :horizons => ByRow(>(3)))
#DataFrames.subset(df, :horizons => ByRow(isodd) .& :horizons => ByRow(==(3)))

filter(:horizons => x -> x % 2 == 0, df) #even
# filter(:horizons => x -> x % 2 == 0  :bf_id => e -> e .> 1000, df) #error
# filter((:horizons => x -> x % 2 == 0) .& (:bf_id => e -> e .> 1000), df) #error

filter(:bf_id => x -> x .> 1000, df)
df[(df.horizons .> 3) .& (300 .< df.bf_id .< 999), :]

using Shapefile
fn = "/mnt/d/Bodendaten/buek200_2020/bk200_sfmerge.shp" 
shp = Shapefile.Handle(fn)

shp.shapes[1:3]|>Plots.plot
shp.crs
tr = readras(r"albe")
cut = mask(tr; with=shp.shapes[1:end])
#Plots.plot(cut)

Pkg.add(["Shapefile", "Clipper"])


crop_poly = map(Matrix, dims(tr, (X, Y)))

#ext = extend(tr)
xmin = dims(tr, (X, Y))[1]|>first
xmax = dims(tr, (X, Y))[1]|>last
ymin = dims(tr, (X, Y))[2]|>first
ymax = dims(tr, (X, Y))[2]|>last


min(dims(tr, (X)))
tr.data[:,:,1][1]

findall(tr.data[:,:,1] .>= 1)
tr.data[:,:,1]|>Plots.plot


#trp = points(A::AbstractRaster; dims=(YDim, XDim), ignore_missing) 
trp = points(tr; ignore_missing=true)|>first 
points(tr; ignore_missing=true)|>last

# Roughly cut out New Zealand from the evenness raster
crop_poly = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]]
rr = tr[t=1]
#nz_bounds = X(165 .. 180), Y(-50 .. -32)
#nz_evenness = evenness[nz_bounds...]
bnds = X(xmin .. floor(xmax)), Y(ymin .. floor(ymax))
nz = rr[bnds...]
#Plots.plot(nz)

sh = shp.shapes[bnds...]

ceil(xmin)

# Crop range to match evenness
nz_range = crop(rnge; to=nz_evenness)
plot(nz_range)



typeof(trp)
#crop_poly = Shapefile.Polygon()
#Shapefile.Polygon(crop_poly) wrong.

#fpolygon = Point.(tofloat.(crop_poly, 3,3))
clipped = Clipper.Clip(shp.shapes[1], crop_poly)
cropped_shapes = []
using Clipper
for shape in Shapefile.shapes(shp)
    clipped = Clipper.Clip(shape.coordinates[1], crop_poly)
    if !isempty(clipped)
        push!(cropped_shapes, Shapefile.Polygon(clipped))
    end
end



sub=mask(shp; with=tr)
m=cut.dims|>DataFrame #|>Matrix
#m[1:2,:]|>Matrix
xm = (cut.dims[1],cut.dims[2],cut.dims[3])
typeof(xm)
map(parent,xm)


"/mnt/d/Wasim/regio/rcm200/rcm.art_bfid"|>rplot
"/mnt/d/Wasim/regio/rcm200/rcm.art_mglo"|>rplot

fs = df[!,3]
Plots.histogram(fs)

s="/mnt/d/Wasim/regio/rcm200/fine_soil_table.tsv"
fs = CSV.read(s,DataFrame,
    #types = [String],
    normalizenames=true,drop=(i, nm) -> i == 3,
    delim="|")

"/mnt/d/Wasim/main/tmp/sb05rcm_1000.mit.nc"|>rplot


using Shapefile
#fn = "/mnt/d/Bodendaten/buek200_2020/bk200_sfmerge.shp" 
fn = "/mnt/d/Wasim/regio/rcm200/ezg.shp"
shp = Shapefile.Handle(fn)
Plots.plot(shp)
rr = readras("/mnt/d/Wasim/main/tmp/sb05rcm_1000.mit.nc")
cut = trim(mask(rr; with=shp.shapes[1:11]))
Plots.plot(cut)

#shp|>Plots.plot

pt="/mnt/c/Users/chs72fw/Downloads/iknmi14_tas_Amon_ECEARTH23_rcp85_6-12E_45-55N_n_++_a.txt"
ms = "#"
df = CSV.read(pt,DataFrame,missingstring=ms,
    ntasks=4, limit=typemax(Int), ignorerepeated=true,
    skipto = 39,
    delim=" ", #    comment="#",
    types = Float64, #    types = [Float64,Float64],
    silencewarnings=true,normalizenames=true)
    select!(df,1:2)
    rename!(df,1=>"date",2=>"t")


#,drop=(i, nm) -> i == 4) |> dropmissing
@df df Plots.plot(:date,:t)

data = dropmissing(df[!,1:2])

d1 = readdf()

st = linearfitxy(data.t, isplot=true, ratio=:auto)


using Plots

bardf("/mnt/d/Wasim/regio/out/rc200/r1/so_precipitation.r1.2017")

pr = readdf("/mnt/d/Wasim/regio/out/rc200/r1/so_precipitation.r1.2017")
pr = yrsum(pr)
gr()
Plots.bar(  pr[:,2],
             title = "precipitation",
             color = ifelse.(pr[:, 2] .> 0, "cornflowerblue", "coral2"),
             xlabel = "time",
             legend = false, 
             ylabel="Precipitation [mm]")
xticks!(pr[:, :year])

annotate!(pr[:,2] ./ 1.5, round.(pr[:, 2], digits = 2), textsize = 0.66)

pr[:,2] == pr[!,2]

hline!([0])

Plots.plot(1:10)
Plots.bar(pr[:,2])
annotate!([2,5], [6,3], ["text at (2,6)", "text at (5,3)"])


using Gadfly
using Cairo, Fontconfig
px = Gadfly.plot(pr, x = :year, y = Symbol.(names(pr)[2]), Geom.bar);
draw(PNG("plot.png", 6inch, 4inch), px)
#Plots.savefig(px,"plot.png")

pt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/sun_rel.txt"
df = readdf(pt)
#df = map(x->ismissing(x),df)
eltype(df.date) 
typeof(df.date) 
filter(df,-9999.0)

nd = findall(df[!,1] .>= 0)
findall(df[!,Not(end)] != -9999.0)

#filter(x->occursin(-9999.0,x),df)
vgjl("findall")

pt="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/sonnenstunden_1970.txt"
df = readdf(pt)

fn="/mnt/d/remo/qm/tas_delta.nc"
#["time", "longitude", "latitude", "tas"]
tt = Tuple(("time", "longitude", "latitude", "tas"))
rr = Raster(fn,dims=tt)
#ncpdq --reorder 'time','longitude','latitude','tas' $a tas-reorder.nc
ncks --trd $a tas-reorder.nc
fn="/mnt/d/remo/qm/tas-reorder.nc"
rr = Raster(fn,dims=(:lat,:lon,:time))
rr
using NCDatasets
nc = NCDatasets.Dataset(fn)
println(nc) 
# Read dimensions
lon = nc["longitude"][:]
lat = nc["latitude"][:]
time = nc["time"][:]
# Load all data of a variable
tas = nc["tas"][:]
Dict(nc.attrib)
# Load the attributes
xx = nc["tas"]
ad = Dict(xx.attrib)
# Define day and depth level of interest
date = DateTime(2000, 12, 15,12)
# Extract the indices 
time_indices = findall(date .== time)[1]   
# Generate the plot 
Plots.heatmap(lon, lat, tas[:,:,time_indices]', 
        xlabel = "Longitude", ylabel = "Latitude", 
        plot_title = string("tas on ", date))

#tas[:,:,time_indices]|>Plots.plot
tas[:,:,time_indices]'|>Plots.plot

#rr = Raster(tas[lon,lat,time_indices])

"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/c6"|>cd
dfp(r"gwst")
dfp(r"sb05")
dfp(r"ts_av")
# results soil temperature thaw depth or active layer thickness as average value for subbasins
bardf(r"ts_loc")
bardfm(r"win")
xd = readdf(r"win")
# rxp qd__fab_1100.sum.nc 0
# bmerge 
# bmerge 8
# hydlat
# bmerge
# bmerge 7
# bmerge 6
# bmerge 
# bmerge 5
kge_read("qout")
qgg()
rr=Raster("/mnt/d/Wasim/regio/rcm/rcm.use2";missingval=-9999)
plotlyjs()
Plots.plot(rr)

"/mnt/d/Wasim/regio/out/rc200/r2"|>cd

#readdf("SCNWIND_rcm.r2.2017")
dras = readallras(".")
filterplot()

"/mnt/d/Wasim/regio/out/rc200/r2"
find . -iname "*$1*" -type f -not -regex '.*.\(png\|svg\|txt\|html\|ftz\|ftz_0\|log\|list\|nc\|xml\|list\|sh\|grd\|yrly\)$' | while read f; do
        perl -pi -e 's/\x00//g' $f;done;

find . -maxdepth 1 -type f -not -name "so*" -not -regex '.*.\(html\|xlm\|ftz_0\|nc\|txt\)$' -print |\
    xargs -I{} perl -pi -e 'print "YY\tMM\tDD\tHH\t_1_\t_2_\t_3_\t_4_\t_5_\t_6_\t_7_\t_8_\t_9_\t_10_\t_11_\t_12_\t_13_\t_14_\t_15_\t_16_\t_17_\t_18_\t_19_\t_20_\t_21_\t_22_\t_23_\t_24_\t_25_\t_26_\t_27_\t_28_\t_29_\t_30_\t_31_\t_32_\t_33_\t_34_\t_35_\t_36_\t_37_\t_38_\t_39_\t_40_\t_41_\t_42_\t_43_\t_44_\t_45_\t_46_\t_47_\t_48_\t_49_\t_50_\t_51_\t_52_\t_53_\t_54_\t_55_\t_56_\t_57_\t_58_\t_59_\t_60_\t_61_\t_62_\t_63_\t_64_\t_65_\t_66_\t_67_\t_68_\t_69_\t_70_\t_71_\t_72_\t_73_\t_74_\t_75_\t_76_\t_77_\t_78_\t_79_\t_80_\t_81_\t_82_\t_83_\t_84_\t_85_\t_86_\t_87_\t_88_\t_89_\t_90_\t_91_\t_92_\t_93_\t_94_\t_95_\t_96_\t_97_\t_98_\t_99_\t_100_\t_101_\t_102_\t_103_\t_104_\t_105_\t_106_\t_107_\t_108_\t_109_\t_110_\t_111_\t_112_\t_113_\t_114_\t_115_\t_116_\t_117_\t_118_\t_119_\t_120_\t_121_\t_122_\t_123_\t_124_\t_125_\t_126_\t_127_\t_128_\t_129_\t_130_\t_131_\t_132_\t_133_\t_134_\t_135_\t_136_\t_137_\t_138_\t_139_\t_140_\t_141_\t_142_\t_143_\t_144_\t_145_\t_146_\t_147_\t_148_\t_149_\ttot_average\n" if $. == 1' {} && \   
find . -maxdepth 1 -type f -name "so_*" -not -regex '.*.\(html\|xlm\|ftz_0\|nc\|txt\)$' -print | xargs -I{} perl -pi -e 'print "YY\tMM\tDD\tHH\tlyr1\tlyr2\tlyr3\t\n" if $. == 1' {} && echo "headers done...try h so*"

#jetzt
dfs = readall(".")
filterplot("wind",dfs)
ct("qgk")
df=readdf(r"qgk")
dfp(r"qges")

gr()
x = randn(1000) # generate some random data
StatsPlots.plot(x, seriestype = :density) # plot a histogram
StatsPlots.qqnorm!(x)
StatsPlots.plot!(x, seriestype = :qq)

vgjl("Gadfly")
vgr("explorer")

using Gadfly
#Gadfly.Theme()
myt=Theme(panel_fill=colorant"white", panel_stroke=colorant"white")
Gadfly.push_theme(myt)
#Gadfly.get_theme(myt)
Gadfly.current_theme()
p1 = Gadfly.plot(x = x, Geom.histogram) # plot a histogram
p2 = Gadfly.plot(x = x, Geom.density); # plot a density plot
p3 = Gadfly.plot(x = x, Geom.bar); # plot a quantile-quantile plot
p = [p1,p2,p3]
#gridstack(reshape([[render(p[i]) for i in 1:3], canvas()],2,2))
gridstack([p1 p2; p3 p1])

Gadfly.push_theme(myt)
p1 = Gadfly.plot([sin,cos], 0, 2pi)
p2 = Gadfly.plot([sin,cos], 0, 2pi, style(line_width=2mm, line_style=[:dash]))
fig = hstack(p1,p2)
Gadfly.pop_theme()

"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/qbas"|>cd
dfs=readall(".")
filterplot("q",dfs)
mm = mall(dfs)

using RDatasets
Gadfly.with_theme(:dark) do
    Gadfly.plot(RDatasets.dataset("datasets", "iris"), x=:SepalLength, y=:SepalWidth, color=:Species)
end
iris=RDatasets.dataset("datasets", "iris")
Gadfly.plot(iris, x=:SepalLength, y=:SepalWidth, Geom.point, Geom.line)

mammals = RDatasets.dataset("MASS", "mammals")
p1 = Gadfly.plot(mammals, x=:Body, y=:Brain, label=:Mammal, Geom.point, Geom.label);
p2 = Gadfly.plot(mammals, x=:Body, y=:Brain, label=:Mammal, Geom.point, Geom.label,
     Scale.x_log10, Scale.y_log10);
hstack(p1, p2)
vstack(p1, p2)

vgr("baseflow")


Gadfly.plot(iris, x=:SepalLength, y=:SepalWidth,Stat.binmean(n=10))

df = dfs[end]
Gadfly.plot(df, x=:tot_average, y=:_2,Stat.binmean(n=50))
D=df
p1 = Gadfly.plot(D,x=:tot_average, y=:_1, Geom.point, intercept=[0], slope=[1], Geom.abline(color="red", style=:dash));
abline = Geom.abline(color="red", style=:dash)
p2 = Gadfly.plot(D, 
    x=:tot_average, 
    y=:_1,  
    Geom.point,  
    Scale.x_asinh, Scale.y_log,
     intercept=[148], slope=[-0.5], abline);
hstack(p1, p2)

homg
fd
mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;

cd("/mnt/d/Wasim/Goldbach/revision/qbfiles")
#dfs = readall(".")
dfs = readall(r"fab150")

#mm=mall(dfs)
#gridstack(reshape([[render(p[i]) for i in 1:3], canvas()],2,2))
function dfsplog(dfs::Vector{DataFrame})
    "plots and adds"
    df = dfs[1]
    s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    p = @df df Plots.plot(:date,
            cols(s),
            yaxis = :log,
            legend = false)
            #legend = :bottom)
    for i in 2:length(dfs)
        nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
        println("adding $nm")
        s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
        @df dfs[i] Plots.plot!(:date,cols(s),
        label="$nm") #geht, wenn oben legend true ist.
        # label="$nm",
        # legend = false)
        # Plots.annotate!(0.5, 0.5, text(nm, 14))
    end
    display(p)
end

function dfsp(dfs::Vector{DataFrame})
    "plots and adds"
    df = dfs[1]
    s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    p = @df df Plots.plot(:date,
            cols(s),
            legend = false)    #legend = :bottom)
    for i in 2:length(dfs)
        nm=DataFrames.metadata(dfs[i])|>only|>last|>basename
        println("adding $nm")
        s = Symbol.(filter(x->!occursin(r"date|year",x),names(dfs[i])))
        @df dfs[i] Plots.plot!(:date,cols(s),
        label="$nm") #geht, wenn oben legend true ist.
    end
    display(p)
end


dfsp(dfs)
mx=readall(r"fab30")
dfsp(mx)

#sfo 45|xargs head -1
mx=readall(r"v45")
dfsp(mx)
map(describe,mx)
mx=readall(r"fab23")
dfsp(mx)
mx=readall(r"fab150")
dfsp(mx)
dm = map(describe,mx)
#xm = mall(mx)

mx=readall(pwd())
dfsp(mx)




p = Plots.plot(rand(10), title="Multiple plots")
for i in 2:3 Plots.plot!(rand(10), label="Plot $i") end
display(p)

df=mx[1]
s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
@df df Plots.plot!(:date,cols(s),label="lsl")

M = Matrix(df[!,Not(:date)])
corrplot(M)
@df df corrplot(cols(s), grid = false)


homes
fd
mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
pwc

cd("/mnt/d/temp/saale/output/qbfiles")
dfs = readall(".")
#dfs = readall(r"fab150")
#dfsplog(dfs)
dfsp(dfs;save="allbas")

homreg
cd out
mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
pwc
cd("/mnt/d/Wasim/regio/out/qbfiles")
dfs = readall(".")
#dfs = readall(r"fab150")
#dfsplog(dfs)
dfsp(dfs;save="allbas")
ll()
dfs = readall(r"b26")
dfsp(dfs;save="allbas")

#Der Hauptunterschied zwischen den beiden Backends ist, dass Plotly() eine Internetverbindung und einen Plotly-Account erfordert, um die Plots zu erstellen und zu hosten, während PlotlyJS() dies nicht tut. 
plotly()
dfs = readall(r"c8")
dfsp(dfs)

mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
pwc
cd("/mnt/d/Wasim/Testdata/qbfiles")
dfs = readall(".")
dfsp(dfs;save="allbas")

mx=mall(dfs)
names(mx)
vio(mx)

hombr
mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
pwc
cd("/mnt/d/Wasim/Tanalys/DEM/brend_fab/qbfiles")
dfs = readall(".")
dfsp(dfs)
dfsp(dfs;save="allbas")

mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
pwc
cd("/mnt/d/Wasim/Tanalys/qbfiles")
dfs = readall(".")
dfsp(dfs;save="allbas")

dfs = readall(r"100")
dfsp(dfs)

dfs = readall(r"2000")
dfsp(dfs;save="bas2000")

readall(r"1991")|>dfsp
readall(r"br50")|>dfsp
readall(r"qbasbr50.m1997")|>dfsp
dfp("qbasbr50.m1997_DEM-HAMONTEST.wa")
dfp("qbasbr50.m1997_DEM-PENMANTEST-penman_w_routing.wa")
"qbasbr50.m1997_DEM-out_v3_vali-v4.wa"|>dfp
"qbasbr50.m1997_DEM-HAMONTEST-hamon_param_GW.wa_qbfiles.wa"|>dfp

"/mnt/d/Wasim/Tanalys/DEM/HAMONTEST"|>cd
dfs=readall(pwd())
dfsplog(dfs)
filterplot("ra",dfs)
ddx = getdf("ra",dfs)|>yrsum
getdf("ra",dfs)|>dfp
bardf(r"rain")

ve = ct("qbas")

md = readall(ve)
dfsp(md)
rplot(r"qd")
rplot(r"wind")
rplot(r"tsoil")
ct("q")

rplot(r"qi")



mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
pwc
cd("/mnt/d/Wasim/Docker/qbfiles")
dfs = readall(".")
dfsp(dfs;save="allbas")


dfs = readall(r"qbasu100.200")
dfsp(dfs)

getnames(dfs)

mx = mall(dfs)
dfp(mx)

fn="qbasbr50.pest2000_clusterbrend-brendpest-brend_pest_v2.wa"
dfp(fn)
lplot(fn)
"qbasbr70.ki2-lake-2008_clusterbrend-output_br70-v2.wa"|>lplot

"qbasu100.2008_clusterufra-qbas.wa"|>lplot

mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
"/mnt/d/Wasim/Docker/clusterufra/"|>cd
dfs = readall(pwd())
dfsp(dfs;save="qbplot")

vgjl("rasters-win.so") #file addressable by alt+klick

([23,33]) .*2

vgr("pyx")

cd("/mnt/d/Wasim/Goldbach/revision/fab150/v4")
dfs=readall(r"qb")
dfsp(dfs)
#recursive_glob_prfx(pwd(),"qb")
route_from_dir(pwd())

v = recursive_glob_prfx(pwd(),"qb")
dfs = map(x->readdf(x),v)
#dfsp(dfs)
typeof(dfs)
function readalldf(pattern::AbstractString)
    """
    reads recursivley and stores to Vector{DataFrame}
    """
    v::Vector{String} = recursive_glob_prfx(pwd(),pattern)
    dfs::Vector{DataFrame} = map(x->readdf(x),v)
    return(dfs)
    pts = map(x->dirname(x),v)
    nm=getnames(dfs)
    # map((x, y) -> joinpath(x,y), v, nm)
    # broadcast((x, y) -> joinpath(x,y), v, nm)
    # #(+).([1, 2, 3], [4, 5, 6])
    # (joinpath).(v, nm)
    # #d = (joinpath).(v, nm)
    #map(x->x*" stored from",nm)
    d = DataFrame(hcat(v,nm), :auto)
    rename!(d,1=>"path",2=>"filename")
    printstyled(d , color=:green    )
    #printstyled("$d\n" , color=:green    )
end

kk=readalldf("qge")
getnames(kk)

# - One way to use map with 2 arguments in Julia is to use the **anonymous function** syntax¹. You can define a function that takes two arguments using the form `(x, y) -> f(x, y)` where `f` is any expression involving `x` and `y`¹. You can then pass this function as the first argument to map, followed by two collections of values for `x` and `y`¹. For example, you can write `map((x, y) -> x + y, [1, 2, 3], [4, 5, 6])` to get `[5, 7, 9]`¹.
# - Another way to use map with 2 arguments in Julia is to use the **dot syntax**². You can apply any 
#function to two collections of values element-wise by adding a dot after the function name². 
#This is equivalent to using map with an anonymous function². For example, you can write 
    #`(+).([1, 2, 3], [4, 5, 6])` to get `[5, 7, 9]`².
# - A third way to use map with 2 arguments in Julia is to use the **broadcast** function³. You can use broadcast to apply any function to two or more collections of values element-wise³. This is similar to using map or dot syntax, but allows more flexibility in handling different shapes and sizes of collections³. For example, you can write `broadcast(+, [1, 2, 3], [4; 5; 6])` to get a 3x3 matrix of sums³.
# Quelle: Unterhaltung mit Bing, 20/04/2023(1) map function with arguments in Julia - Stack Overflow. https://stackoverflow.com/questions/57260592/map-function-with-arguments-in-julia Zugegriffen 20/04/2023.
# (2) map function with multiple arguments to each row of Array in Julia .... https://stackoverflow.com/questions/57263343/map-function-with-multiple-arguments-to-each-row-of-array-in-julia Zugegriffen 20/04/2023.
# (3) Functions · The Julia Language. https://docs.julialang.org/en/v1/manual/functions/ Zugegriffen 20/04/2023.


unique(d)
unique(pts)
split("/").(pts)
x=pts[1]
split(x,"/")

unique(map(x->split(x,"/"),d))
collect(Set(pts))
v
deleteat!(v, occursin(v, v[2:end]))
# - One way to remove duplicates from a Vector{String} in Julia is to use the 
# **unique** function¹. You can pass a vector of strings to unique and 
# it will return a new vector with only the distinct elements¹.
#  For example, you can write `unique(["a", "b", "a", "c"])` to get `["a", "b", "c"]`¹.
# - Another way to remove duplicates from a Vector{String} in Julia is to use 
# the **deleteat!** function². 
# You can use deleteat! to remove elements from a vector at specified indices². 
# You can find the indices of the duplicate elements using the **findin** function²,
#  which returns a vector of indices where a value or a collection of values occurs
#  in another collection². For example, you can write `deleteat!(v, findin(v, v[2:end]))` to remove all duplicates except the first occurrence of each element in vector `v`².
# - A third way to remove duplicates from a Vector{String} in Julia is 
# to use the **Set** type³. You can convert a vector of strings to a set
#  using the Set constructor, which will eliminate any duplicates³. 
#  You can then convert the set back to a vector using the collect function³. 
#     For example, you can write `collect(Set(["a", "b", "a", "c"]))` to get `["a", "b", "c"]`³.
# Quelle: Unterhaltung mit Bing, 20/04/2023(1) Julia: Detect and remove duplicate rows from array?. https://stackoverflow.com/questions/46240706/julia-detect-and-remove-duplicate-rows-from-array Zugegriffen 20/04/2023.
# (2) julia - Remove element from vector while iterating over it - Stack Overflow. https://stackoverflow.com/questions/44874229/remove-element-from-vector-while-iterating-over-it Zugegriffen 20/04/2023.
# (3) Strings · The Julia Language. https://docs.julialang.org/en/v1/manual/strings/ Zugegriffen 20/04/2023.


pt="/mnt/d/Wasim/regio/out/rc200/r2-spin/route_clean"

df = CSV.read(pt,DataFrame,header=false)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
cmd = copy.(eachrow(dout))
##parse.(String, dout[:,1])
string.(dout[:,1])
string.(eachcol(dout))
replace(dout[!,:x1], r"Haf" => "x ")

map(x->replace("pyx "*x,r" " => ""),dout[!,1])
d="qfile"
e = map(x->replace(x*"_qout",r"# " => ""),dout[:,3])
cmd = DataFrame(a="pyx ",b=dout[!,1],c=dout[!,2],d=d,e=e)
cd("/mnt/d/Wasim/regio/out/rc200/r2-spin")
writedf("routecmd",cmd)

#oder inside!
#using PyCall
#python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py
Sys.cpu_info()

a="python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py"
cmd = DataFrame(a=a,b=dout[!,1],c=dout[!,2],d=d,e=e)
cd("/mnt/d/Wasim/regio/out/rc200/r2-spin")
writedf("routecmd",cmd)
onam="routecmd.sh"
xx=" "
repeat(xx, 3)
#Union{Bool, Vector}
nms = ("#/bin/bash","","","","")
Union(1:length(nms),nms)
collect(1:length(nms),nms)

rr=readras(r"qd")
#cntplt(rr)
plot(rr)

#rename(df,("#/bin/bash","","","",""))


# using DelimitedFiles
# hdr = "#/bin/bash"
# open(onam, "w") do io
#     writedlm(io, hdr)
# end
# writedlm(onam, hdr)
#CSV.write(onam,hdr,delim=" ")

hdr = "#/bin/bash\n"
write(onam, hdr)
#CSV.write(onam,hdr,delim=" ",quotestrings=false)
CSV.write(onam,cmd,
header=false,
#quotestrings=false,
#openquotechar=Char(' '),
#closequotechar=Char(' '),
transform = (col, val) -> something(val, missing),
append=true,
delim=" ")  

#cmd
#Base.unlink(onam)

# leider geht das nicht mit qoutechar.
#deshalb in bash:
perl -i -nlwe 'tr/"//d; print if length' routecmd.sh 
. routecmd.sh 
jlkge qout|sort -nr|tee kge.log
#pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/kge_fread.jl"
#include(pt)
cd("/mnt/d/Wasim/regio/out/rc200/r2-spin")
kge_fread()
kge_read(pwd(),"qout")

ds = kge_df("qout")
vgjl("corr")
vgjl("findall")
findall(ds[!,2] .<= -10)

findall(in(ds[!,2]), 2)
ds[(ds.NSE > -10)]

df = ds
vgjl("subset")
DataFrames.subset(ds,:NSE => ByRow(>(-20.)))
ds = filter(:NSE => x -> x > -20, df) 
#ds[ [x in [2,3] for x in ds[!,1]] ,:]

@df ds corrplot(cols(2:3), grid = true)
#ds = ds .* -1
@df ds plot(cols(2:3),xlabel=cols(1))
@df ds plot(cols(2:3),yaxis=:log)

#qqplot(ds)

#DataFrames.eachcol(df)
#(+).(3,DataFrames.eachcol(df))



cd("/mnt/d/Wasim/Tanalys/DEM/Output_Batch_50m")
mkdir("qbfiles")
cd("qbfiles")

cd "/mnt/d/Wasim/Tanalys/DEM/Output_Batch_50m/qbfiles/"
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
# mycommand = "find ../ -type f -size +5k -name 'qbas*'"
# run(mycommand)
dfs = readall(pwd())
dfsp(dfs;save="qbplot")
plotlyjs()
dfsp(dfs)

xdd = map(describe,dfs) #-> transpose

#xm=mall(readall(r"qbasbr50.m"))
xm=(readall(r"qbasbr50.m"))
dfsp(xm)
map(describe,xm)

xm|>first|>dfp

cdb()
nd = rglob("qbas")
dfp(nd[2])
nd[2]|>readdf|>vio


pt=/mnt/d/Wasim/Goldbach/revision/control/
cd $pt
vg "D:/Wasim/Goldbach/" ctl
vg "set \$mainpath" ctl

cd("/mnt/d/Wasim/Goldbach/control/")
vg("tanalys","ctl")


#######streu
"/mnt/d/Wasim/streu/qbfiles"
mkcd qbfiles
find ../ -type f -size +5k -name "qbas*" -not -regex '.*.\(nc\|png\|txt\|sh\)$' -exec \
bash -c 'f="{}"; cp -v "$f" "$(echo $(basename $f))_$((dirname "${f:3}")|tr '/' '-').wa"' \;
#jl
"/mnt/d/Wasim/streu/qbfiles/"|>cd
dfs = readall(pwd())
dfsp(dfs;save="qbplot")

xm=mall(readall(r"v112"))
dfpjs(xm)

xdd = map(describe,dfs) #-> transpose

#mean baseflow 1.394585 Q (m3/s)
ve = rglob(r"v112")
dfpjs("./qbasstr.v112017_out-v11.wa")
dfpjs(ve[2])

dfs=readall(r"v112")
getnames(dfs)

s="/mnt/d/Wasim/streu/out/v9/qbasstr.v92017"
dfpjs(s)
lplot(s)
"/mnt/d/Wasim/streu/out/v9/"|>cd
rplot(r"tsoi")
dfp(r"ts_lo")


#df=readdf("/mnt/d/Wasim/streu/meteo/v2/p_1970.txt")
"/mnt/d/Wasim/streu/out/v15"|>cd
dfs=readall(pwd())
dfsplog(dfs)

filterplot("qg",dfs)
filterplot("qbas",dfs)
filterdf("qbas",dfs)|>lplot
filterdf("sb",dfs)|>lplot

filterdf("sb1",dfs)|>dfl
filterdf("rge",dfs)|>dfl!

x0 = filterdf("rge",dfs)
yrmean(x0)

#@df df Plots.plot(:month,cols(propertynames(nd,Not(:month))))

df = readdf(r"rge")|>monsum
df = readdf(r"rge")|>monmean
str = [ @sprintf("%02i", x) for x in (df.month) ];
s = Symbol.(filter(x->!occursin("month",x),names(df)))
@df df StatsPlots.bar(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=true)
#@df df StatsPlots.violin(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=true)
readdf(r"rge")|>vio
barp(readdf(r"rge"))

"/mnt/d/Wasim/Goldbach/out/sglo_v44"|>cd
ve = rglob(r"rgexggb")
ve = filter(x->!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",x),ve)
dfs=readall(ve)
dfsp(dfs)

plotlyjs()
using PlotlyJS
fig = dfsp(dfs);
# PlotlyJS.relayout!(fig, title = "Time Series with Custom Date-Time Format",
#     xaxis_tickformat = "%d %B (%a)<br>%Y")

# fig

dfsjs(dfs)

vgr("pandas")
("tst.wa")|>vgjl
("tst.wa")|>vgr
("E+34")|>vgjl
(".kml")|>vgjl
(".kml")|>vgpy
("E+34")|>vgr
("tst.wa")|>vgr
("tst.wa")|>vgpy

Sys.cpu_info()

fn="/mnt/d/ClimateExplorer/radiation/glob_10.6293___49.7787_30.txt"# import datatable as dt
# df = dt.fread(file=fn,skip_to_line=3,
#               na_strings=["3e+33","3e+3333333","3e+34","3e+35","-9999","-9999.0","-999","-999.0"],
#             header=True,nthreads=12,fill=True,skip_blank_lines=True)
vgjl("clime")
ms=["3e+33","3.0e33","3e+34","3e+35","-9999","-9999.0","-999","-999.0"]
df = CSV.read(fn,
    DataFrame,
    header=false,
    missingstring=ms,
    ntasks=4, limit=typemax(Int), ignorerepeated=true,
    skipto = 4,delim=" ",normalizenames=true)
    
rename!(df,1=>"date")
df.date = map(x->Date(string.(x),"yyyymmdd"),df.date)

#ds = filter(x->map(y => x -> x > 0, propertynames(df)[2:end]),df) 

ds = df[findall(df[!,2] .< 1000),:]
#replace!(x -> isodd(x) ? 2x : x, [1, 2, 3, 4])

replace!(x -> isnan(x) ? 2x : x,df[!,i])

for i in propertynames(df)[2:end]
    #dx = DataFrames.subset(df,i => ByRow(<(1000.)))
    df[findall(df[!,i] .< 1000),:]
end

ds = df[findall(df[!,2] .< 1000),:]
ds = ds[findall(ds[!,10] .< 500),:]
dfp(ds)


#perl -i -pe 's/0.3000000E+34/-9999/g' glob_10.6293___49.7787_30.txt

julia
x="glob_10.6293___49.7787_30.txt"
file = open(x, "r")
# Read in the contents of the file
contents = readlines(file)
# Close the file
close(file)

# Replace the string "0.3000000E+34" with "-9999"
for i in 1:length(contents)
    contents[i] = replace(contents[i], "0.3000000E+34" => "-9999")
end
# Open the file for writing
file = open("out.txt", "w")
# Write the modified contents back to the file
for line in contents
    write(file, line*"\n")
end
# Close the file
close(file)
#rm("out.txt")

file=open("out.txt", "r")
cs = readlines(file)
close(file)
# cp(x,"tmp.txt")
# tx="tmp.txt"
# for line in eachline(tx)
#     tmp = read(line)
#     replace!(tmp, "0.3000000E+34" => "-9999")
#     open("my_file.txt", "w") do io
#         write(io,tmp)
#     end
# end

fn="/mnt/d/ClimateExplorer/radiation/out.txt"
ms=["3e+35","-9999","-9999.0","-999","-999.0"]
df = CSV.read(fn,
    DataFrame,
    header=false,
    missingstring=ms,
    ntasks=8, 
    limit=typemax(Int), ignorerepeated=true,
    skipto = 4,
    delim=" ",normalizenames=true)

#first(df,10)
rename!(df,1=>"date")
df.date = Date.(string.(df.date),"yyyymmdd")
last(df,10)
cd("/mnt/d/ClimateExplorer/radiation/")
ll()
using GeoJSON
json_str = read("radcrds.json", String)
# Parse the string as a GeoJSON object
first(json_str,20)
geojson_obj = GeoJSON.parse(json_str)
# Print the object
println(geojson_obj)

using GeoDataFrames
s = raw"D:\Wasim\main\tdx\whsed.shp"
#p=GeoDataFrames.read(s)
using Shapefile
riv = Shapefile.Handle(s)
using Plots
riv.shapes[1:end]|>Plots.plot
pyplot()
riv.shapes[5:end]|>plot

using DataFrames
using GeoDataFrames
using PyPlot
#PyPlot.plot(riv.shapes[5])

#PyCall is currently configured to use the Python version at:
#C:\Users\chs72fw\AppData\Local\miniforge3\python.exe


using GeoDataFrames

pts=GeoDataFrames.read("climexp2083545.kml")
#Plots.plot(pts)
pts|>names
#ab = pts[!,Cols(2,1)]
ab = pts[!,Cols(end)]

fn="radcrds2.json"
pts=GeoDataFrames.read(fn)
pts|>names
Plots.plot(pts)
#Shapefile.read
using Shapefile
Shapefile.write("points.shp", pts)
xy = Shapefile.Handle("points.shp")


fn=raw"D:\Wasim\regio\out\rc200\x24\ei__rcm_Layer_1.2017.nc"
xrp(fn)



#https://docs.juliaplots.org/latest/generated/unitfulext_plots/
using Unitful,Plots
group = rand(map((i->"group $(i)"), 1:4), 100)
plot(rand(100)*u"km", layout=@layout([a b;c]), group=group, 
linetype=[:bar :scatter :steppre], linecolor=:match)


plot(rand(100)*u"mm", layout=@layout([a ;b c]), group=group, 
linetype=[:bar :scatter :steppre], linecolor=:match)


"PkgCleanup"|>vgjl
PkgCleanup.artifacts()
PkgCleanup.manifests()

using GeoDataFrames
x="/mnt/d/ClimateExplorer/precipitation/precipitation_climexp2509005.kml"
pts=GeoDataFrames.read(x)

"kml"|>vgpy

rplot("/mnt/d/Wasim/regio/out/rc200/v4/thawrcm.2012.nc")

"Explor"|>vgjl


"/mnt/d/Wasim/regio/out/rc200/v4"|>cd
rplot(r"qd")
using Grep
# Grep.grep("qd",readdir())
# rglob("qd")
# glob("qd")
@time glob("qd")
@time rglob("qd")
@time Grep.grep("qd",readdir()) #fastest.

Grep.grep|>methods


@time Grep.grep(r"^q",readdir())
@time Grep.grep(regand("q","nc"),readdir()) #eqally fast

plotlyjs()
Grep.grep(regand("q","nc"),readdir())[1]|>rplot
rplot(r"qsmer")

ll(x="qges")
ll(x="qout")

ds = kge_df("_qout")
vgjl("sort")
vgjl("order")

sort(ds[!, :KGE],rev=true)
sort(ds, order(:KGE, rev=true))
vgjl("cor")
@df ds corrplot(cols(2:3), grid = true)



"/mnt/d/Wasim/Goldbach/out2/v16/v16-2023"|>cd

lplot(r"rgex")
dfp(r"rgex")
bardfm(r"rgex")
lplot(r"qout")
dfp(r"radiat")

dfp(r"rgex")
dfp(r"qbas")
rplot(r"wi")
rp(r"wind")
rpm(r"wind";msk=.7)
rpm(r"wind";msk=1.1)
rpm(r"wind";msk=.7,gt=true)
rpm(r"tem";msk=10.1,gt=true)


raster tdiff
# pot = readras(r"etp")
# real = readras(r"etr")
pot  = rglob(r"etp")
real = rglob(r"etr")
filter(x->occursin(regand("Layer","nc"),x),rglob(r"etr"))
pot = filter(x->occursin(regand("Layer","nc"),x),pot)|>last|>readras
real= filter(x->occursin(regand("Layer","nc"),x),real)|>last|>readras

td = pot-real
rplot(td;ex=1)

vgjl("rebuild")
td=Rasters.rebuild(td;missingval=-9999)
##for overwriting
write("tdiff.nc",td;force=true)

#########tdiff windows test
pot  = rglob(r"etp")
real = rglob(r"etr")
pot = filter(x->occursin(regand("Layer","nc"),x),pot)|>last|>readras
real= filter(x->occursin(regand("Layer","nc"),x),real)|>last|>readras
td = pot-real
#rplot(td;ex=1)
#td=Rasters.rebuild(td;missingval=-9999) #not really nec
##for overwriting force
plot(td)
write("tdiff.nc",td;force=true)


"/mnt/d/Wasim/regio/out/lowres/c9"|>cd
dfp(r"rgex")
qgk()


"/mnt/d/Wasim/Goldbach/revision/fab150/pestpp/v5/pestout"|>cd

dfpall(r"rgex")
rpm(r"rad_fab";msk=1.0)

v = r"rad"|>rglob
rplot(v[3])

v[broadcast(x->endswith(x,"nc"),v)]|>last|>rplot

nconly("tem")|>only|>rplot
nconly("tsoi")|>only|>rplot

dfp(r"avg")
r"ts_"|>rglob
r"soil"|>rglob
r"soil"|>rglob|>first|>dfp

"soil"|>rglob


"/mnt/d/Wasim/Goldbach/revision/fab150/pestpp/v5/pestout/so-files"|>cd

pr=readdf(r"prec")
dfs=loadalldfs(r"prec")

(r"prec")|>glob

dfs=readall(r"prec")
dfs=alldf(r"prec")
#ja, weil halt txtfiles dabei sind.

oo = mall(dfs)
vio(pr)
vio(oo)
dfp(oo)


"/mnt/d/Wasim/regio/out/wrn_c0"|>cd
precip=readdf(r"pre")
discharge=readdf(r"qges")
# Erstellen Sie ein Diagramm mit zwei Achsen

plot(precip[!,1], discharge[!,1], xlabel="Niederschlag (mm/timestep)", 
ylabel="Flussabfluss (mm/timestep)", label="Wern", 
seriestype=:scatter)

dv=hcat(precip[!,1], discharge[!,1])
qqplot(precip[!,1], discharge[!,1], qqline = :fit)
qpl(DataFrame(dv,:auto))

twinx!()
plot!(grd.LATC, grd.LONC, grd.bathymetry', clims=(-6000, 0), color=:blues, colorbar=false, label="", legend=:topleft)

# cut -f0-4 Ettleben |h
# cut -f0-4 Ettleben |hd
# cut -f1-4 Ettleben |hd
# cut -f1-4 Ettleben > 0.sub
# awk '{row[FNR] = (FNR==NR? $NF: row[FNR] FS $NF)}END {for (i=1;i<=length(row);i++) print row[i]}' *.sub > whdrs

a=readdf(r"Arns")
e=readdf("Ettleben")
s=readdf(r"Sachs")
df = reduce((left, right) -> 
      innerjoin(left, right, on = :date,makeunique=true), 
      [a,e,s])
##reorder columns
df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
###rename!(df,sort(names(df))) ##wrong!
qpl(df)
writewa("wernspec",df)


"/mnt/d/Wasim/regio/out/wrn_c1"|>cd

"pytd"|>glob|>second|>vio
"pytd"|>glob|>second|>vibx
plotlyjs()
"ts"|>glob|>first|>rp
"va"|>glob
"tdiff"|>vgjl

#########tdiff 
function tdifnc()
    """
    path="."
    read non-recursivley and plots tdiff
    glob<->rglob
    stores tdiff.nc

    """
    pot = filter(x->occursin(regand("Layer","nc"),x),glob(r"etp"))|>last|>readras
    real= filter(x->occursin(regand("Layer","nc"),x),glob(r"etr"))|>last|>readras
    td = pot-real
    nm = filter(x->occursin(regand("Layer","nc"),x),glob(r"etp"))|>last
    ti = split(nm,".")|>x->x[end-1]
    plot(td;
    title = "Tdiff[mm] of "*ti,
    c=cgrad(:matter))
    ##for overwriting force
    #write("tdiff.nc",td;force=true)
end



ll(;reg=true) |> third|>rp
ll(;reg=true) |>x->x[end-4]|>facets

ll(;reg=true) |>y -> filter(x->occursin(r"temp",x),y)

ll(;reg=true) |>y -> filter(x->occursin(r"temp",x),y)|>only|>facets





cd("old")
kd = fz()

# sort!(kd, [order(:total,rev = true)])

# a = kd[!,end]|>sum
# b = kd[!,end-1]|>sum

# m=Dict(:a=>a,:b=>b)
# hcat(DataFrame(m),bla=2342)
# jdd()

get_folder_size(".")/1000^2

dfp(r"rge")

cd("/mnt/d/Wasim/regio/out/wrn_c1")
"sobw.wrn_c1"|>dfp
r"global"|>dfp

tdifnc()

r"vapo"|>facets
r"so_vapor_pressure"|>dfp

"/mnt/d/Wasim/regio/out/rc200/v5/vaporcm.2012.nc"|>rplot
rplot!(r"wind")

cd("/mnt/d/temp/saale/out_smf200/cl_cdx_v4")
rplot(r"vapo")
"/mnt/d/wasim/tanalys/dem/input_v2/meteo/ts-0/"|>cd
ll(;reg=true,x="pr")
fz()
df = dfread("pressure_ce.txt")

x="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/cloudcover_ce.txt"
df=waread(x)
dfx = describe(df)|>permutedims
#map(String,dfx[1,:])|>DataFrame
#map(collect,dfx[1,:])
row = dfx[1, :] # first row of df
vec = collect(row) # convert to vector
rename!(dfx,vec)

wu = df[!,Cols(3,:date)]
metadata!(wu, "filename", names(wu)[1], style=:note);
dfl(wu|>dropmissing)
dropmissing!(wu)

sowu=waread("/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/sonnenstunden_1970.txt")
names(sowu)
xwu = sowu[!,Cols(r"Wuerz",:date)]
bo = innerjoin(xwu,wu,on=:date)
"reorder"|>vgjl
"glm"|>vgjl

reorder!(bo)

using DataFrames, StatsPlots, GLM, CSV, Dates
#df = rename(nd,1=>"modis",2=>"wasim")
#model = lm(@formula(modis ~ wasim), df)
df = bo
model = lm(@formula(Wuerzburg ~ WURZBURG), df)
gr()
data = df[!,Not(Cols(1,:date))]
rename!(data,1=>"x",2=>"y")
r = lm(@formula(x ~ y), data)
pred = predict(r, data, interval = :confidence, level = 0.95)
pd = @df data Plots.scatter(:x, :y, leg = false)
# sort data on x
pred_s = pred[sortperm(data[!,:x]), : ]
x_s = sort(data[!, :x])
Plots.annotate!([(7,7,text("hey", 14, :left, :top, :green))])
Plots.annotate!([(10,3,"some R²...")])
Plots.plot!(pd, x_s, pred_s.prediction, linewidth = 2,
        ribbon = (pred_s.prediction .- pred_s.lower, pred_s.upper .- pred_s.prediction))

using LinearFitXYerrors, DataFrames
linearfitxy(data.x, data.y, isplot=true, ratio=:auto)
        

###########with preci ###########################
x="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/pre_ce.txt"
y="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/preci_1970.txt"
dfce=waread(x)
dfdwd=waread(y)
wu = dfce[!,Cols(r"WURZ",:date)]
wu = wu[!,Cols(1,:date)]
metadata!(wu, "filename", names(wu)[1], style=:note);
dropmissing!(wu)
xwu = dfdwd[!,Cols(r"Wuerz",:date)]
xwu = xwu[!,Cols(2,:date)]
bo = innerjoin(xwu,wu,on=:date)
bo = reorder!(bo)
data = bo[!,Not(Cols(:date))]
rename!(data,1=>"x",2=>"y")
using LinearFitXYerrors, DataFrames
linearfitxy(data.x, data.y, isplot=true, ratio=:auto)


########### Radiation ###########################
x="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/rad_ce.txt"
y="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/radiation_2007.txt"
dfce=waread(x)
dfdwd=waread(y)
wu = dfce[!,Cols(r"WURZ",:date)]
metadata!(wu, "filename", names(wu)[1], style=:note);
dropmissing!(wu)
dfp(wu)
dfdwd[1,:]
xwu = dfdwd[!,Cols(r"Rander",:date)]
bo = innerjoin(xwu,wu,on=:date)
bo = reorder!(bo)
bo.Randersacker = bo.Randersacker/24 ##ooh

##subset by DATE without ByRow !!!
filter(:date => x -> month(x) == 12, bo) #
data = filter(:date => x -> Dates.year(x) >= 2012, bo)
data = data[!,Not(Cols(:date))]
rename!(data,1=>"x",2=>"y")
#using LinearFitXYerrors, DataFrames
#gr()
#linearfitxy(data.x, data.y, isplot=true, ratio=:auto)
dropmissing!(data)
r = lm(@formula(y ~ x), data)
pred = predict(r, data, interval = :confidence, level = 0.95)
pd = @df data Plots.scatter(:x, :y, leg = false)
# sort data on x
pred_s = pred[sortperm(data[!,:x]), : ]
x_s = sort(data[!, :x])
Plots.plot!(pd, x_s, pred_s.prediction, linewidth = 2,
        ribbon = (pred_s.prediction .- pred_s.lower, pred_s.upper .- pred_s.prediction))

Plots.annotate!([(7*7*7*2,5,text("W/sqm", 12, :left, :top, :green))])
#Plots.annotate!([(1024,300,"some R²...")])

plotlyjs()
dfpjs(bo)


#################
"/mnt/d/Wasim/regio/out/wrn_c1"|>cd
r"tsoil"|>facets
r"vapo"|>facets
Plots.savefig("sin.png")

methods(Plots.savefig)

pyplot() # activate PyPlot as backend
using PyPlot
#findfont: Font family 'Computer Modern' not found.
r"vapo"|>rp
PyPlot.savefig("sin_pyplot.png", dpi=800)


glob("qout")
d = kge_df("qout.py")


"rad"|>ct|>first|>dfp

"/mnt/d/Wasim/streu/out/coarse/v4"|>cd

rp(r"sb")
rp(r"gwst")
facets("win")
facets("soil")
("nc")|>glob

r"rad"|>dfp
r"rgex"|>dfp
r"qb"|>dfp
r"qges"|>dfp
r"prec"|>bardfm
r"sb1"|>bardfm
r"sb1"|>dfp
r"prec"|>bardf
r"qgko"|>dfp
r"vapo"|>glob|>first|>bardf
r"vapo"|>glob|>first|>dfp
r"vapo"|>dfp

r"ts_"|>glob|>only|>dfp

"humistr"|>glob|>second|>rp

"humistr"|>glob|>second|>rp



#("temp",  "prec",  "wind",  "humi",  "vapo",  "rad_")


ncs = filter(x->occursin(r"temp|prec|wind|humi|vapo|rad_",x),readdir())
#ncs = map(x->nconly(x),ncs)
ncs = filter(x->endswith(x,".nc"),ncs)

for x in ncs
    describe(Raster(x))
end
run(`cdo -V`)
for file in ncs
    cmd = `cdo infon $file`
    run(cmd)
end

cmd=`cdi $ncs[2]`
run(cmd)
##but!
for file in ncs
    cmd = `cdi $file`
    run(cmd)
end



"/mnt/d/Wasim/Tanalys/DEM/brendpest/out_v7/spin"|>cd
r"brend"|>glob
r"brend"|>dfp
r"qbas"|>dfp
"brend-spin"|>ftp

tdifnc()

theplot("brend-spin_qout")
bardf("brend-spin_qout")

vg("control file ","xml") 

function ctlx()
    """
    xmlgrep #snippet::AbstractString
    """
    # [d for d if isdir(d) in readdir()]
    cwd = pwd()
    nwd_list = []
    xmls = []
    ext=".xml"
    for (root, dirs, files) in walkdir(cwd)
        for dir in dirs
            if isdir(dir)
                push!(nwd_list,joinpath(cwd,dir))
            end
        end
        nwd_list = length(nwd_list)==0 ? pwd() : nwd_list
        for f in files
            if (isfile(f)) & (endswith(f, ext))
                push!(xmls,joinpath(cwd,f))
            end
        end
    end

    #nwd_list = length(nwd_list)==0 ? pwd() : nwd_list
    snippet="control file "
    xmls


    for nwd in nwd_list
        #cd(nwd)
        #println("greps from *ctl from $nwd...")
        files = filter(f -> endswith(f, ext), readdir())
        for file in files
            counter = 0 # Zähler initialisieren
            open(file) do f
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line, snippet)
                        printstyled("$counter:$nwd", color=:light_red) 
                        printstyled("$file:\t", color=:light_magenta, underline=true, blink=false, bold=true) 
                        printstyled("$line\n", color=:green, bold=true) 
                    end
                end
            end
        end
    end
    #cd(owd)
end

ctlx()

"/mnt/d/Wasim/Tanalys/DEM/brendpest/out_v7/spin"|>cd



cdb()



function ctlx(ext=".xml",pattern="compiling symbols in control file")
    # use glob to find all xml files recursively in the directory
    #files = Glob.glob("*.xml", directory; recursive=true)
    cwd = pwd()
    # nwd_list = []
    xmls = []
    for (root, dirs, files) in walkdir(cwd)
        # for dir in dirs
        #     if isdir(dir)
        #         push!(nwd_list,joinpath(cwd,dir))
        #     end
        # end
        # nwd_list = length(nwd_list)==0 ? pwd() : nwd_list
        for f in files
            if (isfile(f)) & (endswith(f, ext))
                #push!(xmls,joinpath(cwd,f))
                push!(xmls,relpath(f))
            end
        end
    end
    # loop through each file
    l = []
    for file in xmls
      # open the file for reading
      open(file, "r") do f
        # loop through each line in the file
        for line in eachline(f)
          # use re.search to find the pattern in the line
          match = occursin(pattern, line)
          # if there is a match
          if match
            # print the file name and the line
            #println(file * ": " * line)
                # remove substrings
            sub = "\t<line level=\"4\" description= \"compiling symbols in control file "
            line = replace.(line, sub => "")
            sub2 = "\"/>\n"
            line = replace.(line, sub2 => "")
            push!(l, file * ": " * line)
          end
        end
      end
    end
    return l
  end



ctlx()



r="D:/temp/saale/in_mf/smf2/uf.tif"|>readras
r|>plot
r="D:/temp/saale/output/v6/tsoilsmf180_stack.2017.nc"|>readras
r[t=3]|>plot
"surfa"|>vgjl
using GeoArrays
x = r[t=3]
rp3(x)


vgjl("mos")
vgpy("mos")


raw"D:\Wasim\ecad"|>cd
v=glob("QQ")
##eacd reader
CSV.read(v[1],DataFrame;
    missingstring=["-9999"],
    header=20,
    skipto=21,
    types = Dict(3=>Date,4=>Float64),
    select=[3,4],
#    ignorerepeated = true,
    validate = true,
    stripwhitespace=true,
    dateformat = "yyyymmdd",
    limit = typemax(Int) )

function eacd_read(x::String)
    df = CSV.read(x,DataFrame;
            missingstring=["-9999"],
            header=20,
            skipto=21,
            types = Dict(3=>Date,4=>Float64),
            select=[3,4],
            validate = true,
            stripwhitespace=true,
            dateformat = "yyyymmdd",
            limit = typemax(Int) )
    nm = CSV.read(x,Tuple; skipto=17,limit=1,stripwhitespace=true)|>last|>only
    DataFrames.metadata!(df, "filename", nm, style=:note);
    nx = replace(nm,"("=>"",")"=>"",": "=>"_"," "=>"_")
    rename!(df,1=>"date")
    rename!(df,2=>nx)
    #dropmissing!(df)
    return df
end

dfs = broadcast(eacd_read,v)
dfs = filter(x->nrow(x)>0,dfs)
dfp(dfs[5])

dsub = map(u->filter(x->(year.(x.date)>1969),u),dfs)
mx = mall(dsub)
wawrite(mx,"qq_2023.wa")

#dfm("qq_2011-2023.wa")
od = readdf("qq_2011-2023.wa")
od = readdf("qq_2023.wa")
wa.dfm(od;ann=false)
wa.dfm(od;fun=yrsum,ann=false)
wa.dfm(od;ann=false)


lk=raw"D:\remo\cordex\eobs\rr_ens_spread_crop.nc"
r=Raster(lk)
describe(r)
r
r[Ti=3]|>plot

dx = dfr(raw"D:\remo\cordex\eobs\rr_eobs.sommerach")
wa.dfm(dx;fun=yrsum,mode=:bar)

lk=raw"D:\remo\cordex\eobs\qq_ens_spread_crop.nc"
r=Raster(lk)
r[Ti=3]|>plot

qd = wa.ncdf(r)
baryrsum(qd)
wa.dfm(qd;fun=yrmean,mode=:steppost,ann=false)


r = readras(r"tso")
rsu = r[t=2:10]
plot(rsu;axis=false,legend=false)


using PyPlot,CSV,DataFrames
n=raw"D:\Wasim\regio\out\rc200\x24\etr_rcm.x24_Layer_1.2000"
df = readdf(n)
pygui(true)
pyplot_df(df)


using PyPlot,DataFrames

begin 
    n=raw"D:\Wasim\regio\out\rc200\x24\qbasrcm.x24.2017"
    df = wread(n) #dlm
    pygui(true)
    pyplot_df(df)
end


#dateranges in julia:
collect(Date(2000):Month(1):Date(2004,2,1))

function dtrange(start, stop, step)
    """
    dtrange(Date(2004),Date(2005),Month(4))
    """
    collect(start:step:stop)
end

dtrange(Date(2004),Date(2005),Month(1))|>Plots.plot

dtrange(Date(2004),Date(2005),Dates.Day(10))
dtrange(Dates.Second(1),Dates.Second(1000),Dates.Second(10))
a = dtrange(cos(3),sin(3),0.001)
b = dtrange(Dates.Second(1),Dates.Second(size(a,1)),Dates.Second(1))
Plots.plot(b,a)
Plots.plot(cos.(a))
Plots.plot(cos.(1:1234:100000))




######concatenated map funcs ##################

vector = map(fsize,map(basename,nconly("stack")))
#vector = map(fsize,map(x->replace(x,".nc"=>""),nconly("stack")))

# DataFrame erstellen
df = DataFrame(Name = String[], Size = Float64[])
# Vektor in DataFrame umwandeln
for i in 1:length(vector)
    for j in 1:length(vector[i])
        push!(df, (vector[i][j][1], vector[i][j][2]))
    end
end
pretty_table(df;header=names(df))

for x in eachrow(df)
    println(replace(x[1],".nc"=>""))   
end
#same
map(x->replace(x[1],".nc"=>""), eachrow(df))

###inplace errors on DFs.
for x in eachrow(df)
    (replace!(x[1],".nc"=>""))   
end

df[!,1] .= map(x->replace(x[1],".nc"=>""), eachrow(df))
df



wa.dfp("D:/remo/qm/corgrids/wind_rcm85.wa2")
df = wa.dfr("D:/remo/qm/corgrids/pre_rcm85.wa2")
wa.dfm(df;fun=yrsum,mode=:bar,ann=false)
dfp(df[!,Cols(r"w|date")])
@edit wa.tline(df[!,Cols(r"w|date")],:date)


df = wa.dfr("D:/remo/qm/corgrids/tas_rcm85.wa2")

dfp(df[!,Cols(r"W|date")])
wa.tline!(df[!,Cols(r"W|date")])


se = df[!,Cols(r"W|date")]|>yrsum
se = df|>yrsum
dfm(se;mode=:box,fun=false,ann=false,leg=false)
wa.tline!(se;date_col=:year)

@rsubset df :date => x -> Dates.year(x) == 2017

filter(:date => x -> year(x) > 2050, df)
# @transform(df, :date = year.(:date)) no
@rsubset(df, year(:date) > 2050) #yes

fu = @rsubset se :year > 2050 


println(names(df))



#docstings have to be obove the function
module xa

    """
    Calculates the sum of two numbers.
    
    Args:
        a (Int): The first number.
        b (Int): The second number.
    
    Returns:
        Int: The sum of the two numbers.
    """
function sumx(a::Int, b::Int)
    return a + b
end 
end

@doc xa.sumx

using Base.Docs
@doc Base.Docs.md"xa.sumx"

x = @chain [1, 2, 3] begin
    filter(!=(2), _)
    sqrt.(_)
    sum
end


n="D:/Wasim/Testdata/out2/tsoilcols500.grd.nc"
pyjl.xrfacets(n)

df = dfr(raw"D:\Wasim\Testdata\Testdata2\Output\qout")
df = select(df,[2,1,3])
ftp(df)

# R"""
# pak::pkg_install("https://github.com/hzambran/hydroGOF/releases/tag/v0.4-0") 
# """
"D:/temp/saale/out_smf200/v6/loc/"|>cd
vgjl("@layout")