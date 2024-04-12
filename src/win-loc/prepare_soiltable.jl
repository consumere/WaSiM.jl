cd(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8")
# time: 2023-11-30 15:54:22 W. Europe Standard Time
	df=CSV.read("Leitprofile_BFID.txt",DataFrame)
# time: 2023-11-30 15:55:06 W. Europe Standard Time
	sel=select(df,[1,2,6])
#
# using Shapefile
# s=Shapefile.Handle("D:/Bodendaten/buek200_2020/merged_buek200.shp")
# using GeoDataFrames
# gd = read("D:/Bodendaten/buek200_2020/merged_buek200.shp")
# println(names(gd))

using RCall
@rimport LWFBrook90R as lwf
v = lwf.hydpar_wessolek_tab("Su4")
#R to julia
M=rcopy(v)
@rimport utils as ut
ut.help("hydpar_wessolek_tab")
ut.help("dataframe")
@rimport terra as tr
gd = tr.vect("D:/Bodendaten/buek200_2020/merged_buek200.shp")
println(gd.names)


	lwf
# time: 2023-11-30 15:55:30 W. Europe Standard Time
	bst = sel.BOART
# time: 2023-11-30 15:55:48 W. Europe Standard Time
	@rput bst
# time: 2023-11-30 15:56:07 W. Europe Standard Time
	ot = lwf.hydpar_wessolek_tab(bst)
# time: 2023-11-30 15:56:28 W. Europe Standard Time
	odf = rcopy(ot)
# time: 2023-11-30 15:56:56 W. Europe Standard Time
	odd = hcat(odf, bst)
# time: 2023-11-30 15:57:06 W. Europe Standard Time
	odd = hcat(odf, boart=bst)
# time: 2023-11-30 15:57:10 W. Europe Standard Time
	odd
# time: 2023-11-30 15:57:38 W. Europe Standard Time
	rename!(odd,8=>"bart")
# time: 2023-11-30 15:58:03 W. Europe Standard Time
	lwf.hydpar_wessolek_tab("Sl2")
# time: 2023-11-30 15:58:14 W. Europe Standard Time
	using DataFramesMeta
# time: 2023-11-30 15:58:44 W. Europe Standard Time
	@rsubset odd :8=="Su2"

v=["Lu",  "Lu",  "Lt2",  "Lu"]
ot = lwf.hydpar_wessolek_tab(v)
ot = rcopy(ot)
#[mm/d] to [m/s]  (meters per second)
#aus dem R script
fact = 1/(1000*24*60*60)
ot[!,:ksat] .= ot[!,:ksat]*fact
# aus der Control File
f = [9.56944446588e-06 9.56944446588e-06 7.237384275471e-06 9.56944446588e-06 9.56944446588e-06]
scatter(ot[!,:ksat],vec(f))
cv = vec(f)[1:4]
diff = ot[!,:ksat] .- cv
plot(diff,label="diff",seriestype=:bar)
plot!(ot[!,:ksat],label="ksat",seriestype=:scatter)
plot!(cv,label="cv",seriestype=:scatter)


# time: 2023-11-30 16:07:11 W. Europe Standard Time
	writedf(odd,"wessolek-ptfs.txt")
# time: 2023-11-30 16:07:21 W. Europe Standard Time
	odd
# time: 2023-11-30 16:07:26 W. Europe Standard Time
	df
# time: 2023-11-30 16:08:16 W. Europe Standard Time
	all = innerjoin(odd, df, on=[:BOART, :bart])
# time: 2023-11-30 16:08:26 W. Europe Standard Time
# mode: help
	innerjoin
# time: 2023-11-30 16:09:02 W. Europe Standard Time
	rename!(df, :BOART=>:bart)
# time: 2023-11-30 16:09:11 W. Europe Standard Time
	all = innerjoin(odd, df, on=:bart)
# time: 2023-11-30 16:09:42 W. Europe Standard Time
	writedf(all, "soildata_wessolek.txt")


df = CSV.read("soildata_wessolek.txt", DataFrame)
tg = CSV.read("tblZuordnungGLE_BLE.csv", DataFrame)

println(names(df))
ot = innerjoin(df, tg, on=:TKLE_NR,makeunique=true)

writedf(ot, "soildata_wessolek_genid.txt")

## änderung und nochmal
"tblZuordnungGLE_BLE.csv" :TKLE_NR -> "tblBlattlegendeneinheit.csv"
## "tblZuordnungGLE_BLE.csv" :GEN_ID -> "tblProfile.csv" :BF_ID 1 -> ∞ "tblHorizonte.csv"
"tblZuordnungGLE_BLE.csv" :GEN_ID -> "tblProfile.csv" :GEN_ID "tblHorizonte.csv"

zuo = CSV.read("tblZuordnungGLE_BLE.csv", DataFrame)
leg = CSV.read("tblBlattlegendeneinheit.csv", DataFrame)
pro = CSV.read("tblProfile.csv", DataFrame)
hor = CSV.read("tblHorizonte.csv", DataFrame)
select!(hor,1:10)

using DataFramesMeta
@rsubset hor :BF_ID==3212
@rsubset hor :BF_ID==10479
@rsubset hor :BF_ID==4340
@rsubset hor :BF_ID==1008
@rsubset hor :BF_ID==1010
@rsubset hor :BF_ID==10461

ot = lwf.hydpar_wessolek_tab(["O", "Lu", "Lu"])|>rcopy
ot[!,:ksat] .= ot[!,:ksat]*fact


#select!(hor,Not(:SW_HORIZ))
@rsubset hor :BOART=="Tu4"
#map(x->replace(x,missing => "Tu4"),hor.BOART)
map(x->replace(x,missing => "Tu4"),hor.BOART)
v = ismissing.(hor.BOART)
#count([x==true for x in v])
last(hor)
@rsubset hor ismissing(:9)
@rsubset hor ismissing(:BOART)
@rsubset hor !ismissing(:GROBBOD_F)
@rsubset hor ismissing(:BOART)
@rsubset hor  :UTIEF .>= 21
for i in eachrow(hor)
    if any(ismissing, i)
        i.BOART = "Tu4"
    end
end
#@rsubset hor ismissing(:BOART)
@rsubset hor :BOART=="Tu4" #stark schluffiger Ton

#load only Saale bfids.
setup()
ras = readras("D:/Wasim/regio/rcm200/v12/rcm.art-bfid")
plot(ras)
using GeoDataFrames
using GeoInterface
const gdf = GeoDataFrames
n="D:/Wasim/regio/rcm200/v12/catchment.shp"
g = gdf.read(n)
xm = mask_trim(ras,g.geometry)
plot(xm)

bfids = round.(Int64, unique(vec(xm.data))) #69 without 0  

#bfids = parse.(Int64,unique(vec(ras.data))) #error.
#bfids = round.(Int64, unique(vec(ras.data))) #94
subset_hor = hor[in.(hor.BF_ID, Ref(bfids)), :]
subset_hor = subset_hor[!, Not(10:end)]
println(names(subset_hor))
any(ismissing.(subset_hor.BOART))

df2 = innerjoin(subset_hor,pro, on=:BF_ID, makeunique=true)

#Ref(unique(df2.GEN_ID))
df2.GEN_ID in df2.GEN_ID_1
unique(df2[!,Cols(r"^GEN")])
unique(df2[!,Cols(r"^BOD")])
df3 = leftjoin(df2,zuo,on=:GEN_ID, makeunique=true)
out = leftjoin(df3,leg,on=:TKLE_NR, makeunique=true)

unique(out[!,Cols(r"^GEN")])
soils = unique(df2[!,Cols(r"^BOD|^BF_ID")])
sort!(soils,:4)
pwd()
open("s200-table.txt", "w") do io
    pretty_table(io, unique(soils,:3); 
    header=uppercasefirst.(names(soils)), backend = Val(:text))
end

bfids == sort(unique(out.BF_ID))

#m ids
mi = sort(unique(out.BF_ID))
mr = sort(bfids[2:end])           #exclude 0
mr = mr[2:end] #now exclude -9999
mi == mr          ##fits!

select!(out,Not(Cols(r"^Column|^Neig|^E|^S|^M")))
names(out)|>println
#out.BO_SUBTYP_TXT|>unique
writedf(out, "s200_tkle.txt")
#append vg parameters
ot = lwf.hydpar_wessolek_tab(out.BOART)|>rcopy
ot[!,:ksat] .= ot[!,:ksat]*fact
plot(ot.alpha, ot.ksat, seriestype=:scatter, label="vg")
#wa.qplot()
#x,y = map(x->replace(x,missing=>0), [ot.alpha, ot.ksat])
# x, y = map(x -> collect(skipmissing(x)), [ot.alpha, ot.ksat])
# r2 = round(cor(x, y)^2, digits=3)
# p = qqplot(x, y, qqline = :fit)
# annotate!(p,:bottomright, text("R² = "*string(r2), 
#     :black, :left, 16))
wa.qplot(ot.alpha, ot.ksat)
ot
ot.alpha|>plot
vgout = hcat(out,ot)

#replace Missvals of O to 
# classes after wösten 1999
# Topsoils
# Coarse
# Medium
# Mediumfine
# Fine
# Very Fine
lwf.hydpar_hypres_tab(texture = ["C","M","MF","F","VF"], topsoil = fill(true, 5))
lwf.hydpar_hypres_tab(texture = ["C","M","MF","F","VF"], topsoil = fill(false, 5))

mf = lwf.hydpar_hypres_tab(texture = ["MF"], topsoil = false )|>rcopy
fact = 1/(1000*24*60*60)
mf.ksat = mf.ksat*fact
mf
@rsubset vgout :BOART=="O"
#now replace 
for i in eachrow(vgout)
    if isequal(i.BOART, "O")
        i.ths   = only(mf.ths)
        i.thr   = only(mf.thr)
        i.ksat  = only(mf.ksat)
        i.alpha = only(mf.alpha)
        i.npar  = only(mf.npar)
        i.mpar  = only(mf.mpar)
        i.tort  = only(mf.tort)
    end
end

any(ismissing.(vgout.ksat))
findmin(vgout.ksat)
vgout[6,Cols(r"^B")]
@rsubset hor :BF_ID==10461
@rsubset vgout :BF_ID==10461

#now replace last soil in 10461 with CLAY
mf = lwf.hydpar_wessolek_tab("Tt")|>rcopy
fact = 1/(1000*24*60*60)
mf.ksat = mf.ksat*fact

for i in eachrow(vgout)
    if ismissing(i.BOART) #&& isequal(i.BF_ID, 10461)
        i.BOART = "Tt"
        i.ths   = only(mf.ths)
        i.thr   = only(mf.thr)
        i.ksat  = only(mf.ksat)
        i.alpha = only(mf.alpha)
        i.npar  = only(mf.npar)
        i.mpar  = only(mf.mpar)
        i.tort  = only(mf.tort)
    end
end
@rsubset vgout :BF_ID==10461
x =  @rsubset vgout :BOART=="O"
x[!,Cols(r"^t")] 

mf = lwf.hydpar_wessolek_tab("Hh")|>rcopy

any(ismissing.(vgout.ksat))
vgout[findmin(vgout.ksat)[2],Cols(r"^B")]
vgout[ismissing.(vgout.ksat),Cols(r"BF_ID|^HO|thr|^BODSYSTEINH_TXT|BOART")]
x = @rsubset vgout :BF_ID==5438
plot(x.ksat,x.ths,seriestype=:scatter, label="vg")

vgout.thickness .= vgout.UTIEF .- vgout.OTIEF
vgout.thickness .= vgout.thickness .* 0.1       #cm to m
writedf(vgout, "s200_tkle_vg.txt")

gr = groupby(vgout, :BF_ID)
#plot cumsum thickness for each group
p = plot()
for i in gr
    plot!(p,
    i.thickness,
    cumsum(i.thickness), 
    seriestype=:scatter,
    #label=string.(unique(i.BF_ID)))
    label=false)
end
display(p)

#which soil has the most layers?
findmax(vgout.HOR_NR)
vgout[findmax(vgout.HOR_NR)[2],:]
vgout[findmax(vgout.thickness)[2],:]

gr[2].thickness|>cumsum|>plot
a,b,anns = [],[],[]
for i in gr
    push!(a,i.thickness|>last)
    push!(b,cumsum(i.thickness)|>last)
    push!(anns, string.(unique(i.BF_ID)...))
end
# Create the annotations
anns_text = [text(ann, 7, :center, :red, rotation=-30) for ann in anns]
plot(a,b,
    seriestype=:scatter,
    markersize=.5,
    annotations=(a,b,anns_text),
    label=false)


@doc pop! #See also: popfirst!, popat!, delete!, deleteat!, splice!, and push!.
# Remove an item in collection and return it.

x = @rsubset vgout :BOART=="O"
plot(x.OTIEF, x.UTIEF, seriestype=:scatter, label="vg")

plot(boxplot(vgout.ksat|>skipmissing, label="ksat"),
    boxplot(vgout.ths|>skipmissing, label="ths"),
    boxplot(vgout.thr|>skipmissing, label="thr"),
link=:x,xticks=false)

#@df vgout boxplot(cols(:ths, :thr, :ksat),link=:x,xticks=false)


@doc leftjoin
#A left join includes all rows from df1.
supro = leftjoin(subset_hor,pro, on=:GEN_ID, makeunique=true)
df2 = rightjoin(subset_hor,pro, on=:BF_ID, makeunique=true)
#:BF_ID onetomany
dropmissing!(df2, :BOART)
sel = @rsubset df2 :BF_ID==3212
select!(sel,Not(Cols(r"^Column|^Neig")))
pretty_table(sel[!,17:end])




println(names(supro))
supro = leftjoin(supro,zuo,on=:GEN_ID, makeunique=true)
println(names(supro))
out = leftjoin(supro,leg,on=:TKLE_NR, makeunique=true)
select!(out,Not(Cols(r"^Column")))
println(names(out))
plot(out.BF_ID, out.BF_ID_1, seriestype=:scatter, label="bfids")

wes = CSV.read("soildata_wessolek_genid.txt", DataFrame)
out2 = leftjoin(out,wes,on=:BF_ID, makeunique=true)
plot(out2.ths, out2.ksat, seriestype=:scatter, label="vg")
out2 = filter(row -> !ismissing(row.ksat), out2)
println(names(out2))
writedf(out2, )
a = unique(out2.BOART)
b = unique(out2.bart)
plot(boxplot(out2.UTIEF,label="1st"), 
    boxplot(out2.UTIEF_1,label="2nd"), 
    boxplot(out2.UTIEF_2,label="3rd"))
select(out2,r"BF")
select(out2,r"UTIE")

#all_same = unique(out2[:, :UTIEF]) == unique(out2[:, :UTIEF_2]) == unique(out2[:, :UTIEF_1])
#println(all_same ? "All values are the same" : "Values are not all the same")

writedf(out2, "vgdata_big.txt")




using PyCall
gpd = pyimport("geopandas")
gdf = gpd.read_file("D:/Bodendaten/buek200_2020/merged_buek200.shp")
shp = gpd.read_file("D:/Wasim/regio/rcm200/v12/catchment.shp")
gdf = gpd.clip(gdf, shp)
gdf.TKLE_NR
cd("D:/Wasim/regio/rcm200/v12/")
gdf.to_file("v12_merged_buek200.shp")
#tk = convert(Array,(gdf.TKLE_NR)) #Julia breaks.
tk = convert(Array,gdf.TKLE_NR) 
using RCall
@rimport terra as tr
Z = tr.vect("D:/Bodendaten/buek200_2020/merged_buek200.shp")
da = tr.subset(Z,select="TKLE_NR")


pd = pyimport("pandas")
fn="D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/vgdata_clean.txt"
# vdf = pd.read_csv(fn,delim_whitespace=true,usecols=[0,1,2,4,5,6],
#     skip_blank_lines=true, encoding="cp1252")
#vdf=pd.read_csv(fn,delim_whitespace=True,
#     usecols=list(range(0,60,1)),
# skip_blank_lines=True,encoding="cp1252")
vdf=pd.read_table(fn,usecols=0:40,
    skip_blank_lines=true,encoding="cp1252")

#vdf=pd.read_table(fn,   usecols=list(range(0,60,1)),
#skip_blank_lines=True,encoding="cp1252")


#dat = gpd.sjoin(gdf, vdf, how="inner")
#ValueError("'right_df' should be GeoDataFrame, got <class 'pandas.core.frame.DataFrame'>"
pd.merge(gdf, vdf, how="inner", on="TKLE_NR")



lk="D:/Bodendaten/buek200_2020/merged_buek200.shp" 
using GeoDataFrames, Shapefile
# Load the shapefile
#shp = Shapefile.Table(lk)
# Convert the shapefile to a GeoDataFrame
gdf = GeoDataFrames.read(lk)
# Assuming `out` is the DataFrame you want to join with `gdf`

# Perform the spatial join
joined = sjoin(out, gdf, :inner, :intersects, on=:TKLE_NR)
using Conda
Conda.add("geopandas")

x = "D:/Wasim/regio/out/rc200/x22/f19/Wolfsmuenster-qoutjl"
df = waread(x)
kge(df)

for yr in unique(year.(df.date))
    println(year)
    dx = filter(row -> year(row.date) == yr, df)
    dx[!, :year] = year.(dx[!,:date]);
    #dm = DataFrames.combine(groupby(dx, :year), y .=> sum .=> y);
    kge2(dm)|>println
end

dm = yrsum(df)
#DataFrames.combine(groupby(dm, :year), y .=> kge .=> y)
cb(dm)


for yr in dm.year #eachrow(dm)
    println(yr)
    x = filter(row -> year(row.date) == yr, df)
    println(kge(yr[!,3],yr[!,2]))
end

#kge2(dm[!,3],dm[!,2])

# Add a year column to the DataFrame
df[!, :year] = year.(df[!, :date])
# Group the DataFrame by year and calculate the KGE for each group
grouped_df = groupby(df, :year)
vs = DataFrames.combine(grouped_df) do group
    simulated, observed = vec(Matrix(group[!,Cols(1)])),vec(Matrix(group[!,Cols(2)]))
    return (kge2(simulated, observed))
end

"""
daily df input, yearly kge nse ve as DF
"""
function byear(x::DataFrame)
    df = copy(x)
    df[!, :year] = year.(df[!,:date]);
    grouped_df = groupby(df, :year)
    DataFrames.combine(grouped_df) do group
        simulated, observed = vec(Matrix(group[!,Cols(1)])),vec(Matrix(group[!,Cols(2)]))
        kge = wa.kge2(simulated, observed)
        ve = wa.vef(simulated, observed)
        nse = wa.nse(simulated, observed)
        #grouping key is returned as first column
        dout = DataFrame(kge=kge,nse=nse,ve=ve)
        return dout
    end
end

x = "D:/Wasim/regio/out/rc200/x22/f19/Wolfsmuenster-qoutjl"
x = "D:/Wasim/Tanalys/DEM/brend_fab/out/m8/Schweinhof-qoutjl"
df = waread(x)  
cdof(df)
@time setup()

wa.byear(r"qoutjl")
wa.byear(r"Sch+.*qoutjl")
ed = r"^evar"|>dfr
wa.byear(dropmissing(ed))
wa.byear(r"^tem")

for z in dfonly(r"qoutjl$")
    wa.byear(z)|>println
end

s="D:/Wasim/regio/out/rc200/x22/f18/Wolfsmuenster-qoutjl"
wa.byear(s)


# af = groupby(df, :year)
# #transform(af, x -> (x=kge(x.C22,x.Wolfsmuenster)), keepkeys=false)
# DataFrames.transform(df, names(df)[1:end-2] .=> mean)
# #DataFrames.transform(df, names(df)[1:2] .=> kge)

s = (raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\vgdata_big.txt")
ot = CSV.read(s, DataFrame)
s = "D:/Wasim/regio/rcm200/v12/tkle_v12.txt"
tk = CSV.read(s, DataFrame;header=true)

joi = rightjoin(ot, tk, on=:TKLE_NR, makeunique=false)

un = unique(ot.BF_ID)
unique(ot.GEN_ID)

tr = unique(ot.thr)
ts = unique(ot.ths)

ag = groupby(joi, :GEN_ID)
tg = copy(ag[2])
@rsubset(tg, :HOR_NR .∈ Ref(unique(tg.HOR_NR)))
#not same
@rsubset(tg, :HOR_NR in Ref(unique(tg.HOR_NR)))
unique(tg, :HOR_NR)
out=[]
for i in ag
    push!(out, unique(i, :HOR_NR))
end
map(size,out)
cdof(s)
xo = reduce(vcat, out)
cls = ["GROBBOD_F", "BF_ID_1", "BOF_NR", "STATUS", "BODTYP", "BO_SUBTYP", "BO_SUBTYP_TXT", "BODSYSTEINH", "BODSYSTEINH_TXT", "SUBSTRSUBTYP", "NeigDomin", "NEIGMax", "NEIGMin", "EXPOS", "RLFORM", "KULTUR", "EROSI", "EGRAD", "HUFORM", 
"GWS", "SPEZGW", "VNGRAD", "OEKFEU", "MHGW", "MGW", "MNGW", "FLANT_SPANNE", "FLANT_MITTELW", "TKLE_NR", "TK", "LE_NR", "TK_1", "LE_NR_1"]
select!(xo,Not(cls))
xo = xo[1:nrow(xo)-1,:]

cd()
cd(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8")
println(names(xo))
writedf(xo, "vgdata_2.txt")
k = select(xo,Cols(r"^TKL"))
#rename!(k,Dict(:TKLE_NR_1=>:TKLE_NR))
rename!(k,Dict(:TKLE_NR_2=>:TKLE_NR))
innerjoin(k,tk,on=:TKLE_NR,makeunique=true)
innerjoin(ot,tk,on=:TKLE_NR,makeunique=true)

using DataFrames

function remove_duplicate_columns(df::DataFrame)
    cols_to_keep = Bool[]
    colnames = names(df)
    
    for i in 1:length(colnames)
        if !(any(isapprox.(df[!, i], df[!, j]) for j in 1:i-1))
            push!(cols_to_keep, true)
        else
            push!(cols_to_keep, false)
        end
    end
    
    return df[:, cols_to_keep]
end
sa = remove_duplicate_columns(xo)


ot[!,[:OTIEF_1, :OTIEF]]
ot[!,Cols(r"^OTIEF")]

nd = unique(ot,:BF_ID)
nd[!,Cols(r"^OTIEF")]

plot()
for z in out
    bar!([findmax(z.UTIEF)[1]],label=unique(z.GEN_ID))
end
plot!(legend=:bottomright)

z=out[6]
bar([findmin(z.UTIEF)[1],findmax(z.UTIEF)[1]],label=unique(z.GEN_ID))
bar([findmax(z.UTIEF)[1]],label=unique(z.GEN_ID))


ds = select(xo,1:9)
transform(ds, :UTIEF => ByRow(findmax) => :max)

#group and find max
ds = []
for group in groupby(xo, :GEN_ID)
    max_UTIEF = maximum(group.UTIEF)
    t = filter(row -> row.UTIEF == max_UTIEF, group)
    push!(ds, t)
end

#plot(xo.thickness,xticks=xo.GEN_ID,group=xo.GEN_ID)

ud = unique(xo,:thickness)
scatter(ud.thickness, 
    annotations=(1:length(ud.thickness),
    ud.thickness,unique(ud.GEN_ID)),
    legend=false)

# #annotations = (ds.name,ds.KGE, ann, :top),
# bar(ud.thickness,
#     annotations=(1:length(ud.thickness),
#         ud.thickness,
#         Plots.text(
#             string(ud.LE_KURZ_1)[18:30],        
#             7, :black, :center, 
#         halign=:center, rotation=-35.0)),
#     legend=false)
bar(ud.thickness,
    annotations=(1:length(ud.thickness),
        ud.thickness,
        string(ud.LE_KURZ_1)[18:30],7),
    legend=false)

#group and find max
ds = []
for group in groupby(xo, :GEN_ID)
    max_t = maximum(group.ths)
    t = filter(row -> row.UTIEF == max_t, group)
    push!(ds, t)
end

using StatsPlots
# Concatenate the data frames
all_data = vcat(ds...)
# Plot the grouped bar chart
@df all_data groupedbar(
    #xticks = :GEN_ID,
    1:nrow(all_data), 
        cols(:OTIEF), group = :UTIEF)

ctlg("D:/Wasim/regio/control/","v12")

setup()

# Load the numpy array
lk="C:/Users/chs72fw/Documents/EFRE_GIS/Bodendaten/PTF/Rosetta3_update_KK/Rosetta3_code/names.df"
#nm = CSV.read(lk, Vector;)
#nm = CSV.read(lk, DataFrame;)
nm = readdlm(lk, String;header=false)
# Read the CSV file
results = CSV.read("C:/Users/chs72fw/Documents/EFRE_GIS/Bodendaten/PTF/Rosetta3_update_KK/Rosetta3_code/output/ka5_results.csv", DataFrame;header=false)
#(nm|>transpose)
nm = Symbol.(String.(nm[2,:]))
# Rename the columns
rename!(results, nm)
        # OUTPUT
        # theta_r [cm3/cm3]
        # theta_s [cm3/cm3]
        # alpha  [1/cm]
        # n
        # Ks in [cm/day]
# Convert cm/day to m/s
fact = 1.15740741e-7
results.ks = results.ks .* fact
# Read the second CSV file
fb = CSV.read("C:/Users/chs72fw/Documents/EFRE_GIS/Bodendaten/PTF/Rosetta3_update_KK/Rosetta3_code/input/KA5_Feinboden.csv", DataFrame)
# Convert 1/cm to 1/m
results.alpha = results.alpha .* 100
# Add a new column
results.bart = fb.nam

results
writedf(results, "ka5_rosetta.txt")
se = @rsubset vgout :BOART=="Ss"
wa.qplot(se.ths, se.ksat, seriestype=:scatter, label="vg")
plot(se.ksat, results.ks, seriestype=:scatter, label="vg")
diff = se.ksat[1] - results.ks[1]
bar([se.ksat[1] ,results.ks[1] , diff])
#x = se.ksat[1:length(results.ks)]
#wa.qplot(x, results.ks)

se.ksat[1] / results.ks[1]   #factor 19.38 what? 

se = @rsubset vgout :BOART=="Tt"

se.ksat[1] / results.ks[end]   #factor 1.9
diff = se.ksat[1] - results.ks[1]
kd = DataFrame(,:auto)

bar(Dict(["wessolek" ,"rosetta", "diff"].=>
[se.ksat[1] ,results.ks[1] , diff]),
    legend=false,
    xaxis=true,
    labels=true)

se = unique(vgout,:BOART)
se.BOART in results.bart

boa_results = results.bart
findall(x->x in boa_results, se.BOART)
sx = se[findall(x->x in boa_results, se.BOART),:]
rename!(sx, Dict(:BOART=>:bart))
sx = innerjoin(sx, results, on=:bart, makeunique=true)
#startswith(names(sx),"k")
select(sx,(Cols(r"^k")))
@df sx groupedbar(sx.bart, cols(:ks,:ksat), seriestype=:scatter, label=["a","b"])
#plot(sx.ksat, results.ks, seriestype=:scatter, label="vg")
@df sx groupedbar(sx.bart, cols(:ks,:ksat),
    yaxis = :log,
    xrotation = -35,  # Rotate x-axis labels for better readability
    group = sx.bart, label=false)

@df sx[1:4,:] groupedboxplot(:bart, cols(:ks,:ksat),
#    yaxis = :log,
    xrotation = -35,  # Rotate x-axis labels for better readability
    group = :bart, label=false)

@df sx[1:4,:] groupedbar(:bart, cols(:ks,:ksat),
    #    yaxis = :log,
        xrotation = -35,  # Rotate x-axis labels for better readability
        group = :bart, label=false)


sx.diff = sx.ksat .- sx.ks
@df sx bar(sx.bart, cols(:ks,:ksat,:diff),
    yaxis = :log,
    xrotation = -35,  # Rotate x-axis labels for better readability
    group = :bart)

@show plotattr(:Axis)
@show plotattr(:Series)

anns_text = [text(ann, 7, :center, :red, rotation=-30) for ann in string.(collect(unique(sx.bart)))]


plot(sx.ks, label="rosetta")
plot!(sx.ksat,label="wessolek")
plot!(sx.diff,seriestype=:scatter,markershape=:diamond,
    markersize=5,color=:blue,label="diff",yaxis=:identity,xaxis=false)
title!("PTF methoden - Ksat [m/s]")
#xlabel!(1:nrow(sx),string.(collect(unique(sx.bart))))
for i in 1:nrow(sx)
    annotate!(i,sx.diff[i],text(string(sx.bart[i]), 9, 
        :top, :blue, rotation=-20, halign=:left))
end
plot!(legend=:topright)
savefig("ptf_ksat.png")