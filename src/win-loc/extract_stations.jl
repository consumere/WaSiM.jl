


rpt = "D:/remo/qm/corgrids/pre/pre_cor.nc"
r = Raster(rpt;key= "pre")
xd = ncdf(r)
dubplot(xd) #no dubs!!!!
dfp(xd)
#theme()
dfm(xd;fun=yrsum,ann=false,mode=:bar)
dfm(xd;fun=monsum,ann=false,mode=:bar) #mit obs vergleichen!
hydromon(xd)
bargroup(xd;fun=monmean)
obsrast = "d:/remo/cordex/eobs/v28/pre/pre_rcm_obs.nc"
or = Raster(obsrast;key= "pre")
xdo = ncdf(or)
xdo.date .= xdo.date .+ Dates.Hour(12)
map(x->x.pre .= Float64.(x.pre),[xd,xdo])
mx = mall(xd,xdo)
bargroup(mx;fun=monsum)
dfm(mx;fun=monsum,ann=false,mode=:bar)
pwd()
wa.cmplot()

#### extract stations from netcdf
@vv "extract(st"
od = wa.lpro("D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt")
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
#st = Raster("D:/remo/cordex/eobs/v28/biascorr/pre_cor.nc";key="pre") #doesnt work in Rasters, due to dimension declaration
#st = Raster("D:/remo/cordex/eobs/v28/pre/pre_cor.nc";key="pre")
st = r
# Extract values from raster
dou = collect(Rasters.extract(st, pnts;atol=.1)|>skipmissing)
dx = DataFrame(dou)
#GeoDataFrames.getgeometrycolumns(dx)
first(dx)
#make station df

cvar = Symbol("pre")
k = dx[!, cvar][1]
df = DataFrame([Rasters.lookup(k,1),
	Float64.(k.data)],[:date,cvar])

df.date .= parse.(DateTime, string.(df.date))
@df df plot(df.date,cols(cvar))

#make station df loop
df = []
for i in 1:nrow(dx)
	k = dx[!, cvar][1]
	tmp = DataFrame([Rasters.lookup(k,1),
	Float64.(k.data)],[:date,cvar])
	tmp.date .= parse.(DateTime, string.(tmp.date))
	push!(df, tmp)
end
df|>first

##
ds = innerjoin(unique.(df, :date)..., on = :date, makeunique=true)
nn = map(x->replace(x, 
    "\xfc" => "ue","\xdf" => "ss",
    #r"[(-]" => "-",
    r"[(]|[)]" => "",
    r"[--]" => "-"),
    # r",.*" => "",
    # r"-.*" => ""),
    string.(od.name))
nn[map(x->startswith(x,"W"),nn)]

string.(od.name)[map(x->startswith(x,"W"),string.(od.name))]

rename!(ds, names(ds)[2:ncol(ds)] .=> nn)

cd(raw"D:\remo\qm\corgrids\pre")
wa.writewa("pre_stations.wa",ds)
first(ds)
##append headers.
crd = DataFrame(dx[!,1])
crd.name .= map(y->(y[1:2]*y[end-1]) ,nn)
rename!(crd,1=>:lon,2=>:lat)
crdf = permutedims(crd[!,[:name,:lon,:lat]])
#writedf("cords.wa",crdf)

fn="D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
station_crds = readlines(fn)[1:5]
# Split each line into fields
station_crds = [split(line, '\t') for line in station_crds]
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
pre_stations_lines = open("pre_stations.wa", "r") do f
    readlines(f)
end
# Write station_crds and all but the first line of pre_stations_lines to "pre_stations2.wa"
open("pre_stations2.wa", "w") do f
    for line in [station_crds; pre_stations_lines[2:end]]
        println(f, line)
    end
end
#or append raw header to pre_stations.wa 25832
station_crds = readlines(fn)[1:5]
open("pre_stations.wa", "w") do f
	for line in [station_crds; pre_stations_lines[2:end]]
		println(f, line)
	end
end


tb = selt(ds,2)
rglob("taub")
tb2 = fread(".\\cor_old\\Taube_i_pre_cor.tsv")
tb2.date = parse.(DateTime, string.(tb2.date))
tb2.date = tb2.date .+ Dates.Hour(12)
tbm = mall(tb,tb2)
qplot(tbm)
pwd()
savefig("qplot_tbb.png")
wa.bargroup(tbm;fun=monsum)


###########tas ###########################

"D:/remo/qm/corgrids/tas"|>cd
r = Raster("tas_cor_raw.nc";key="tas")
xd = ncdf(r)
dubplot(xd) #no dubs!!!!
first(xd)
#crds for extraction
fn="D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
od = wa.lpro(fn)
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
st = r
# Extract values from raster
dou = collect(Rasters.extract(st, pnts;atol=.1)|>skipmissing)
dx = DataFrame(dou)
first(dx)
#make station df
cvar = Symbol("tas")
k = dx[!, cvar][1]
df = DataFrame([Rasters.lookup(k,1),
	Float64.(k.data)],[:date,cvar])

df.date .= parse.(DateTime, string.(df.date))
@df df plot(df.date,cols(cvar),label=string(cvar),legend=:topleft)
#make station df loop
df = []
for i in 1:nrow(dx)
	k = dx[!, cvar][1]
	tmp = DataFrame([Rasters.lookup(k,1),
	Float64.(k.data)],[:date,cvar])
	tmp.date .= parse.(DateTime, string.(tmp.date))
	push!(df, tmp)
end
df|>first
##
ds = innerjoin(unique.(df, :date)..., on = :date, makeunique=true)
nn = map(x->replace(x, 
    "\xfc" => "ue","\xdf" => "ss",
    #r"[(-]" => "-",
    r"[(]|[)]" => "",
    r"[--]" => "-"),
    # r",.*" => "",
    # r"-.*" => ""),
    string.(od.name))
#check
nn[map(x->startswith(x,"W"),nn)]
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
wa.writewa("$cvar-stations.wa",ds)
first(ds)
##append headers.
crd = DataFrame(dx[!,1])
crd.name .= map(y->(y[1:2]*y[end-1]) ,nn)
rename!(crd,1=>:lon,2=>:lat)
crdf = permutedims(crd[!,[:name,:lon,:lat]])
#writedf("cords.wa",crdf)

station_crds = readlines(fn)[1:5]
# Split each line into fields
station_crds = [split(line, '\t') for line in station_crds]
x = round.(Vector(crdf[2,:]),digits=2)
y = round.(Vector(crdf[3,:]),digits=2)
# Replace fields in lines 3 and 4 in station_crds with fields from cords_wa
station_crds[3][5:end] = string.(x)
station_crds[4][5:end] = string.(y)
# Join the fields back together
station_crds = [join(fields, '\t') for fields in station_crds]
# Write the modified contents back to "station_crds.txt"
pre_stations_lines = open("$cvar-stations.wa", "r") do f
    readlines(f)
end
# Write station_crds and all but the first line of pre_stations_lines to "pre_stations2.wa"
open("$cvar-stations2.wa", "w") do f
    for line in [station_crds; pre_stations_lines[2:end]]
        println(f, line)
    end
end
#or append raw header to pre_stations.wa 25832
station_crds = readlines(fn)[1:5]
open("$cvar-stations.wa", "w") do f
	for line in [station_crds; pre_stations_lines[2:end]]
		println(f, line)
	end
end

##no check with obs.
tb = selt(ds,2)
fn = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/temp_1970.txt"
dfobs = waread2(fn) #not fn
tbobs = selt(dfobs,1)
tbobs.date .= parse.(DateTime, string.(tbobs.date)) .+ Dates.Hour(12)
mx = mall(tb,tbobs)|>dropmissing
qplot(mx)
wa.hyeval(mx)
@doc wa.hyeval(mx;freq="Y")
wa.hyeval(mx;freq="Y",yscale=:identity,fun=mean,ylab="[degC/year]")
wa.hyeval(mx;freq="Q",yscale=:log2,fun=mean)


############12 2023######################
sf = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
od = lpro(sf)
wa.stplot(sf)
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
rasterpath= "d:/remo/qm/corgrids/pre/pre-cor.nc"
st = Raster(rasterpath;name=:pre)
plot(st[Ti=5])
dou = collect(extract(st, pnts;atol=0.1))
#dx = DataFrame(dou)|>dropmissing
dx = dropmissing(DataFrame(dou),:pre)
using GeoDataFrames
GeoDataFrames.getgeometrycolumns(dx)
#make station df loop
myvar = propertynames(dx)[2]
df = []
for i in 1:nrow(dx)
	k = tovec(dx,2)[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,myvar]))
end
df
ds = innerjoin(df..., on= :date, makeunique=true)
nn = map(p->replace(p, 	"\xfc" => "ue",    "\xdf" => "ss",
r",.*" => "",r"-.*" => "_"),string.(od.name))
#nn = map(y->(y[1:5]*y[end]),string.(nn))
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
ds.date .= parse.(DateTime, string.(ds.date))
obs = waread2(sf)
obs.date = parse.(DateTime, string.(obs.date)) .+ Dates.Hour(12)
#mx = mall(selt(ds,2),selt(obs,1)) 
mx = mall(selt(ds,Cols(r"Sand")),selt(obs,Cols(r"Sand"))) 
wa.hyeval(mx;freq="Y",yscale=:identity,fun=sum,ylab="[mm/year]")
wa.dfm(mx;fun=yrsum,ann=false,mode=:bar)
#@doc wa.hyeval(mx;freq="Q",yscale=:log2,fun=mean)
ys = yrsum(mx)|>dropmissing
@doc vef(ys)
dfbar(ys)
vef(tovec(ys,2),tovec(ys,3))
vef(tovec(mx,1),tovec(mx,3))

cdof(rasterpath)
pwd()
outname = string(myvar)*"_ext.wa"
wa.writewa(outname,ds)
# op()

xo = fread(raw"D:\remo\qm\corgrids\pre\pre_stations2.wa")
xo.date = parse.(DateTime, string.(xo.date)) .+ Dates.Hour(12)
mx = mall(xo,ds)
qplot(mx)
tline(mx)
############now eobs######################

sf = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/preci_1970.txt"
od = lpro(sf)
#wa.stplot(sf)
first(od.geometry,5)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))

#rasterpath = "D:/remo/cordex/eobs/v28/rcm_eobs/rr_wa.nc"
#plot(st[t=5])
rasterpath = "D:/remo/cordex/eobs/v28/rr_ens_mean_v28.0e_crop.nc"
st = Raster(rasterpath;name=:rr)
plot(st[Ti=5])
cdof(rasterpath)
dou = collect(extract(st, pnts;atol=0.1))
dx = dropmissing(DataFrame(dou),:rr)
#make station df loop
myvar = propertynames(dx)[2]
df = []
for i in 1:nrow(dx)
	k = tovec(dx,2)[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,myvar]))
end
ds = innerjoin(df..., on= :date, makeunique=true)
nn = map(p->replace(p, 	"\xfc" => "ue",    "\xdf" => "ss",
r",.*" => "",r"-.*" => "_"),string.(od.name))
#nn = map(y->(y[1:5]*y[end]),string.(nn))
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)
ds.date .= parse.(DateTime, string.(ds.date)) + Dates.Hour(12)


# ds = reorder_df(ds)
# obs = reorder_df(obs)

@doc wa.corrbar(ds,obs)
obs = waread2(sf)
obs.date = parse.(DateTime, string.(obs.date)) .+ Dates.Hour(12)
ds.date
#a,b = map(x->select(x,[5,10,12,15,17,20,ncol(x)]),[ds,obs])
a,b = map(x->selt(x,[5,10,12,15,17,20,ncol(x)]),[ds,obs])
# map(x->println(names(x)),(a,b))
# rename!(a, names(a) .=> names(b))
# map(x->println(names(x)),(a,b))
# b
# wa.corrbar(a,b)

mx = mall(selt(ds,Cols(r"Sand")),selt(obs,Cols(r"Sand"))) 
wa.hyeval(mx;freq="Y",yscale=:identity,fun=sum,ylab="[mm/year]")
wa.dfm(mx;fun=yrsum,ann=false,mode=:bar)

mx = mall(ds,obs)
wawrite(mx,"joined_eobs.wa")

pwc()

function xbar(a::DataFrame, b::DataFrame)
    # Load the DataFrames
    df_A = copy(a)
    df_B = copy(b)
    ti = try 
        nma = DataFrames.metadata(df_A)|>only|>last|>basename
        nmb = DataFrames.metadata(df_B)|>only|>last|>basename
        ti = "$nma vs $nmb"        
    catch
        @info    "no metadata in $a or $b !"
        ti = "$a vs $b"
    end

    # # Remove last column from df_B #has to have the same names
    # select!(df_B, Not(ncol(df_B)))
    # select!(df_A, propertynames(df_B)) 
    #subset by time
    # tcol_A = select(df_A, Cols(r"date|year|month|day"))
    # tcol_B = select(df_B, Cols(r"date|year|month|day"))
    df_A = yrsum(df_A)
    df_B = yrsum(df_B)

    # Define a function to replace missing with NaN in a vector
    replace_missing_with_NaN(v::AbstractVector) = replace(v, missing => NaN)

    #typeof(df_B[!,1])
    # Apply the function to each numeric column in df_B
    #df_B[!, :] .= [col isa Vector{Union{Missing, Float64}} ? replace_missing_with_NaN(col) : col for col in eachcol(df_B)]
    df_A = DataFrame([col isa Vector{Union{Missing, Float64}} ? replace_missing_with_NaN(col) : col for col in eachcol(df_A)], names(df_A))
    df_B = DataFrame([col isa Vector{Union{Missing, Float64}} ? replace_missing_with_NaN(col) : col for col in eachcol(df_B)], names(df_B))

    # la,lb = map(nrow,[tcol_A,tcol_B])
    # if la > lb
    #     df_A_subset = leftjoin(df_A, df_B, on = propertynames(tcol_A)|>first)
    #     df_B_subset = df_B
    # else
    #     df_A_subset = df_A
    #     df_B_subset = leftjoin(df_B, df_A, on = propertynames(tcol_B)|>first)
    # end


    # # Find the common date range
    # common_dates = intersect(df_A.date, df_B.date)
    # # Subset df_A and df_B to the common date range 
    # df_A_common = filter(x->x.date ∈ common_dates,df_A)
    # df_B_common = filter(x->x.date ∈ common_dates,df_B)



    common_dates = intersect(df_A.year, df_B.year)
    # Subset df_A and df_B to the common date range 
    df_A_common = filter(x->x.year ∈ common_dates,df_A)
    df_B_common = filter(x->x.year ∈ common_dates,df_B)

    #remove time column and copy back to df_A and df_B
    df_A = select(df_A_common, Not(Cols(r"date|year|month|day")))
    df_B = select(df_B_common, Not(Cols(r"date|year|month|day")))
    # if size not equal, return error
    if size(df_A) != size(df_B)
        @error("DataFrames must have the same size!")
        return
    end
    
    
    # Calculate correlations and replace NaN with 0
    #correlations = cor.(eachcol(df_A_subset), eachcol(df_B)).^2
    correlations = Vector{Float64}(undef, size(df_A, 2))
    for i in 1:size(df_A, 2)
        correlations[i] = cor(df_A[:, i], df_B[:, i])^2
    end
    replace!(correlations, NaN => 0)

    p0 = Plots.bar(1:size(df_A, 2), correlations,
        legend = false,
        title = ti,
        fillcolor = ifelse.(correlations .> 0.35, 
            "cornflowerblue", "coral2"),
        xticks = (1:size(df_A, 2), propertynames(df_A)),
        xrotation = 45,
        xlabel = "",
        ylabel = "Correlation",
        fontfamily="Computer Modern",
        left_margin = 10mm,
        bottom_margin = 2mm);

    ann = map(x->string.(round(x;sigdigits=2)),correlations)

    for i in 1:size(df_A, 2)
        Plots.annotate!(i,correlations[i],
        Plots.text(ann[i],9,:center,:top,:black))
        #println("R² "*ann[i]*" of Basin "*names(df_A)[i])
        #show("R² "*ann[i]*" of Basin "*names(df_A)[i]*" added")
    end

   return p0
end

names(mx)
mx[!,Cols(r"_$")]
xbar(selt(mx, r"_1"),selt(mx, r"_$"))

c,d = selt(b,r"Ki"),selt(a,r"Ki")
c
d

xbar(selt(b,r"K"),selt(a,r"K"))
xbar(selt(b,r"Ki"),selt(a,r"Ki"))

#@doc wa.hyeval(mx;freq="Q",yscale=:log2,fun=mean)
ys = yrsum(mx)|>dropmissing
dfbar(ys)
vef(tovec(ys,2),tovec(ys,3))
vef(tovec(mx,1),tovec(mx,3))








function correlation_heatmap(a::DataFrame, b::DataFrame)
    # Load the DataFrames
    df_A = copy(a)
    df_B = copy(b)
    ti = try 
        nma = DataFrames.metadata(df_A)|>only|>last|>basename
        nmb = DataFrames.metadata(df_B)|>only|>last|>basename
        ti = "$nma vs $nmb"        
    catch
        @info    "no metadata in $a or $b !"
        ti = "$a vs $b"
    end

    #subset by time
    # tcol_A = select(df_A, Cols(r"date|year|month|day"))
    # tcol_B = select(df_B, Cols(r"date|year|month|day"))
    df_A = yrsum(df_A)
    df_B = yrsum(df_B)

    common_dates = intersect(df_A.year, df_B.year)
    # Subset df_A and df_B to the common date range 
    df_A_common = filter(x->x.year ∈ common_dates,df_A)
    df_B_common = filter(x->x.year ∈ common_dates,df_B)

    #remove time column and copy back to df_A and df_B
    df_A = select(df_A_common, Not(Cols(r"date|year|month|day")))
    df_B = select(df_B_common, Not(Cols(r"date|year|month|day")))
    # if size not equal, return error
    if size(df_A) != size(df_B)
        @error("DataFrames must have the same size!")
        return
    end

    # Define a function to replace missing with NaN in a vector
    replace_missing_with_NaN(v::AbstractVector) = replace(v, missing => NaN)

    #typeof(df_B[!,1])
    # Apply the function to each numeric column in df_B
    #df_B[!, :] .= [col isa Vector{Union{Missing, Float64}} ? replace_missing_with_NaN(col) : col for col in eachcol(df_B)]
    dfa = DataFrame([col isa Vector{Union{Missing, Float64}} ? replace_missing_with_NaN(col) : col for col in eachcol(df_A)], names(df_A));
    dfb = DataFrame([col isa Vector{Union{Missing, Float64}} ? replace_missing_with_NaN(col) : col for col in eachcol(df_B)], names(df_B));

    # Calculate correlations
    correlations = cor(Matrix(dfa), Matrix(dfb))^2

    # Replace missing values with 0
    replace!(correlations, NaN => 0)

    # Create labels for the heatmap
    x_labels = names(dfa)
    y_labels = names(dfb)

    # Create the heatmap
    p = heatmap(x_labels, y_labels, correlations,
            aspect_ratio=1,
            color=:balance,
            clims=(-1, 1),
            xrotation=90,
            yticks=(1:length(y_labels), y_labels),
            xticks=(1:length(x_labels), x_labels),
            title="Correlation Heatmap",
            xlabel="DFA",
            ylabel="DFB",
            framestyle=:box)

    # Add annotations
    for i in 1:size(correlations, 1)
        for j in 1:size(correlations, 2)
            #annotate!(j, i, 
            annotate!(j + 0.5, i + 0.5,
                text(round(correlations[i, j], digits=2), 
                    8, :white))
        end
    end
    return(p)
end

#a,b = map(x->select(x,[1,2,5,10,12,15,17,20,ncol(x)]),[ds,obs])
a,b = map(x->selt(x,r"wass|sand|kiss|fladung"i),[ds,obs])
correlation_heatmap(a,b)

dtf = innerjoin(a,b,on=:date,makeunique=true)
wawrite(dtf,"obs_eobs.wa")

df = dtf
md = cor(Matrix(df[:, Not(:date)])).^2
replace!(md, 1.0 => missing)

wa.heat(dtf)
ls()
kd = fread(r"qq_ext")
wa.heat(kd)
wa.heat(r"qq_ext")
fn="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
fn="D:/Wasim/Tanalys/DEM/Input_V2/meteo/saale_spec_ordered.09"
wa.heat(fn;ann=false)
wa.heat(fn;)
da = fread(fn)
hydro(da)
hydromon(da)
dfp(da)
dfm(da|>dropmissing;fun=monsum)
mbx(da|>dropmissing)
tbx(da|>dropmissing)
yr = yrsum(da)

wa.mbx(da|>dropmissing;fun=sum)



####2024########
sf = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/temp.2021"
od = lpro(sf)
wa.stplot(sf)
first(od.geometry,5)
dropmissing(od, :geometry)
# Convert observations to points note: first, lat, then lon
pnts = collect((ArchGDAL.gety(o,0),ArchGDAL.getx(o,0)) for o in (od.geometry) if !ismissing(o))
#nams = collect(o for o in (od) if !ismissing(o))
rasterpath="D:/remo/qm/corgrids/tas/tas_cor_raw.nc"
st = Raster(rasterpath;name=:tas)
Plots.plot(st[Ti=5])
#Plots.plot(st[Y(Near(45.6)),X(Near(12))])
dou = collect(extract(st, pnts;atol=0.1))
dx = dropmissing(DataFrame(dou),:tas)
using GeoDataFrames
GeoDataFrames.getgeometrycolumns(dx)
myvar = propertynames(dx)[2]
dx.nam = od.name
df = []
for i in 1:nrow(dx)
	k = tovec(dx,2)[i]
	push!(df, DataFrame([Rasters.lookup(k,1),
		Float64.(k.data)],[:date,myvar]))
end
df
ds = innerjoin(df..., on= :date, makeunique=true)
nn = map(p->replace(p, 	"\xfc" => "ue",    "\xdf" => "ss",
r",.*" => "",r"-.*" => "_"),string.(od.name))
#nn = map(y->(y[1:5]*y[end]),string.(nn))
ds.date .= parse.(DateTime, string.(ds.date))
size(nn)
size(ds)
rename!(ds, names(ds)[2:ncol(ds)] .=> nn)

ds

obs = waread2(sf)
obs.date = parse.(DateTime, string.(obs.date)) .+ Dates.Hour(12)

mx = mall(selt(ds,Cols(r"Sommerach")),selt(obs,Cols(r"Sommerach"))) 
println(names(obs))
mx = mall(selt(ds,Cols(r"Ettleben")),selt(obs,Cols(r"Ettleben"))) 
dropmissing!(mx)
wa.hyeval(mx;freq="Y",yscale=:identity,fun=mean,ylab="[°C]")
wa.hyeval(mx;freq="Q",yscale=:identity,fun=mean,ylab="[°C]")

cd()
m2=dfr("tas-stations2.wa")
m2.date = parse.(DateTime, string.(m2.date)) .+ Dates.Hour(12)
z = mall(selt(ds,Cols(r"Bischbrunn")),
    selt(m2,Cols(r"Bischbrunn"))) 
hyeval(z;freq="Y",yscale=:identity,
    fun=mean,ylab="[°C]")

hyeval(z;freq="M",fun=mean,ylab="[°C]")

@pyjl

xr  = pyimport("xarray")
sh="d:/remo/qm/rsds/simh.nc"
simh=xr.open_dataset(sh)
#sh="d:/remo/qm/rsds/simh.nc"
sh="d:/remo/qm/rsds/rsds-cor.nc"
simp=xr.open_dataset(sh)
tl="D:/remo/cordex/eobs/v28/rsds/rsds_eobs.nc"
obsh=xr.open_dataset(tl)
pyjl.doyplot(simh,simp,obsh)

##rain
sh="d:/remo/qm/prec/simh.nc"
simh=xr.open_dataset(sh)
#sh="d:/remo/qm/rsds/simh.nc"
sh="d:/remo/qm/prec/pre_cm.nc"
sh="D:/remo/cordex/eobs/v28/pre/pre_cor.nc"
simp=xr.open_dataset(sh)
tl="D:/remo/cordex/eobs/v28/pre/pre_rcm_obs.nc"
obsh=xr.open_dataset(tl)
pyjl.doyplot(simh,simp,obsh;tosum=true)

##-> im ergebnisteil für jeden param so ein plot.




wa.dfm(mx;fun=yrsum,ann=false,mode=:bar)
ys = yrsum(mx)|>dropmissing
#@doc vef(ys)
dfbar(ys)
vef(tovec(ys,2),tovec(ys,3))
vef(tovec(mx,1),tovec(mx,3))

cdof(rasterpath)
pwd()
outname = string(myvar)*"_ext.wa"
wa.writewa(outname,ds)
