#lanuse reader
fn = raw"D:\Wasim\regio\control\rcm_adj_v5_loc.ctl"
lutable = open(fn) do io
    a = readbetween(io, "[landuse_table]", "\$JDVegReset")
    return(join(a[3:end-1]," ")) #rem first 2 and last lines
end


"""
controlfile::String = "",
match1::String = "[landuse_table]",
match2::String = "\$JDVegReset",
lutable::String = ""
see also landuse.jl for plotfuncs
"""
function lutab(;
    controlfile::String = "",
    match1::String = "[landuse_table]",
    match2::String = "\$JDVegReset",
    lutable::String = "")
    if controlfile != ""
        lutable = open(controlfile) do io
            a = readbetween(io, match1, match2)
            return(join(a[3:end-1]," ")) #rem first 2 and last lines
        end
    end

    entries = split(lutable, "}")
    entries = strip.(entries)
    entries = filter(s -> !isempty(s) && length(s) > 2 && !occursin(r"^(?i)#",s) && !occursin(r"^[[]",s), entries)
    map!(x->replace(x,r"\t" => " "),entries,entries)
    
    dictionary = Dict{Int, DataFrame}()
    
    for entry in entries
        # Extract key
        m = match(r"(\d+)\s+(\w+)\s+\{", entry)
        if m === nothing
            println("No match found in string: $entry")
            continue
        end
        key = parse(Int, m.captures[1]|>strip)
        nm = m.captures[2]|>strip
        println("check: $key $nm")
        
        # Extract value
        value_match = match(r"\{(.*)", entry)
        if value_match === nothing
            continue
        end
        value = value_match.captures[1]|>strip
        #value = replace(value, "method = MultipleHorizons;  EvapMaxDepth = 0.15;" => "")
        value = replace(value, r";" => ",")
        value = replace(value, r"\$" => "")
        value = replace(value, r"#" => "")
        value = replace(value, r"\s+" => " ")
        value = replace(value, r"^\s+|\s+$" => "")
        value = replace(value, r",$" => "")
        pairs = split(value, ",")
        pairs = lstrip.(pairs)
        pairs = filter(s -> !isempty(s), pairs)
        
        data = [split(strip(pair), " = ", limit = 2) for pair in pairs]
        #data = vcat([["class",nm]], data) #prepend to first pos
        data = vcat([["key",key]], data) #prepend to first pos
        #push!(data, ["class", nm]) #append last pos
        df = DataFrame(Param=first.(data), Value=last.(data))
        rename!(df, :Value => Symbol(nm))
        # Store key-value pair in the dictionary
        dictionary[key] = df
    end
    
    #return dictionary
    return [dictionary[i] for i in keys(dictionary)]
end

md = lutab(lutable)
md[7]
vgjl("lutab")

dfs = lutab(controlfile=fn)

#dfs = vcat([df[i] for i in keys(df)]...)
#dfs = [df[i] for i in keys(df)]
@rsubset dfs[5] :Param=="LAI"


function luplot(dfs::Vector{DataFrame};par::String="LAI")
    # Initialize an empty plot
    p1 = Plots.plot();
    # Iterate over all dataframes in dfs
    for (i, df) in enumerate(dfs)
        # Filter rows where :Param is "LAI"
        #lai_rows = @rsubset df :Param==par
        lai_rows = filter(row -> occursin(Regex(par,"i"), row[:Param]), df)
        key = @rsubset df :Param=="key"
        key = first(key[!,2])
        
        # Check if any rows were found
        if nrow(lai_rows) > 0
            # Parse "Wetland" column values into an array of integers
            lvs = parse.(Float64, split(first(lai_rows[!,2])))
            lab = replace(last(names(lai_rows)), r"_" => " ")

            # Add the values to the plot
            #Plots.plot!(p, lvs, label = "$i $lab key: $key")
            Plots.plot!(p1, lvs, label = "$lab [$key]",
                #seriestype=:bar,
                grid = false,
                #xlims=(1,12),
                xticks = (1:12, 
                [ monthabbr(x) for x in 1:12 ]),
                xrotation = 35
                ) 
        end
    end
    # Display the plot
    return p1
end

luplot(dfs;par="VCF")

dfs[1].Param|>println
luplot(dfs;par="AltDep")
luplot(dfs;par="rsc")
fn=raw"D:\Wasim\regio\control\rcm200_x22-genre-loc.ctl"
isfile(fn)|>println
dfs = lutab(controlfile=fn)

df = dfs[1]
@rsubset df :Param==r"VCF"i
filter(row -> occursin(r"VCF"i, row[:Param]), df)
luplot(dfs;par="rsc")
luplot(dfs;par="vcf")

# dat=vcat(
#     map(df->filter(row -> occursin(r"VCF"i, row[:Param]), 
#         df),dfs)...)

dat=map(df->filter(
    row -> occursin(r"VCF"i, row[:Param]), 
                df),dfs)
               
ad = DataFrame(
    class=last.(names.(dat)),
    vcf=map(z->parse.(Float64,
    split(z[!,2]|>first)),dat)
    )

stack(ad,:vcf)

ad2 = mapreduce(df -> DataFrame(
        class = last(names(df)),
        vcf = parse.(Float64, split(first(df[!, 2])))
    ), vcat, dat)

vcf_var = r"VCF"i

function luvars(dfs::Vector{DataFrame};par::String="LAI")
    reg_var = Regex(par,"i")
    ad = mapreduce(df -> begin
        filtered_df = filter(row -> occursin(reg_var, 
            row[:Param]), df)
        if nrow(filtered_df) > 0
            return DataFrame(
                class = last(names(filtered_df)),
                par = parse.(Float64, 
                split(first(filtered_df[!, 2])))
            )
        else
            return DataFrame(class = String[], 
            par = Float64[])
        end
    end, vcat, dfs)
    rename!(ad, :par => Symbol(par))
    return ad
end


ad = luvars(dfs;par="lai")
ad = luvars(dfs;par="vcf")

Plots.plot(ad[!, :lai], label = "lai", seriestype=:box)


ad = luvars(dfs;par="Z0")
cls = propertynames(ad)|>last
@df ad groupedbar(1:nrow(ad),cols(cls), 
        legend = :outertopright,
        xticks = false,
        grid = false,
        group = ad.class)
        #xrotation = 45) #,
        #xlabel = unique(ad.class))
        #ylabel = "[mm]", title = ti)

using Plots
# Create a box plot
Plots.boxplot(ad[!, :class], ad[!, cls], 
xrotation = 45, labels = false) #,orientation = :horizontal)
# Add annotations for each class
unique_classes = unique(ad[!, :class])
for (i, class) in enumerate(unique_classes)
    # Calculate the median Z0 value for this class
    median_z0 = median(filter(row -> row[:class] == class, ad)[!, :Z0])
    # Add an annotation for this class
    Plots.annotate!(i, median_z0, Plots.text(class, :center, 10))
end

# Display the plot
current()

#this is ET
fn="C:/Users/chs72fw/Documents/EFRE_GIS/Landcover/20181016_Deutschland_LC_clip_for_Ullmann/Urban/Urban-MOD16A2-006-results.csv"
cdof(fn)
glob("csv")
fn = "Urban-MCD15A2H-006-results.csv"
df = CSV.read(fn,DataFrame)
println(names(df))
grep(r"lai"i,names(df))
grep(r"FparLai_QC"i,names(df))
df.MCD15A2H_006_FparLai_QC|>first
#MCD15A3H.006	FparLai_QC	0	Good quality (main algorithm with or without saturation)
#filter qualitiy
ds = filter(row -> row[:MCD15A2H_006_FparLai_QC] .== 0, df)
#ds = filter(row -> row[:MCD15A2H_006_Lai_500m] .< 50, ds)
@df ds Plots.scatter(:Date, cols(:MCD15A2H_006_Lai_500m))


cdu()
fdi()
rglob(r"wood"i)
rglob(regand("wood","csv"))
#fn = rglob(r"wood"i)|>first
#cdof(fn)
#glob("csv")
fn = rglob(regand("wood","csv"))[3]
df = CSV.read(fn,DataFrame)
ds = filter(row -> row[:MCD15A2H_006_FparLai_QC] .== 0, df)
#ds = filter(row -> row[:MCD15A2H_006_Lai_500m] .< 50, ds)
@df ds Plots.scatter(:Date, cols(:MCD15A2H_006_Lai_500m))

ind=findmax(ds.MCD15A2H_006_Lai_500m)
ds[ind[2],:]
#@doc ArchGDAL.createpoint(ds[ind[2],:][4],ds[ind[2],:][3])
ArchGDAL.createpoint(ds[ind[2],:][4],ds[ind[2],:][3])
   
# a = gdf.DataFrame(
#     geometry=ArchGDAL.createpoint.(ds.Longitude, ds.Latitude), 
#     lai=ds.MCD15A2H_006_Lai_500m)
# a = unique(a, :geometry)
#besser,weil weniger punkte..
db = filter(x->year.(x.Date) .== 2018, ds)
a = DataFrame(
    geometry=gdf.createpoint.(db.Longitude, db.Latitude), 
    lai=db.MCD15A2H_006_Lai_500m)
plot(a.geometry, markersize=5, color=:red, legend=false)

pt = "D:/Wasim/regio/rcm200/v13/v13_4326.json"
gp = gdf.read(pt)
ex = reduce(gdf.union,gp.geometry)
plot!(ex, color=:blue, legend=false)
#subset all in ex
@doc gdf.within
@rsubset a 
isa(ex,ArchGDAL.AbstractGeometry)
isa(a.geometry,ArchGDAL.AbstractGeometry)
isa(a.geometry|>first,ArchGDAL.AbstractGeometry)
#Returns true if g1 is contained within g2.
#gdf.within(a.geometry, ex)
v = [ gdf.within(i, ex) for i in a.geometry ]
ptin = a[v,:]
plot!(ptin.geometry, 
    markersize=5, color=:green, legend=false)

##all in one.
fn = raw"C:\Users\chs72fw\Documents\EFRE_GIS\Landcover\20181016_Deutschland_LC_clip_for_Ullmann\Coniferous.Woodland\Coniferous-Woodland-MCD15A2H-006-results.csv"
#subset quality and select year 2018
ds = filter((x -> x[:MCD15A2H_006_FparLai_QC] .== 0 && year.(x.Date) .== 2018), CSV.read(fn,DataFrame))
a = DataFrame(
    geometry=gdf.createpoint.(ds.Longitude, ds.Latitude), 
    lai=ds.MCD15A2H_006_Lai_500m)

dsub = a[[ gdf.within(i, ex) for i in a.geometry ],:]
plot(ex,fillcolor=:transparent,title="Coniferous Woodland 2018")
plot!(dsub.geometry, 
    markersize=5, color=:green, legend=false)

@df dsub plot(:lai,xaxis=false, seriestype=:box, alpha=.66, title="Coniferous Woodland 2018",legend=false)


fn=raw"C:\Users\chs72fw\Documents\EFRE_GIS\Landcover\gee\lai-ee-chart.csv"
df = eeread(fn)
dfp(select(df,[1,4]))
dfp(select(df,[1,2]))
#df2 = CSV.read(fn,DataFrame)

fn = raw"D:\Wasim\regio\control\rcm200_x22-f18-cl.ctl"
fn = raw"D:\Wasim\regio\control\rcm200_x22-f18-spin-loc4.ctl"
da = lutab(controlfile=fn)
wa.luplot(da;par="LAI")

ad = luvars(da;par="lai")
ad = groupby(ad,:class)
#add month for each class
ad = DataFrames.combine(ad) do sdf
    transform(sdf, :class => (x -> 1:12) => :month)
end
@rsubset ad :month==6



m1 = r"^[[landuse_table]]"
m2 = r"[^\s+|\s+$]"
m2 = r"^JulDays"
m1 = r"^[[]landuse_table[]]"
m2 = r"^[[]"
a = readbetween(open(fn), m1, m2) #|>first

za = wa.lutab(controlfile=fn)


"D:/Wasim/regio/out/rc200/x22/spin/"|>cd
qgk()
kernelplot("route.txt")
klog("route.txt")


klog(a;lin=true)

df = readroute("route.txt") 
cols = Cols(3,5)
dsu = qba()
M = parse.(Float64,dsu[!,cols])|>y->subset(y,1=> ByRow(>(0.)))|>Matrix
mns = names(dsu[!,cols])
B = kde(M)
Plots.plot(B,legend=false,xlabel=mns[1],ylabel=mns[2])
rename!(dsu,:basin=>:sim)
dsu = innerjoin(dsu,df,on=:sim)
rename!(dsu,:name=>:basin)

wa.klog("route.txt";ann=false)
wa.kernelplot("route.txt")
findlog()

#dsu.basin  .= df.name
#dsu = innerjoin(dsu,df,on=:basin=>:sim)
    #renamecols=(:basin=>:name))
#name_mapping = df    
msk = @rsubset(parse.(Float64,dsu[!,cols]) .> 0)
nmn = dsu[msk[!,1],:]
nmn.basin=map(x->replace(x,r"_" => " "),nmn.basin)
Plots.annotate!([(M[i,1], M[i,2], 
    Plots.text(
    nmn.basin[i], 
    8, :left, 
    halign=:center, 
    rotation=-15.0)) for i in 1:size(nmn, 1)])
Plots.plot!(legend=false)


fx=raw"D:\Wasim\regio\control\srcm250.c2.ctl"
rd2=lutab(controlfile=fx)
luplot(rd2)
luplot(rd2;par="rs_evaporation")
## SOIL surface resistance in s/m (for evaporation only)
evp = luvars(rd2;par="rs_evaporation")
stack(evp)
ad = evp
cls = propertynames(ad)|>last
@df ad groupedbar(1:nrow(ad),cols(cls), 
        legend = :outertopright,
        xticks = false,
        grid = false,
        group = ad.class)

#again 
fx = raw"D:\Wasim\Testdata\Control\sample_control_file.ctl"
rd2 = lutab(controlfile=fx)
luplot(rd2)
luplot(rd2;par="rs_evaporation")
ad = luvars(rd2;par="rs_evaporation")
cls = propertynames(ad)|>last
@df ad groupedbar(1:nrow(ad),cols(cls), 
        legend = :outertopright,
        xticks = false,
        grid = false,
        group = ad.class)

@view ad[findmax(ad.rs_evaporation)[2],:]
extrema
z = luvars(rd2;par="Z0")
transform!(z, :Z0 => (x -> extrema(x)) => :extrema)
#@rusbset z :Z0 .> 0
a = select(rd,:ksat)
transform!(a, :1 => (x -> extrema(x)) => :extrema)
#first(a.extrema)|>bar

df = dfr("D:/Wasim/Tanalys/DEM/brend_fab/out/allwurz")

setup()
fn = raw"D:\Wasim\Tanalys\DEM\brend_fab\control\fab_c1.ctl"
df = lutab(controlfile=fn)
luplot(df;par="VCF")
ros = [findindf(i,"rc") for i in df]
vcat(ros...)
cdof(fn)
grep("set \$Def",readlines(fn))
cd(raw"D:\Wasim\Tanalys\DEM\brend_fab\out\c1\spin7")
dfp("wurz")
cdu()
dfs = rglob(r"wurzfab")
dfs = map(waread2,dfs[Not(endswith.(dfs,"nc"))])
#mall(dfs)

luplot(df;par="RootDist")

luvars(df;par="RootDist")