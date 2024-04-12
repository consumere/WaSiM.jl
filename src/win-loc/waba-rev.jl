#reworked waba 

af = filter(x -> occursin(r"^so_", x), readdir(pwd()))
if length(af) <= 2 || any(ismissing.(af))
    error("match failed \n ... abort.\nno special output files present!")
    exit(86)
end

println("calculating yearly water balance of WaSiM special output data..\n
bwvr = rain + snow + uprs - perc - qb - qd - qi - etr_ - ei_ - etrs_\n")


re = Regex("preci|snow_storage_tota|Capi|Perc|baseflo|directflo|interflo|real_evap|real_tran|interception_evaporatio|snow_evaporatio","i")
my = filter(x -> occursin(re, x), af)
printstyled("loading...\n $my\n",color=:green)

if (length(my) .!= 11)==true
    lng=length(my)
    printstyled("found only $lng files...\n $my\n",color=:yellow)
    error("\nfiles are missing!\ncheck so files...")
    exit(86)
end

using DataFrames, CSV, Dates 
#Distributions, Printf, Statistics

function so_read(x::AbstractString)
    "--- reader with drop exept of first col ---"
    ms=["-9999","lin","log"]
    df::DataFrame = CSV.read(x,DataFrame,
    missingstring=ms,
    types = Float64,
    delim="\t",
    silencewarnings=true,
    normalizenames=true,
    drop=(i, nm) -> i == 4) |> dropmissing
    df.YY=map(x ->Int(x),df.YY);
    df.MM=map(x ->Int(x),df.MM);
    df.DD=map(x ->Int(x),df.DD);
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    #df=df[:,Not(1:3)]
    df=df[:,Cols(4,end)]
    DataFrames.metadata!(df, "filename", x, style=:note);
end
#df=so_read(x)
#df=readdf(x)

rain = filter(x -> occursin("precip", x), af)|>only|>so_read
snow = filter(x -> occursin("snow_storage_total", x), af)|>only|>so_read
uprs = filter(x -> occursin("Capil", x), af)|>only|>so_read
perc = filter(x -> occursin("Perco", x), af)|>only|>so_read
qb = filter(x -> occursin("baseflow", x), af)|>only|>so_read
qd = filter(x -> occursin("directflow", x), af)|>only|>so_read
qifl = filter(x -> occursin("interflow", x), af)|>only|>so_read
etr = filter(x -> occursin("real_evapo", x), af)|>only|>so_read
etrans = filter(x -> occursin("real_trans", x), af)|>only|>so_read
ei = filter(x -> occursin("interception_evaporation", x), af)|>only|>so_read
etrs = filter(x -> occursin("snow_evaporation",x), af)|>only|>so_read

l = [rain, snow, uprs, perc, qb, qd, qifl, etr, etrans, ei, etrs]
#typeof(l)
nm = [names(l[i])[1] for i in 1:size(l, 1)]
#same:
#map(x->names(x)[1],l)

#println(string.(nm))
println("loaded dataframes:\n$nm")
#collect(nm)

#if any(length.(l) .!= 1)
    #missing_files = 
    #error("$(missing_files) is missing!\ncheck so files...")

if (length(l) .!= 11)==true
    error("files are missing!\ncheck so files...")
end

function mall(files::Vector{DataFrame})
    "reduces + merges by date"
    df = reduce((left, right) -> 
      innerjoin(left, right, on = :date,makeunique=true), 
      files)
    return(df)
end

function writedf(file, table)
    CSV.write(file, table, transform = (col, val) -> something(val, missing),delim="\t")  
    nothing
end

function writewa(file::AbstractString, df::DataFrame)
    dout = df
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 0
    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    CSV.write(file, dout, transform = (col, val) -> something(val, missing),delim="\t")  
    nothing
end

d = mall(l)
xd=copy(d)
writewa("waba-input.wa",xd)

function yrsum(x::DataFrame)
    df = x
    y = filter(x->!occursin(r"date|year",x),names(df))
    s = map(y -> Symbol(y),y)
    df[!, :year] = year.(df[!,:date]);
    df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
    return(df_yearsum)
end

dyr = yrsum(d)

pos = d[!,Cols(r"date|^(prec)|^(snow_stora)|^(Capi)")]
pos = yrsum(pos)
# calculate the sum of each row
psum = DataFrame(
    possums = [sum(eachrow(pos)[i]) for i in 1:size(pos, 1)],
    year=pos[!,:year]
)

neg = d[!,Not(Cols(r"^(prec)|^(snow_stora)|^(Capi)"))]
neg = sum.(yrsum(neg))

nsum = DataFrame(
    negsums = [sum(eachrow(neg)[i]) for i in 1:size(neg, 1)],
    year=neg[!,:year]
)


bw = innerjoin(psum, nsum, on=:year)
bw.bw = bw[!,:possums] .- bw[!,:negsums]

##########plotting ########################################
#using StatsPlots, Plots.PlotMeasures

using Plots
using StatsPlots: @df
using Plots.PlotMeasures: mm

fact=.95
plotsize = (1600*fact,800*fact)

ti="water-balance of "*basename(pwd())

#theme(:vibrant)
Plots.theme(:wong)
#theme(:wong2)
#theme(:wong2)
# showtheme(:dao)
ann = map(x->string.(round(x;sigdigits=3))*" [mm]",bw.bw)

p1 = @df bw Plots.plot(
    :year,:bw,
    legend = false, 
    seriestype=:bar,
    xticks = bw.year,
    xtickfont = 12,
    xlabel = "",
    ylabel = "[mm]",
    ytickfont = 12,
    title = ti,
    fillcolor = ifelse.(bw.bw .> 0, "cornflowerblue", "coral2"),
    size=plotsize,
    left_margin = 10mm,
    bottom_margin = 10mm, 
    xrotation = 45)
    
for i in 1:length(bw.year)
    Plots.annotate!(bw.year[i],bw.bw[i],(ann[i],10,:center,:top,:black))
    println(ann[i]*" added")
end

outname="waba-jl.png"
writedf("waba-input.yrly",dyr)
writedf("waba-year.yrly",bw)

Plots.savefig(p1,outname)
printstyled("$outname saved!\n",color=:green)