#water balance julia.
#"/mnt/d/Wasim/regio/out/rc200/v4"|>cd
#"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/c6/new"|>cd
#"/mnt/d/temp/saale/output/v6"|>cd
#"/mnt/d/Wasim/Goldbach/revision/fab150/pestpp/v5/pestout/so-files/wabasel"|>cd
#"/mnt/d/Wasim/Tanalys/DEM/brendpest/out_v7/spin"|>cd

af = filter(x -> occursin(r"^so_", x), readdir(pwd()))
if length(af) <= 2 || any(ismissing.(af))
    error("match failed \n ... abort.\nno special output files present!")
    exit(86)
end

println("calculating yearly water balance of WaSiM special output data..\n
bwvr = rain + snow + uprs - perc - qb - qd - qi - etr_ - ei_ - etrs_\n")


#perl -lne 'print substr($&, 1, -2)while/(?:(["'\''])(??{q<(?:[^\\\\>.$1.q<]|\\\\.)*>.$1}))/gx;' kk|tr '\n' '|' |cb


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
#default(show = true)
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

#names(d)
#xd = hcat(d[!,:date],d[!,Not(:date)])
#writedf("waba-input.wa",d)

# d2=d

# pos = d2[!,Cols(r"^(prec)|^(snow_stora)|^(Capi)")]
# neg = d2[!,Not(Cols(r"^(prec)|^(snow_stora)|^(Capi)|date"))]

# #neg[!,Not(1)] .* -1
# neg = neg .* -1

# xd = hcat(d2[!,:date],pos)
# xd = hcat(xd,neg)
# names(xd)
# rename!(xd,1=>"date")
# bardf(xd)

# dlist = Dict{Any,Any}()
# for i in l
#     dlist[i] = so_read(i)
# end
# dlist|>show

function yrsum(x::DataFrame)
    df = x
    y = filter(x->!occursin(r"date|year",x),names(df))
    s = map(y -> Symbol(y),y)
    df[!, :year] = year.(df[!,:date]);
    df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
    return(df_yearsum)
end

dyr = yrsum(d)

# #posbal = findfirst(x -> occursin(r"^(Perc)", x), names(dyr))
# #propertynames(dyr)
# pos = dyr[!,Cols(r"^(prec)|^(snow_stora)|^(Capi)|year")]
# neg = dyr[!,Not(Cols(r"^(prec)|^(snow_stora)|^(Capi)"))]

# xn = permutedims(neg[!,Not(1)])
# rename!(xn,Symbol.(neg[:,1]))
# #collect(eachcol(xn))|>sum
# col_sums = sum.(eachcol(xn))
# #col_sums = sum.(eachcol(neg))

# #filterplot("int",l)
# id = filterdf("int",l)
# sum.(id[!,1])
# select(id, names(id)[1] .=> ByRow(sum))


# l[11]|>dfp
# l[9]|>dfp
# l[10]|>dfp
# l[5]|>dfp
# l[1]|>dfp
# for x in 2:length(l)
#     println("adding",names(l[x]))
#     display(dfp!(l[x]))
# end

# df = pos
# names(df)|>print
# y = filter(x->!occursin(r"date|year",x),names(df))
# df_csum = combine(groupby(df, :year), y .=> sum .=> y);
# # Compute cumulative sum for each row of columns A and B
# df = df_csum[!,Not(:year)]
# # calculate the sum of each row
# df_sum = DataFrame(sum_each_row = [sum(eachrow(df)[i]) for i in 1:size(df, 1)])
# # add the new column to the original dataframe
# hcat(df, df_sum)

# df = pos
# y = filter(x->!occursin(r"date|year",x),names(df))
# df_csum = combine(groupby(df, :year), y .=> sum .=> y)
# # Compute cumulative sum for each row of columns A and B
# df = df_csum[!,Not(:year)]


pos = d[!,Cols(r"date|^(prec)|^(snow_stora)|^(Capi)")]
pos = yrsum(pos)
# calculate the sum of each row
psum = DataFrame(
    possums = [sum(eachrow(pos)[i]) for i in 1:size(pos, 1)],
    year=pos[!,:year]
)

#df = neg
# y = filter(x->!occursin(r"date|year",x),names(df))
# df_csum = combine(groupby(df, :year), y .=> sum .=> y);
# # Compute cumulative sum for each row of columns A and B
# df = df_csum[!,Not(:year)]
# calculate the sum of each row

neg = d[!,Not(Cols(r"^(prec)|^(snow_stora)|^(Capi)"))]
#names(neg)
neg = sum.(yrsum(neg))

# sum.(eachcol(permutedims(neg[!,Not(1)])))|>Plots.bar
# sum.(eachcol(permutedims(pos[!,Not(1)])))|>Plots.bar

nsum = DataFrame(
    negsums = [sum(eachrow(neg)[i]) for i in 1:size(neg, 1)],
    year=neg[!,:year]
)


# xd = hcat(d2[!,:date],pos)
# xd = rename(xd ,1=>"date")|>yrsum

# xn = hcat(d2[!,:date],neg)
# xn = rename(xn ,1=>"date")|>yrsum
# xn = nsum
# @df xn Plots.plot(:year,
# cols(propertynames(xn[!,Not(:year)])),
# legend = :topright, seriestype=:bar)

#ln = Symbol.(filter(x->!occursin("year|YY|",x),names(neg)))
# @df neg bar(:year,cols(propertynames(neg)[2:end-4]),
# bar_position=:stack,orientation = :h)

# @df neg Plots.plot(:year,cols(propertynames(neg)[2:end-4]),
# seriestype=:bar,bar_position=:group,orientation = :h)



bw = innerjoin(psum, nsum, on=:year)
bw.bw = bw[!,:possums] .- bw[!,:negsums]

# select!(bw, [:year, :bw])
# ti=DataFrames.metadata(df)|>only|>last|>basename




##########plotting ########################################
using StatsPlots, Plots.PlotMeasures



fact=.85
plotsize = (1600*fact,800*fact)

ti="water-balance of "*basename(pwd())

# uppercase(basename(pwd())) #capitalize

# set the border thickness and color
#Pkg.add("PlotThemes")
#https://docs.juliaplots.org/latest/generated/plotthemes/
#theme(:dao)
theme(:vibrant)
# using Plots
# showtheme()
ann = map(x->string.(round(x;sigdigits=3))*" [mm]",bw.bw)

# M = (hcat(1:nrow(bw),ann))
# M = (hcat(bw.year,ann))
# D = Dict(M[:, 1] .=> M[:, 2])


# p1 = @df bw Plots.plot(
#     :year,:bw,
#     #annotations =(bw.year, bw.bw, ann, :center),
#     #annotations =(bw.year, bw.bw .-200, ann, :top),
#     annotations =(bw.year, bw.bw, ann, :top),
#     #annotations =(D, :center),
#     legend = false, 
#     seriestype=:bar,
#     xticks = bw.year,
#     xlabel = "",
# #    ylab = "bwvr[mm]",
#     ylabel = "[mm]",
#     title = ti,
#     fillcolor = ifelse.(bw.bw .> 0, "cornflowerblue", "coral2"),
#     size=plotsize)

p1 = @df bw Plots.plot(
    :year,:bw,
    #annotations =(bw.year, bw.bw, ann, :top),
    #annotations = (bw.year,bw.bw,(ann,10,:center,:top,:black)),
    #annotations = (bw.year, bw.bw, ann, 8, :left, :top, :black),
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
    #xrotation = 60);
    left_margin = 10mm,
    bottom_margin = 10mm, 
    #bottom_margin = 10px, 
    xrotation = 45)
    

for i in 1:length(bw.year)
    Plots.annotate!(bw.year[i],bw.bw[i],(ann[i],10,:center,:top,:black))
    println(ann[i]*" added")
end



outname="waba-jl.png"
#Plots.savefig(p1,outname,width=1600,height=800)

writedf("waba-input.yrly",dyr)
writedf("waba-year.yrly",bw)


Plots.savefig(p1,outname)
printstyled("$outname saved!\n",color=:green)

#names(d)



#display(p1)
#######make image ##############
#jsimage="/home/ubu/.julia/JuliaImages/sys_js_ts_plots.so"
# using PackageCompiler
# using DataFrames, StatsPlots, CSV, Dates
# pt,nm=(pwd(),"ts_statsplots.so")
# pts = joinpath(pt,nm)
# create_sysimage([:StatsPlots, :DataFrames, :CSV, :Dates],sysimage_path=pts)
#statsplotimage="/home/ubu/.julia/JuliaImages/ts_statsplots.so"


#data = df[!,:value]
# data = bw[!,:bw]
# for (i, val) in enumerate(data)
#     annotate!(i, val-0.2, text(string(val), :center, 12))
# end
#push!(plots, barplot)
#@df p1 annotate!("lala")
#[2:nrow(dyr),:]
# df = permutedims(dyr,Symbol.(dyr[!,1]))
# df = hcat(df, propertynames(dyr))
# #grouped = groupby(dyr, :year)
# grouped = groupby(df, :x1)
# # create a barplot for each group
# plots = []
# #grouped[2][!,1:3]
# #grouped[2][!,end]|>only

# for group in grouped
#     data = group[!,1]
#     barplot = bar(data, legend=false)
#     for (i, val) in enumerate(data)
#         #annotate!(i, val+0.1, text(string(val), :center, 8))
#         annotate!(i, val-1.2, text(string(val), :center, 12))
#         annotate!(i, val+1.2, text(string(group[!,end]|>only), :center, 12))
#     end
#     push!(plots, barplot)
# end
# # # display the plots
# plot(plots..., layout=(length(plots), 1))



# using DataFrames, Plots
# # create a sample DataFrame
# df = DataFrame(group=["A", "A", "B", "B"], value=[1, 2, 3, 4])
# # group the data by the "group" column
# grouped = groupby(df, :group)
# # create a barplot for each group
# plots = []
# for group in grouped
#     data = group[!,:value]
#     barplot = bar(data, legend=false)
#     for (i, val) in enumerate(data)
#         #annotate!(i, val+0.1, text(string(val), :center, 8))
#         annotate!(i, val-0.2, text(string(val), :center, 12))
#     end
#     push!(plots, barplot)
# end

# # display the plots
# plot(plots..., layout=(length(plots), 1))

# df = DataFrame(x = 1:5, y = rand(5), labels = ["A", "B", "C", "D", "E"])
# # create the barplot

# @df df bar(df.x, df.y, 
# labels = df.labels,
# annotations = (df.x, df.y, text(string(df.y)), :black, :center, :bold))


# # add values inside the bars
# for i in 1:size(df, 1)
#     StatsPlots.annotate!(df.x[i], df.y[i], 
#     text(string(df.y[i]), :black, :center, :bold))
# end

# for (i, val) in enumerate(df[!,1])
#     annotate!(i, val-0.2, text(string(val), :center, 12))
# end

#round(bw.bw[1])

# ann = map(x->string.(round(x;sigdigits=3))*" [mm]",bw.bw)
# annotations =(1:nrow(bw), bw.bw .- 1, ann, :center),
# @df bw Plots.plot(:year,:bw,
# annotations = (1:nrow(bw), 1:nrow(bw), ann, :center),
#     legend = false, 
#     seriestype=:bar,
#     xticks = bw.year,
#     xlabel = "",
#     #ylabel = "bwvr[mm]",
#     ylab = "[mm]",
#     yguideposition = :top,
#     title = ti,
#     fillcolor = ifelse.(bw.bw .> 0, "cornflowerblue", "coral2"),
#     size=plotsize ./ 2)

# using DataFrames, StatsPlots
# dx = DataFrame(x = ["A", "B", "C"], y = [10, 20, 30])
# @df dx bar(:x, :y, 
# annotations = (1:3, dx.y .- 1, string.(dx.y), :center))

# dx|>show
# bw|>show                     



# @df bw bar(:year, :bw, 
# annotations =(1:3, bw.bw .- 1, string.(bw.bw), :center),
# ylab = "[mm]", yguideposition = :top)



# using PlotlyJS
# trace = PlotlyJS.bar(x = ["A", "B", "C"], y = [10, 20, 30])
# layout = PlotlyJS.Layout(
#     yguideposition = :top,    
# yaxis = attr(title = attr(text = "Value", standoff = 20)))
# PlotlyJS.plot(trace, layout)
                        