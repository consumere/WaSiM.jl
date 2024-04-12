#Climate Diagram
using Plots
using StatsPlots
using DataFrames

pwd()
prec = waread(r"prec")
temp = waread(r"temprc")
prec = monsum(prec)
temp = monmean(temp)

pyplot()
p1 = Plots.bar(prec.month, prec.C160, color=:cornflowerblue, 
    xlabel="Months", xflip=false,
    ylabel="Precipitation [mm]", legend=false, yflip=true);
for i in prec.month
    val = round(prec.C160[i], digits=1)
    annotate!(i, prec.C160[i], text("$(val)",8, :bottom))
end
p2 = twinx(p1);
plot!(p2, temp.C160, xlabel="", 
    ylabel="Temperature [°C]", color=:coral2,
    label=false, linewidth=3);
# # Add annotations for temperature values
# for i in 1:12
#     val = round(temp.C160[i], digits=1)
#     annotate!(i, temp.C160[i], text("$(val)",8, :center))
# end
display(p1)


nm="C160"

##DataFrame -> Vector
temp.C160 == vec(Matrix(select!(temp, nm)))

function climateplot(temp::Regex,prec::Regex,col::AbstractString)
    """
    col = :C160 of interest
    """
    col = col
    prec = waread(prec)
    yrs = year.(prec.date)|>unique|>length
    prec = monsum(prec)
    precvec = vec(Matrix(select(prec, col)))
    precvec = precvec ./ yrs
    
    temp = waread(temp)|>monmean
    tempvec = vec(Matrix(select(temp, col)))
    month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    p1 = Plots.bar(prec.month, precvec, color=:cornflowerblue, 
        xlabel="", 
        #xlabel="Months", 
        xflip=false,
        ylabel="Precipitation [mm]", 
        legend=false, yflip=true);
    xticks!(1:12, month_abbr)
    for i in prec.month
        val = round(precvec[i]; digits=1)
        annotate!(i, precvec[i], text("$(val)",8, :bottom))
    end
    p2 = twinx();
    #ann2 = map(x->string.(round(x;sigdigits=0)),tempvec)
    plot!(p2, tempvec, xlabel="", 
        ylabel="Temperature [°C]", color=:coral2,
        #annotations = (temp.month,tempvec, ann2, :center),
        label=false, linewidth=3);
    # # Add annotations for temperature values
    for i in 1:length(tempvec)
        val = round(tempvec[i]; sigdigits=1)
        annotate!(i, tempvec[i], text("$(val)",7, :center))
    end
    display(p1)
end

"D:/Wasim/regio/out/rc200/x12/loc2/"|>cd
climateplot(r"temprc",r"prec_re","tot_average")
climateplot(r"temprc",r"prec_re","C200")

r = r"temprc"
temp = waread(r)|>monmean
tempvec = vec(Matrix(select(temp, col)))
p2 = plot(tempvec);
for i in 1:length(tempvec)
    val = round(tempvec[i]; sigdigits=1)
    annotate!(i, tempvec[i], text("$(val)",7, :center))
end
display(p2)


ms = r"prec_re"|>glob|>first|>dfr|>monsum

prec = waread(r"prec_re")
year.(prec.date)|>unique|>length


vec(Matrix(select!(ms, nm)))



# Example data
temp = [6.9, 8.1, 11.6, 14.7, 19.1, 23.0, 25.6, 25.3, 22.1, 16.5, 11.1, 7.7]
prec = [71.2, 53.5, 65.7, 65.5, 59.5, 41.7, 11.5, 15.9, 25.5, 54.1, 80.7, 82.8]
months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
data = DataFrame(month=months, temp=temp, prec=prec)

@df data begin
    # Create the bar chart for precipitation
    Plots.bar(:month, :prec, color=:blue, xlabel="Month", ylabel="Precipitation (mm)",
        xticks=(1:12, ["" "Jan" "" "Mar" "" "May" "" "Jul" "" "Sep" "" "Nov" ""]),
        xlims=(0.5, 12.5), ylims=(0, maximum(:prec)+50), legend=false)
        
    # Create the line chart for temperature
    plot!(:month, :temp, color=:red, ylabel="Temperature (°C)",
        yticks=(-20:5:40), ylims=(-20, 40))
end

# Add a title to the plot
title!("Climate Diagram")

using Plots

# Generate some example data
month = collect(1:12)
temp = [4.9, 6.1, 8.8, 11.7, 15.4, 18.4, 20.7, 20.5, 17.8, 13.5, 9.0, 5.7]
precip = [60.7, 58.3, 73.5, 87.2, 122.5, 138.1, 112.7, 112.1, 86.5, 76.7, 70.5, 68.1]

# Create the plot
p1 = Plots.plot(month, temp, xlabel="Month", ylabel="Temperature (°C)", label="Temperature", legend=:topleft)
p2 = twinx(p1) # Create a second x-axis for precipitation
plot!(p2, precip, color=:blue, label="Precipitation", xlabel="Precipitation (mm)", ylabel="", legend=:bottomright, yflip=true) # Plot the precipitation data


using DataFrames
using PlotlyJS

pwd()
cd("/mnt/d/temp/saale/output/thulba/v0")
af = filter(x -> occursin(r"^so_", x), readdir(pwd()))


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

#filter(x -> occursin("temper", x), af)

temp = filter(x -> occursin("temperature.", x), af)|>last|>so_read
rain = filter(x -> occursin("precip", x), af)|>only|>so_read
humid = filter(x -> occursin("humid",x), af)|>only|>so_read
# var1 = df1[:, "tas"] # surface air temperature
# var2 = df2[:, "pr"] # precipitation
# var3 = df3[:, "hurs"] # relative humidity
# concatenate the variables into one dataframe
#df = hcat(var1, var2, var3)
#df = hcat(temp, rain, humid)

df = reduce((left, right) -> 
innerjoin(left, right, on = :date,makeunique=true), 
[temp, rain, humid])

z =monmean(df)
# rename the columns
#rename!(df, ["Surface air temperature", "Date", "Precipitation", "Relative humidity","month"])
#df = df[:,Not(:Date)]
# create a figure with three subplots using PlotlyJS
fig = PlotlyJS.plot(df, 
layout=PlotlyJS.Layout(grid=PlotlyJS.attr(rows=3, columns=1)))
# update the figure layout and save it as a html file
PlotlyJS.relayout!(fig, title="Climate plot", height=800)
#savefig(fig, "climate_plot.html")




df = reduce((left, right) -> 
innerjoin(left, right, on = :date,makeunique=true), 
[temp, rain, humid])

z =monmean(df)

rename!(z, 2=>"temp")
rename!(z, 3=>"prec")
rename!(z, 4=>"hum")

@df z begin
    # Create the bar chart for precipitation
    Plots.bar(:month, :prec, color=:blue, xlabel="Month", ylabel="Precipitation (mm)",
        xticks=(1:12, ["" "Jan" "" "Mar" "" "May" "" "Jul" "" "Sep" "" "Nov" ""]),
        xlims=(0.5, 12.5), ylims=(0, maximum(:prec)+50), legend=false)     
    # Create the line chart for temperature
    plot!(:month, :temp, color=:red, ylabel="Temperature (°C)",
        yticks=(-20:5:40), ylims=(-20, 40))
    plot!(:month, :hum, color=:green, ylabel="Relative humidity",
        yticks=(-20:5:40), ylims=(0, 1))
end

# Add a title to the plot
title!("Climate Diagram")










# import libraries
using Plots
using Rasters
using GeoArrays

# load a 3D Raster file as a GeoArray
#ga = GeoArray("path/to/raster.tif")
ga = GeoArrays.read("/mnt/d/temp/saale/saale_25/thulba/saale.dhm")

x="/mnt/d/temp/saale/saale_25/dem.tif"
ga = GeoArrays.read(x)

# get the values, coordinates and CRS of the GeoArray
values = ga.A # a 3D array of raster values
GeoArrays.coords(ga) # a tuple of x, y and band coordinates
crs = ga.crs # a string of CRS definition

#coords =map(parent,GeoArrays.coords(ga))
coords = (1:2616, 1:1382)

GeoArrays.coords(ga)|>last
cds = GeoArrays.coords(ga)|>size
typeof(cds)
typeof(coords)




function rp3(x)
    ga = GeoArrays.read(x)
    values = ga.A # a 3D array of raster values
    GeoArrays.coords(ga) # a tuple of x, y and band coordinates
    #crs = ga.crs # a string of CRS definition
    t = GeoArrays.coords(ga)|>size
    coords = (1:t[1], 1:t[2]) # a Tuple{UnitRange{Int64}, UnitRange{Int64}}
    surface(coords[1], coords[2], values[:, :, 1]) # plot the first band
    xlabel!("x")
    ylabel!("y")
    zlabel!("value")
    #title!("3D Raster Plot")
    ti=basename(x)
    title!(ti)
end


x=("/mnt/d/temp/saale/saale_25/thulba/saale.dhm")
x=("/mnt/d/temp/saale/saale_25/thulba/saale.fzs")
rp3(x)


using PlotlyJS

# List of file names
files = [
    "./so_air_humidity.v0.2017",
    "./so_albedo.v0.2017",
    # Add more file names here
]

files=af
# Iterate over the files
for file in files
    # Read the data from the file (assuming it's a CSV file)
    #data = CSV.read(file)
    data = readdf(file)
    # Extract the required data columns from the DataFrame
    x = data[:, :date]
    y = data[:, 1]
    # Create the climate diagram plot using PlotlyJS
    layout = PlotlyJS.Layout(
        title = "Climate Diagram for $file",
        xaxis = PlotlyJS.attr(title = "X Axis Label"),
        yaxis = PlotlyJS.attr(title = "Y Axis Label"),
    )

    PlotlyJS.plot = PlotlyJS.scatter(x, y, mode = "lines", 
    layout = layout)
    #display(plot)
end



months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
pwd()
temperature1 = readdf(r"so_temp")|>monmean
precipitation1 = readdf(r"so_preci")|>monsum
# # Repeat for each file
temperature2 = readdf(r"temps")|>monmean
precipitation2 = readdf(r"preci")|>monsum
# Repeat for each file
#glob(r"tem")

temperature2 = temperature2[!,[1,end]]

trace1_temp = scatter(x = months, y = temperature1, mode = "lines", name = "Temperature 1")
trace1_precip = bar(x = months, y = precipitation1, name = "Precipitation 1")

trace2_temp = scatter(x = months, y = temperature2, mode = "lines", name = "Temperature 2")
trace2_precip = bar(x = months, y = precipitation2, name = "Precipitation 2")
#PlotlyJS.bar(precipitation1)
bar(x = precipitation2[!,1], y = precipitation2[!,end], name = "Precipitation 2")

layout = PlotlyJS.Layout(title = "Climate Diagram",
                xaxis_title = "Month",
                yaxis_title = "Value")

#PlotlyJS.
plot([trace1_temp, trace1_precip, trace2_temp, trace2_precip], layout)





temperature1 = readdf(r"so_temp")|>monmean
precipitation1 = readdf(r"so_preci")|>monsum
humid1 = readdf(r"humi")|>monsum
humid1 = humid1[!,[1,end]]

z = reduce((left, right) -> 
innerjoin(left, right, on = :month,makeunique=true), 
[temperature1, precipitation1, humid1])

#dcols = propertynames(z)[2]

rename!(z, 2=>"temp")
rename!(z, 3=>"prec")
rename!(z, 4=>"hum")

@df z begin
    # Create the bar chart for precipitation
    Plots.bar(:month, :prec, color=:blue, xlabel="Month", 
    name="Precipitation (mm)",
        xticks=(1:12, ["" "Jan" "" "Mar" "" "May" "" "Jul" "" "Sep" "" "Nov" ""]),
        xlims=(0.5, 12.5), ylims=(0, maximum(:prec)+50), 
        legend=false, yaxis=:lin)     
    # Create the line chart for temperature
    plot!(:month, :temp, color=:red, name="Temperature (°C)",
        yticks=(-20:5:40), ylims=(-20, 35))
    plot!(:month, :hum, color=:green, name="Relative humidity",
        yticks=(-20:5:40), ylims=(0, 1))
end

xlabel!("" "Jan" "" "Mar" "" "May" "" "Jul" "" "Sep" "" "Nov" "")
title!("Climate Diagram")



#using CSV, DataFrames
# Read the data into DataFrames
temperature1 = readdf(r"so_temp")|>monmean
precipitation1 = readdf(r"so_preci")|>monsum
humid1 = readdf(r"humi")|>monmean
select!(humid1, 1, ncol(humid1))

z = reduce((left, right) -> 
innerjoin(left, right, on = :month,makeunique=true), 
[temperature1, precipitation1, humid1])
# Rename the columns
rename!(z, [:month, :temp, :prec, :hum])


gr()
# Plot the data using Plots
bar(z.month, z.prec, color = :blue, xlabel = "Month", ylabel = "Precipitation (mm)",
xticks=(1:12, ["" "Jan" "" "Mar" "" "May" "" "Jul" "" "Sep" "" "Nov" ""]),
    xlims = (0.5, 12.5), ylims = (0, maximum(z.prec) + 50), legend = false)

plot!(z.month, z.temp, color = :red, ylabel = "Temperature (°C)",
    yticks = (-20:5:40), ylims = (-20, 40))

plot!(z.month, z.hum, color = :green, ylabel = "Relative humidity",
    yticks = (-20:5:40), ylims = (0, 1))


plotlyjs()

z
# Extract the variables of interest

using PlotlyJS.jl
temp = z[:, :temp] # temperature column
prec = z[:, :prec] # precipitation column

trace1 = PlotlyJS.scatter(x = 1:12, y = temp, yaxis = "y1", 
    name = "Temperature (°C)")
trace2 = PlotlyJS.bar(x = 1:12, y = prec, yaxis = "y2", 
    name = "Precipitation (mm)")
layout = PlotlyJS.Layout(
    title = "Climate diagram",
    showlegend=false,
                xaxis_title = "Month",
                yaxis_title = "Temperature (°C)",
                yaxis2_title = "Precipitation (mm)",
                yaxis2_overlaying = "y",
                yaxis2_side = "right";template="seaborn")
fig = PlotlyJS.plot([trace2, trace1], layout)

##############################


tdifnc()



###################################################

cd("/mnt/d/temp/saale/output/v7")
#af = filter(x -> occursin(r"^so_", x), readdir(pwd()))
af = filter(x -> occursin(r"temp|prec", x), readdir(pwd()))

temp = filter(x -> occursin("temperature.", x), af)|>last|>so_read
rain = filter(x -> occursin("precip", x), af)|>last|>so_read
sb = filter(x -> occursin(r"sb05", x), readdir(pwd()))|>second|>so_read
df = reduce((left, right) -> 
innerjoin(left, right, on = :date,makeunique=true), 
[temp, rain, sb])

temp = filter(x -> occursin("temperature.", x), af)|>last|>so_read|>monmean
rain = filter(x -> occursin("precip", x), af)|>last|>so_read|>monsum
sb = filter(x -> occursin(r"sb05", x), readdir(pwd()))|>second|>so_read|>monmean

z = reduce((left, right) -> 
innerjoin(left, right, on = :month,makeunique=true), 
[temp, rain, sb])
rename!(z, 2=>"temp")
rename!(z, 3=>"prec")
rename!(z, 4=>"soil_moisture")

temp = z[:, :temp] # temperature column
prec = z[:, :prec] # precipitation column
sm = z[:, :soil_moisture] # rel. soil_moisture

trace1 = PlotlyJS.scatter(x = months[1:12], y = temp, yaxis = "y2", 
    name = "Temperature (°C)")
trace2 = PlotlyJS.bar(x = months[1:12], y = prec, yaxis = "y1", 
    name = "Precipitation (mm)")
trace3 = PlotlyJS.scatter(x = months[1:12], y = sm, yaxis = "log", 
    name = "rel. soil_moisture (-)")

layout = PlotlyJS.Layout(
    #title = dirname(pwd())
    title = split(pwd(),"/")|>last,
    #showlegend=false,
    showlegend=true,
                xaxis_title = " ",
                yaxis2_title = "Temperature (°C)",
                yaxis_title = "Precipitation (mm)",
                yaxis2_overlaying = "y",
                yaxis2_side = "right";template="seaborn")
fig = PlotlyJS.plot([trace2, trace1,trace3], layout)

PlotlyJS.savefig(fig,"tst.png",scale=2)





#sm / T plot

trace1 = PlotlyJS.scatter(x = months[1:12], y = temp, yaxis = "y2", 
    name = "Temperature (°C)")
# trace2 = PlotlyJS.bar(x = months[1:12], y = prec, yaxis = "y1", 
#     name = "Precipitation (mm)")
trace3 = PlotlyJS.scatter(x = months[1:12], y = sm, yaxis = "log", 
    name = "rel. soil_moisture (-)")

layout = PlotlyJS.Layout(
    #title = dirname(pwd())
    title = split(pwd(),"/")|>last,
    #showlegend=false,
    showlegend=true,
                xaxis_title = " ",
                yaxis2_title = "Temperature (°C)",
                yaxis_title = "rel. soil_moisture (-)",
                yaxis2_overlaying = "y",
                yaxis2_side = "right";template="seaborn")
fig = PlotlyJS.plot([trace1,trace3], layout)

PlotlyJS.savefig(fig,"tst.png",scale=2)



# Create the plot
p1 = Plots.plot(months[1:12], temp, xlabel="Month", 
ylabel="Temperature (°C)", label="Temperature", legend=:topleft)
p2 = twinx(p1); # Create a second x-axis for precipitation
bar!(p2, prec, color=:blue, label="Precipitation", 
xlabel="Precipitation (mm)", xflip=false,
ylabel="", legend=:bottomright, yflip=true) # Plot the precipitation data




##T/P with flipped axis

trace1 = PlotlyJS.scatter(x = months[1:12], y = temp, yaxis = "y2", 
    name = "Temperature (°C)")
trace2 = PlotlyJS.bar(x = months[1:12], y = prec, yaxis = "y1", 
    name = "Precipitation (mm)")
layout = PlotlyJS.Layout(
    title = "Climate diagram",
    xaxis_title = "Month",
    yaxis_title = "Precipitation (mm)",
    yaxis2_title = "Temperature (°C)",
    yaxis2_overlaying = "y",
    yaxis2_side = "right",
    template="seaborn",
    showlegend=false,
    yaxis_autorange="reversed",
)
fig = PlotlyJS.plot([trace1,trace2], layout)

#Layout(yaxis_autorange="reversed"),



#########again
"/mnt/d/temp/saale/output/thulba/v2"|>cd
"/mnt/d/temp/saale/output/thulba/v3"|>cd
af = filter(x -> occursin(r"^so_", x), readdir(pwd()))
temp = filter(x -> occursin("temperature.", x), af)|>last|>so_read|>monmean
rain = filter(x -> occursin("precip", x), af)|>last|>so_read|>monsum
sb = filter(x -> occursin(r"sb05", x), readdir(pwd()))|>second|>so_read|>monmean
z = reduce((left, right) -> 
innerjoin(left, right, on = :month,makeunique=true), 
[temp, rain, sb])
rename!(z, 2=>"temp")
rename!(z, 3=>"prec")
rename!(z, 4=>"soil_moisture")

temp = z[:, :temp] # temperature column
prec = z[:, :prec] # precipitation column
sm = z[:, :soil_moisture] # rel. soil_moisture

months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

##T/P with flipped axis
begin
    trace1 = PlotlyJS.scatter(x = months[1:12], y = temp, yaxis = "y2", 
        name = "Temperature (°C)")
    trace2 = PlotlyJS.bar(x = months[1:12], y = prec, yaxis = "y1", 
        name = "Precipitation (mm)")
    layout = PlotlyJS.Layout(
        title = "Climate diagram",
        xaxis_title = "Month",
        yaxis_title = "Precipitation (mm)",
        yaxis2_title = "Temperature (°C)",
        yaxis2_overlaying = "y",
        yaxis2_side = "right",
        template="seaborn",
        showlegend=false,
        yaxis_autorange="reversed",
    )
    fig = PlotlyJS.plot([trace1,trace2], layout)
end


dfs=readalldfs(".")

dfpjs(dfs)


f="D:/Wasim/Goldbach/revision/fab150/pestpp/v6/output/pestout/qgkofab150.p1.2018"
DataFrame(CSV.File(f,skipto=10,header=1))

using DelimitedFiles
function anydf(f::String)
    fx=DelimitedFiles.readdlm(f)
    DataFrame(fx,:auto)
end

df = anydf(f)
vgjl("df = filter")

dx = filter( [2]=> x -> !any(f -> f(x), (ismissing, isnothing)), df)

#for col in eachcol(data);println(col);end
# Read the text file
data = CSV.read(f, DataFrame, delim="\t",header=false,types=Float64)
dx = filter([1]=> x -> !any(f -> f(x), (ismissing, isnothing)), data)

df = dx

#map(x->parse.(Int,x),[df[!,1],df[!,2],df[!,3]])
map(x->Int(x),[df[!,1],df[!,2],df[!,3]])
#df[!,1]=map(x ->Int(x),df[!,1])
#for col in eachcol(df)
for i in 1:3
    df[!,i]=map(x ->Int(x),df[!,i])
end
df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
df=df[:,Not(1:4)]
names(df)
dfp(df)



x="D:/temp/saale/output/thulba/v3/qgesth.v3.2017"
function readf(x::String)
    # Read the text file, preserve line 1 as header column
    df = CSV.read(x, DataFrame, delim="\t",header=1,
    silencewarnings=true,
    types=Float64)
    dropmissing!(df)
    # filter numeric only / subsetting on first column
    df = filter([1]=> x -> !any(f -> f(x), (ismissing, isnothing)), df)
    #df = filter([6]=> x -> !any(-9999.0,x), df)
    # replace to missing...inplace doesnt work!
    #df = filter([5]=> x -> !any(f -> f(x),replace(-9999.0 => missing)), df)
    for i in 5:size(df,2)
        #replace!(df[!,i],-9999.0 => missing)
        df[!,i]=replace(df[!,i],-9999.0 => missing)
    end
    
    for i in 5:size(df,2)
        replace!(df[!,i],-9999.0 => missing)
    end
    # map to int for dates
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df=df[:,Not(1:4)]
    metadata!(df, "filename", x, style=:note);
end

df = readf("D:/temp/saale/output/thulba/v3/qgesth.v3.2017")

xd = CSV.read(x, DataFrame, delim="\t",header=1,types=Float64)

#bigfile ~ 36Mb
x="D:/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/tmean_ce.txt"
@time df = readf(x)
#1.004485 seconds (2.02 M allocations: 371.289 MiB, 2.70% gc time, 10.30% compilation time)
#0.920516 seconds (1.89 M allocations: 270.995 MiB, 3.28% gc time)
dfp(df)

@time df = readdf(x)
#0.821027 seconds (404.41 k allocations: 132.287 MiB, 3.25% gc time, 25.90% compilation time)
#0.612215 seconds (36.79 k allocations: 103.911 MiB, 4.46% gc time)

names(df)
dfp(df)


@df df bar(df.date,df.BAMBERG)

bardf(df[!,Cols(4,end)])
xd=yrmean(df[!,Cols(4,end)])
@df xd bar(cols(2),xlabel=cols(1))
@df xd plot(xd.year,cols(2))

xd = yrmean(df[!,Cols(1,end)])
px = @df xd plot(xd.year,cols(2),label=names(xd)[end],legend=:bottomleft);

#counter = 1
for i in 12:18
    #counter += 1 #no need
    #px[counter]
    xd=yrmean(df[!,Cols(i,end)])
    #@df xd plot!(xd.year,cols(2),label=names(xd)[end])|>display
    px=@df xd plot!(xd.year,cols(2),
    label=names(xd)[end])
end
display(px)
#range(5,10)


##how to show increasing yearmen temperatures from station data..
x = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/temp_1970.txt"
df = readf(x)
#df = readdf(x)
xd = yrmean(df[!,Cols(1,end)])
px = @df xd plot(xd.year,cols(2),label=names(xd)[end],
        legend=:outerbottomleft,yaxis=:log10)
        #legend=:outertopleft)
        #legend=:bottomleft);
for i in 20:size(df,2)-1
    xd=yrmean(df[!,Cols(i,end)])
    px=@df xd plot!(xd.year,cols(2),
    label=names(xd)[end])
end
display(px)

describe(df)


x = "D:/Wasim/Tanalys/DEM/Input_V2/meteo/ts-0/pre_ce.txt"
@time df = readf(x) #736cols
#3.737430 seconds (2.02 M allocations: 1.352 GiB, 3.56% gc time)
begin
    xd = yrmean(df[!,Cols(1,end)])
    px = @df xd plot(xd.year,cols(2),label=names(xd)[end],
            legend=:outerbottomleft,yaxis=:log10)
    for i in 2:size(df,2)-700
        xd=yrmean(df[!,Cols(i,end)])
        px=@df xd plot!(xd.year,cols(2),
        label=names(xd)[end])
    end
    display(px)
end


begin
    xd = yrmean(df[!,Cols(1,end)])
    px = @df xd plot(xd.year,cols(2),label=names(xd)[end],
            legend=:outerbottomleft,yaxis=:log10)
    for i in 2:size(df,2)-700
        xd=yrmean(df[!,Cols(i,end)])
        px=@df xd plot!(xd.year,cols(2),
        label=names(xd)[end])
    end
    display(px)
end