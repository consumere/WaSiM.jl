using CSV, DataFrames, Statistics

# Load the data from the two files into two separate data frames
# df1 = CSV.read("file1.csv", dateformat="yyyy-mm-dd", DataFrame)
# df2 = CSV.read("file2.csv", dateformat="yyyy-mm-dd", DataFrame)
f1="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
f2="/mnt/d/Wasim/regio/out/rc200/v2/qgesrcm.v2.2012"
df1 = waread(f1)
df2 = waread(f2)

# Remove any columns that are not relevant to your analysis
# df1 = select!(df1, Not(:cat_var1, :cat_var2, :missing_var))
# df2 = select!(df2, Not(:cat_var1, :cat_var2, :missing_var))

df = innerjoin(df1, df2, on = :date,makeunique=true)


merged_df = df[:,Not(:date)]#|>dropmissing
# df = reduce((left, right) -> 
#       innerjoin(left, right, on = :date,makeunique=true), 
#       dfs)

#dropmissing(merged_df,1,disallowmissing=true)
#dropmissing(merged_df,1)
#,names(merged_df))

merged_df = merged_df[:,2:end]
# Calculate the pairwise correlation coefficients and R² values between each column in the merged data frame
results = [(names(merged_df)[i], names(merged_df)[j], cor(merged_df[:, i], merged_df[:, j]), cor(merged_df[:, i], merged_df[:, j])^2) for i in 2:size(merged_df, 2), j in 2:size(merged_df, 2)]

# Drop the rows containing NaN values
# results = filter(row -> !any(isnan.(row)), results)
#filter(row -> !any(ismissing.(row)), results[5:10])
results = filter(row -> !any(ismissing.(row)), results)

#map(x -> coalesce(x, missing), results)
#x = (1.0, NaN, 2.0, NaN, 3.0)
#(1.0, missing, 2.0, missing, 3.0)
#y = map(x -> replace(isnan(x), missing), x)

#results = filter(row -> !any(isnan.(row)), results)




# Identify the highest correlation coefficient and the corresponding column names in each data frame
# Create a DataFrame to store the results
output_df = DataFrame(column_1=String[], column_2=String[], correlation_coefficient=Float64[], r_squared=Float64[])
# Loop through the results and add them to the output DataFrame
for result in results
    push!(output_df, [result[1], result[2], result[3], result[4]])
end

# Print the output DataFrame
display(output_df)

filtered_df = output_df[0.6 .<= output_df.r_squared .<= 0.9, :]


filtered_df
writedf("corrmat.txt",filtered_df)




using DataFrames, Plots
# create a random DataFrame with 4 columns
df = filtered_df

@df df corrplot(cols(propertynames(df)[3:end]), grid = true)

vgjl("heatmap")
@df df StatsPlots.heatmap(cols(propertynames(df)[end]))



# calculate the correlation matrix
cm = cor(df[:,end],df[:,end-1])
# plot the correlation matrix as a heatmap
heatmap(cm, xticks = 1:4, yticks = 1:4, colorbar_title = "correlation")



# Identify the highest correlation coefficient and the corresponding column names in each data frame
max_r_squared = -Inf
max_corr_col1 = []
max_corr_col2 = []
for result in results
    if abs(result[3]) > abs(max_r_squared)
        push!(max_r_squared,result[4])
        push!(max_corr_col1,result[1])
        push!(max_corr_col2,result[2])
    end
end


# Get the data types of each column
# col_types = map(x->eltype(x),)

# Cols(df)|>eltype
# eltype(df[:,1])

# Cols(merged_df)|>typeof

col_types = []
for i in eachcol(df)
push!(col_types,eltype(i))
end

# Calculate the pairwise correlation coefficients and R² values between each column in the merged data frame
results = Tuple{String, String, Float64, Float64}[]
for i in 2:size(merged_df, 2), j in 2:size(merged_df, 2)
    # Check if both columns are numeric before computing the correlation coefficient
    if (col_types[i] <: Number) && (col_types[j] <: Number)
        # Calculate the correlation coefficient and R² value
        r = cor(merged_df[:, i], merged_df[:, j])
        r2 = r^2
        # Add the result to the output array
        push!(results, (names(merged_df)[i], names(merged_df)[j], r, r2))
    end
end

# Drop the rows containing NaN values
results = filter(row -> !any(isnan.(row)), results)

# Create a DataFrame to store the results
output_df = DataFrame(column_1=String[], column_2=String[], correlation_coefficient=Float64[], r_squared=Float64[])

# Loop through the results and add them to the output DataFrame
for result in results
    push!(output_df, [result[1], result[2], result[3], result[4]])
end

# Print the output DataFrame
display(output_df)



f1="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/spec_main_12.txt"
f2="/mnt/d/Wasim/regio/out/rc200/v2/qgesrcm.v2.2012"
df1 = waread(f1)
df2 = waread(f2)
df = innerjoin(df1, df2, on = :date,makeunique=true)
merged_df = df[:,Not(:date)]
merged_df |>dropmissing!
results = [(names(merged_df)[i], names(merged_df)[j], cor(merged_df[:, i], merged_df[:, j]), cor(merged_df[:, i], merged_df[:, j])^2) for i in 2:size(merged_df, 2), j in 2:size(merged_df, 2)]
output_df = DataFrame(Col1=String[], Col2=String[], correlation_coefficient=Float64[], r_squared=Float64[])
# Loop through the results and add them to the output DataFrame
for result in results
    push!(output_df, [result[1], result[2], result[3], result[4]])
end
# Print the output DataFrame
display(output_df)
filtered_df = output_df[0.6 .<= output_df.r_squared .<= 0.95, :]
filtered_df
writedf("corrmat_specdis.txt",filtered_df)
# create a random DataFrame with 4 columns
df = filtered_df
#@df df corrplot(cols(propertynames(df)[3:end]), grid = true)
dfp("/mnt/d/Wasim/regio/out/rc200/v2/Bruennstadt_qout")
lplot("/mnt/d/Wasim/regio/out/rc200/v2/Wuerzburg_Main_qout")
lplot("/mnt/d/Wasim/regio/out/rc200/v2/Wuerzburg_Main_qout")

pt="/mnt/d/Wasim/regio/out/rc200/v5/Atzhausen_qout"
lplot(pt)
u = readdf(pt)
qqp(u)
"corr"|>vgjl

@df u corrplot(cols(1:2), grid = true)


colnames1 = names(df1)
colnames2 = names(df2)
# calculate the correlation coefficient and R-squared
function cor_coef(x, y) 
    cor(x, y)
end
r_squared(x, y) = cor(x, y)^2

# calculate the correlation coefficients for all possible pairs of columns
cor_coefs = Dict{String, Float64}()
r_squares = Dict{String, Float64}()
for col1 in colnames1[1:end-1]
    for col2 in colnames2[1:end-1]
        x = df1[!, col1]
        y = df2[!, col2]
        cor_coef = cor(x, y)
        cor_coefs["$col1 vs $col2"] = cor_coef
        r_squares["$col1 vs $col2"] = cor_coef^2
    end
end
# create a summary data frame
summary_df = DataFrame(
    Columns = keys(cor_coefs),
    Correlation = values(cor_coefs),
    R_squared = values(r_squares)
)

# Plots.heatmap(
#     df[!,1],df[!,2],df[!,4]
# )


"/mnt/d/Wasim/regio/out/rc200/v4/Bad_Brueckenau_qout"|>dfp
"/mnt/d/Wasim/regio/out/rc200/v4/Bad_Brueckenau_qout"|>lplot
"/mnt/d/Wasim/regio/out/rc200/v4/Scheinfeld_qout"|>lplot



using DataFrames, Plots
# create a sample DataFrame with two columns and a date column
df = DataFrame(Date = Date.(["2022-01-01", "2022-01-02", "2022-01-03", "2022-01-04"]),
               Simulated = [1.0, 2.0, 3.0, 4.0],
               Observed = [1.1, 1.8, 3.2, 3.9])

# create a time series plot of the two columns
plot(df.Date, [df.Simulated, df.Observed], label=["Simulated" "Observed"], xlabel="Date", ylabel="Value")
# calculate the R2 score and display it on the plot
r2 = round(cor(df.Simulated, df.Observed)^2, digits=2)
annotate!(last(df.Date), maximum(df.Observed), text("R² = $r2", :black, :right))

x="/mnt/d/Wasim/regio/out/rc200/v4/Scheinfeld_qout"

function dpr(x)
    df=readdf(x)
    v = (names(df[!,1:2]))
    a = reshape(v, 1, 2)
    Plots.plot(df.date,[df[!,2], df[!,1]], 
    #label=["Simulated" "Observed"], xlabel="Date", ylabel="Value")
    label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)
    r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
    annotate!(last(df.date), 0.85*maximum(df[!,1]),
    text("R² = $r2", 10, :black, :right))
    kge = round(kge2(df[!,1], df[!,2]), digits=2)
    annotate!(
        #:topright,
        last(df.date), 0.95*maximum(df[!,1]),
    text("KGE = $kge", 10, :black, :right))
end

dpr("/mnt/d/Wasim/regio/out/rc200/v4/Bad_Brueckenau_qout")

glob(r"qout")|>first|>dpr
glob(r"qout")|>second|>dpr

hcat(names(df[!,1:2]))
["Simulated" "Observed"]|>typeof

glob(r"qout")|>first|>dfp
glob(r"qout")|>second|>dpr!
#even
dpr(r"qout")


"/mnt/d/Wasim/regio/out/rc200/r2"|>cd
v=rglob("rad")
v|>second|>rp



#"D:\Wasim\regio\rcm\ez2.shp"
########## Rasters cut--- > geht nich gut in Julia #####################
using Shapefile
shp = "/mnt/d/Wasim/regio/rcm200/ezg.shp"
ezg = Shapefile.Handle(shp)
ezg.crs

pt="/mnt/d/Wasim/regio/out/rc200/v2/precrcm.sum.nc"
#r = Raster(pt,missingval=-9999)
pt="/mnt/d/Wasim/regio/out/rc200/v4/rainrcm.v4.2010"
r = readras(pt)

rp(r)
ezg.shapes[1:5]|>plot!

ezg.shapes|>size
msk=mask_trim(r, ezg.shapes[120:140])
rp(msk)


mask_trim(raster, poly) = Rasters.trim(Rasters.mask(raster; with=poly); pad=10)
msk=mask_trim(r, ezg.shapes[1:2])

msk=mask_trim(r, ezg.shapes[5:50])
msk=mask_trim(r, ezg.shapes)
msk=Rasters.trim(Rasters.mask(r;with=ezg.shapes[70:end]))
plot(msk)

plot!(r)
ezg.shapes[1:5]|>plot!


########## Rasters cut--- > geht nich gut in Julia #####################
using Shapefile
shp = "/mnt/d/Wasim/regio/rcm/ez2.shp"
ezg = Shapefile.Handle(shp)
ezg.crs
ezg|>plot
annotate!(ezg, text(1:length(ezg.shapes)))
"annotate!"|>vgjl

r = readras("/mnt/d/Wasim/regio/out/rc200/v5/station2/temprcm.2012.nc")

ezg.shapes|>size
msk=mask_trim(r, ezg.shapes[2:10])
rp(msk)


mask_trim


