#
#lk="D:/Wasim/Goldbach/stateini/evaporation_modis_500/2022/ctlcum"
lk="/mnt/d/Wasim/Goldbach/stateini/evaporation_modis_500/2022/ctlcum"
using DataFrames, StatsPlots, GLM, CSV, Dates

function readdf(x::AbstractString)
    "--- main reader ---"
    ms=["-9999","lin","log"]
    df::DataFrame = CSV.read(x,DataFrame,missingstring=ms,
    types = Float64,
    delim="\t",
    silencewarnings=false,
    normalizenames=true,
    drop=(i, nm) -> i == 4) |> dropmissing
    df.YY=map(x ->Int(x),df.YY);
    df.MM=map(x ->Int(x),df.MM);
    df.DD=map(x ->Int(x),df.DD);
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
    DataFrames.metadata!(df, "filename", x, style=:note);
end

nd=readdf(lk)
df = rename(nd,1=>"modis",2=>"wasim")

model = lm(@formula(modis ~ wasim), df)
@df df scatter(:modis , :wasim);
plot!(model)

data = df[!,Not(:date)]
rename!(data,1=>"x",2=>"y")
r = lm( @formula(y ~ x), data)
pred = predict(r, data, interval = :confidence, level = 0.95)
p = @df data scatter(:x, :y, leg = false)
# sort data on x
pred_s = pred[sortperm(data[!,:x]), : ]
x_s = sort(data[!, :x])
plot!(p, x_s, pred_s.prediction, linewidth = 2,
        ribbon = (pred_s.prediction .- pred_s.lower, pred_s.upper .- pred_s.prediction))

using LinearFitXYerrors, DataFrames
st = linearfitxy(data.x, data.y, isplot=true, ratio=:auto)



linearfitxy(df.modis, df.wasim, isplot=true, ratio=1)
plot_covariance_ellipses!

linearfitxy(data.x, data.y; r=0,isplot=true,ratio=:auto)

px=linearfitxy(df.modis, df.wasim, isplot=false, ratio=:auto)
LinearFitXYerrors.plot_linfitxy(px)
LinearFitXYerrors.plot_covariance_ellipses!(px)

outname="xy-plot.jpg"
Plots.savefig(LinearFitXYerrors.plot_linfitxy(st; ratio=:auto),outname)

#save_linfitxy(linearfitxy(df.modis, df.wasim, isplot=false, ratio=:auto), outname)


lnk="/mnt/d/Wasim/Goldbach/stateini/evaporation_modis_500/2022/wasim_etr_modis_v45.dat"
df = readdf(lnk)
describe(df)
linearfitxy(df.wasim, df.modis, isplot=true, ratio=:auto)
