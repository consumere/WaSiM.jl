#The Budyko curve is a graphical representation of the 
#long-term water balance of a catchment. 
# It is typically represented with the ratio of 
# actual evapotranspiration to precipitation (ET/P) 
# on the y-axis and the ratio of 
# potential evapotranspiration 
# to precipitation (PET/P) on the x-axis.

cd("D:/Wasim/regio/rcm200/v12/")
"D:/Wasim/regio/out/rc200/x22/f19/"|>cd
et = dfr(r"^evar")
ep = dfr(r"^evap")
pre = dfr(r"^prec")
dfs = map(yrsum,[et,ep,pre])
dm = mall(dfs;xcol="year")
#xm = select(dm,r"tot")
#xm = @rsubset dm :year==2016
# rename!(xm,["ET","PET","P"])
# subcatchments = xm

#stack(xm)
ET = stack(yrsum(et),value_name=:ET)
PET = stack(yrsum(ep),value_name=:PET)
P = stack(yrsum(pre),value_name=:P)

subcatchments = hcat(ET[!,end],PET[!,end],P[!,end])
sc = DataFrame(subcatchments,:auto)
subcatchments = rename(sc,["ET","PET","P"])

using Plots
# Assuming `subcatchments` is a DataFrame with columns `ET`, `PET`, and `P`
subcatchments.ET_P = subcatchments.ET ./ subcatchments.P
subcatchments.PET_P = subcatchments.PET ./ subcatchments.P

begin
    # Add the Budyko curve
    #x = 0:0.01:2.5
    x = sort(subcatchments.PET ./ subcatchments.P)
    y = x ./ (1 .+ x) # Budyko curve
    plot!(x, y, label = "Budyko curve", color = :black)
    scatter(subcatchments.PET_P, subcatchments.ET_P, label = "Subcatchments")
    plot!(x, y, label = "Budyko curve", color = :black)
    xlabel!("PET/P")
    ylabel!("ET/P")
    title!("Budyko Space")
end

"""
after (Fuh,1981)
ω is a parameter that controls the shape of the curve.
"""
function bcurve(pet,p,et,ω)
    return 1+(pet/p)-(1+(pet/p)^ω)^(-1/ω)
end

begin 
    #x = 0:0.01:.69
    #Aridity Index PET/P
    x = sort(subcatchments.PET ./ subcatchments.P)
    #x = 0:0.01:(length(y)-1)
    y = bcurve.(subcatchments.PET,subcatchments.P,
        subcatchments.ET,1.5)
    sort!(y)
    scatter(subcatchments.PET_P, subcatchments.ET_P, label = "Subcatchments")
    plot!(x, y, label = "Budyko curve", color = :black)
    xlabel!("PET/P")
    ylabel!("ET/P")
    title!("Budyko Space after Paper")
end

# How to KruskalWallisTest
using HypothesisTests, DataFrames, RDatasets
airquality = dataset("datasets", "airquality")
dropmissing!(airquality)
# Assuming `airquality` is your DataFrame
groups = groupby(airquality, :Month)
# Get the Ozone values for each month
ozone_values = [group.Ozone for group in groups]
# Perform the Kruskal-Wallis test
result = HypothesisTests.KruskalWallisTest(ozone_values...)
println(result)

#Kolmogorov-Smirnov test
#solar_values = [group.Solar.R for group in groups]
solar_values = [vec(Matrix(select(group,2))) for group in groups]
result = HypothesisTests.ApproximateTwoSampleKSTest(
    ozone_values[1],solar_values[1])
ApproximateTwoSampleKSTest(
        ozone_values[1],ozone_values[2])

#If you want to test whether the `ozone_values` 
#and `solar_values` come from the same distribution, 
#you can use the two-sample Kolmogorov-Smirnov test. This test is nonparametric and does not assume any specific distribution of the data.

#The Kolmogorov-Smirnov test compares the empirical distribution functions of two samples to see if they come from the same distribution. The null hypothesis is that the two samples come from the same distribution.
p_val = []
for (ozone, solar) in zip(ozone_values, solar_values)
    result = ApproximateTwoSampleKSTest(ozone, solar)
    println(result)
    push!(p_val, pvalue(result))
end
println("The p-values of the Kolmogorov-Smirnov test are: ", 
    pval)
#[f for f in p_val]

# using HypothesisTests
# # Assuming `df` is your DataFrame
# ks_test_result = KruskalWallisTest(df.ths, df.thr)
# # Get the p-value
# p_val = pvalue(ks_test_result)

# using RDatasets
# aq = dataset("datasets", "airquality")
# # KruskalWallisTest("Ozone","Month",aq
# # #BoxPierceTest(y, lag, dof=0)
# # parse.(Float32, out2.ksat)
# # BoxPierceTest(out2.ksat, 1, 0)
