using DataFrames, Dates

function waread(x::String, flag::Bool=true)
    if flag
        df = CSV.read(x, DataFrame,missingstring=["-9999"], delim="\t",comment="-",
        silencewarnings=true)
        
        names(df)
        
        occursin(r"[A-z]|-",AsTable(df.YY))
        
        occursin(r"[A-z]|-",x)

#        df = df[.map(x->occursin(r"[A-z]|-",x), df.YY), :]  # filtering rows that contain letters or "-"
        
        select(df, map(x->occursin(r"[A-z]|-",x)))
        
        df = df[.!occursin(r"[A-z]|-", df.YY), :]  # filtering rows that contain letters or "-"
        
        
        
        df."HH" .= 12
        source_col_loc = findfirst(x -> x == "YY", names(df))
        df.date = join.(Ref("-"), string.(df[!, source_col_loc:source_col_loc+2])...)
        select!(df, Not(names(df)[1:4]))
        df."date" = DateTime.(df."date")
        df = setindex!(df, "date")
        return df
    else
        df = CSV.read(x, delim_whitespace=true)
        return df
    end
end

df[!, source_col_loc:source_col_loc+2]|>first

mx=ct("so_a")
x=mx[2]
df = CSV.read(x, DataFrame,missingstring=["-9999"], delim="\t",
skipto=4,
        #comment="[A-z]",
        comment="-",
        silencewarnings=true) |>dropmissing
df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
df = select!(df, Not(names(df)[1:4]))
plotf(df)

#od = select!(df,(names(df)[1:4]))

da = DateTime.(od)



function nse(predictions, targets)
    return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
end
simfile="qgesrcm.c8.2016"

#la = waread(simfile)

xd = readdf(simfile)
#newname = "Basin " * string((last(names(xd))))
#rename!(xd, last(names(xd)) => newname)
#a = applycols(convert, df, select!(xd, Not(:date)), numeric, Any) #!

function nsef(predictions, targets)
    return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
end

function vef(obs, sim)
    return (1 - (sum(abs.(obs .- sim)) / sum(obs)))
end



#finaldf = a
finaldf = df
ve = groupby(finaldf,year.(finaldf.date)) .|> combine(names(finaldf), sum) |> unique
#df[!, :year] = year.(df[!,:date]);
finaldf[!,:year]=year.(finaldf[!,:date]) 
#ve = combine(groupby(finaldf,:year), nrow => :count, :date => maximum, size => :size, :x1 => sum, :x2 => sum)
ve = combine(groupby(finaldf,:year), nrow => :count,  :_4 => sum)
#a = DataFrame([parse(Float64, s) for s in col] for col in eachcol(xd)) 

ve.year .= year.(ve.date)
rename!(ve, :year => :year, Symbol(newname) => "sim")
ve.sim ./= ve[:, end-1] .* 100
@show ve

out = nse(predictions=finaldf[newname], targets=finaldf[:, end-1])
println("NSE from $simfile:\nPredictions: $newname\nTarget: $(names(finaldf)[end-1])\nNSE: $out")
println("Yearly Grouped NSE: $(@sprintf("%.2f", nsef(targets=ve[:, end-2], predictions=ve[:, end-1])))")
println("Volume Efficiency: $(@sprintf("%.2f", vef(obs=finaldf[:, end-1], sim=finaldf[newname])))")

overall_pearson_r = cor(finaldf, pairwise=true)
@show overall_pearson_r
