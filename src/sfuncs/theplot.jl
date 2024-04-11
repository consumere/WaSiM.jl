# 
function theplot(x::Regex)
        x = first(rglob(x))
        df = DataFrame(CSV.File(x))
        nm=names(df)[end-1] #lastbefore column (qobs)
        ##subset DF by value (all positive vals..)
        df = filter(nm => x -> x > 0, df)
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        ndf = df[!,Not(1:4)]
        rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
        overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
        r2 = overall_pearson_r^2
        #nse(simulations, evaluation)
        nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
        kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
        ve = round(vef(ndf[!, :Simulated],ndf[!, :Observed]), digits=2)
        subs = "RSQ: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))\nVE: $ve"
        ti = first(split(basename(x),"_"))
        p = Plots.plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
        Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
        Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
        Plots.annotate!(
        :bottomright,
        Plots.text("$subs", 10, :black, :right;family="Computer Modern"))
        return p
    end

    ftp(z::AbstractString) = theplot(first(
        filter(x->occursin(Regex(z,"i"),x),
        filter(x->endswith(x,"qoutjl"),readdir())))
        )

    