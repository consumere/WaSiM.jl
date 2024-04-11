# 
function ftplin(df::DataFrame)
        if size(df)[2]!=3
            #throw(@warn "wrong number of columns - using dfp!")
            @warn "wrong number of columns - using dfp!"
            @warn "need :Simulated,:Observed,:Date !"
            display(dfp(df))
            return
        end
        ndf = copy(df)
        # rename!(ndf,3=>"date")
        # @warn "last col renamed! !"
        #reorder
        ndf = hcat(ndf[!,Not(Cols(r"date"i))],ndf[:,Cols(r"date"i)])
        # nm = names(ndf)[end-2] #lastbefore date column (qobs)
        # ##subset DF by value (all positive vals..)
        # ndf = filter(nm => x -> x > 0, ndf)
        #filter([1]=> x -> !any(f -> f(x), (ismissing, isnothing, isnan)), ndf)
        ndf = filter(:date=> x -> !any(f -> f(x), (ismissing, isnothing)), ndf)   
        rename!(ndf, [:Simulated,:Observed,:Date])
        dropmissing!(ndf) ##hmm sketchy..
        overall_pearson_r = cor(ndf[!, :Simulated],ndf[!, :Observed])
        r2 = overall_pearson_r^2
        #nse(simulations, evaluation)
        nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
        kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
        ve = round(vef(ndf[!, :Simulated],ndf[!, :Observed]), digits=2)
        subs = "RSQ: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))\nVE: $ve"
        ti = try
            basename(last(collect(DataFrames.metadata(ndf)))[2])
        catch
        @warn "No basename in metadata!"
            raw""
        end 
        
        p = plot(ndf[!, :Date], ndf[!, :Simulated], 
        title=ti, 
        line=:dash, color=:blue, label="Modeled",
        ylabel="[mm/day]", xlabel="", 
        legend=:topleft)
        plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
        annotate!(
        :topright,
        Plots.text("$subs", 10, :black, :right;family="Computer Modern"))
        return p
    end

    """
    wrapper for kge_read(pwd(),"outjl")
    """
    