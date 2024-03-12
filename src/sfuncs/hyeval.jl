# 
function hyeval(
        x::Union{Regex,DataFrame}; 
        yscale::Symbol = :log,
        fun::Function = sum, 
        freq::String="monthly",
        ylab::String="[mm/$freq]",
        simcol=1,obscol=2)
        if x isa Regex
            x = first(dfonly(x))
            ndf = waread2(x;silencewarnings=true)
        else
            ndf = reorder_df(x)
        end
                
        ndf = select(ndf,Cols(simcol,obscol,r"date|month|year"))
        dropmissing!(ndf)

        @info "aggregation (sum) freq: $freq"
        try 
        # Resample the DataFrame based on the specified freq
            if freq == "daily" || freq == "D" || freq == "day"
                #ndf = ndf[hour.(ndf.date) .== 0, :]
                ndf = ndf
            elseif freq == "monthly" || freq == "M" || freq == "mon" || freq == "month"
                #df[day.(df.date) .== 1, :]
                ndf = qrtr(ndf;fun=fun,agg=month) #monsum(ndf)
            elseif freq == "quarterly" || freq == "Q" 
                #df[day.(df.date) .== 1 && month.(df.date) % 3 == 1, :]
                ndf = qrtr(ndf;fun=fun)
            elseif freq == "seasonal"
                #df[day.(df.date) .== 1 && month.(df.date) % 3 == 1, :]
                ndf = qrtr(ndf;fun=fun)
            elseif freq == "yearly" || freq == "Y" || freq == "yr" || freq == "year"
                #df[day.(df.date) .== 1 && month.(df.date) == 1, :]
                ndf = qrtr(ndf;fun=fun,agg=year)  #yrsum(ndf)
                #DataFrames.combine(groupby(df, year.(df.date)), y .=> mean .=> y);
                #DataFrames.combine(groupby(df, quarterofyear.(df.date)), y .=> mean .=> y);
            end
        catch e
            @error("smth went wrong",e)
            return
        end

        #sim,obs names
        sim,obs = names(ndf[!,Not(Cols(r"date|month|year"))])[1:2]
        ndf = hcat(ndf[!,Not(Cols(r"date|month|year"))],
            ndf[:,Cols(r"date|month|year")])
    
        #printstyled(names(ndf)...,bold=true,color=:red)
        #split(names(ndf)...,' '),bold=true,color=:red)
        printstyled(names(ndf),bold=true,color=:red)
        
        rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
        
        printstyled(" changed to :Simulated,:Observed,:Date ! \n",bold=true,color=:red)
        
        overall_pearson_r = cor(ndf[!, :Observed], 
            ndf[!, :Simulated])
        r2 = overall_pearson_r^2
        #nse(simulations, evaluation)
        nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
        kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
        ve = round(vef(ndf[!, :Simulated],ndf[!, :Observed]), digits=2)
        subs = "RSQ: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))\nVE: $ve"
        

        ti = try 
            first(split(basename(x),"_"))
        catch
            @warn "No basename in metadata!"
            raw""
        end
        
        fr = replace(freq,"ly"=>"")
        p = Plots.plot(title=ti, ylabel=ylab, 
            xlabel="", 
            yaxis = yscale, 
            legend=:topleft)

            #sim,obs
        Plots.plot!(p, ndf[!, :Date], ndf[!, :Simulated], 
                color=:red, label=sim) #label="Modeled")
        Plots.plot!(p, ndf[!, :Date], ndf[!, :Observed], 
        line=:dash, color=:blue, label=obs)
        #label="Observed")
            
        
        Plots.annotate!(
            :topright,
            Plots.text("$subs", 10, :black, :right;
            family="Computer Modern")
        )

        if freq == "quarterly" || freq == "Q" || freq == "qrtr" 
            Plots.xticks!(
                1:4, ["Q1", "Q2", "Q3", "Q4"])
        end
                
        if freq == "monthly" || freq == "M" || freq == "mon" || freq == "month"
            Plots.xticks!(
                1:12, monthabbr.(1:12))
        end
        
        if freq == "yearly" || freq == "Y" || freq == "yr" || freq == "year"
            Plots.xticks!(
                #1:nrow(ndf), 
                ndf.Date, 
                string.(ndf.Date))
        end
        

        return p
    end

    """
    also for pest outputfiles
    """
    