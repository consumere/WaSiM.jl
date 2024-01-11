module wajs

    using PlotlyJS
    using DataFrames, CSV, Statistics, Dates, StatsPlots, Distributions
    using DelimitedFiles, Grep, Printf

    function dfyrs(df::DataFrame;logy=true)
        ti = DataFrames.metadata(df)|>only|>last|>basename
        fact,logy = 1,0
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :year] = year.(df[!,:date]);
        df = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        #df = df[!,Not("date")]
        ln = Symbol.(filter(x->!occursin("year",x),names(df)))
        nrows=size(df)[2]-1
        if nrows == 1
            ln = only(ln)
            fig = 
            PlotlyJS.plot(
            PlotlyJS.scatter(x=df.year, y=df[!,ln],
            name=ln,type="bar")
            );
            PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
        else
            fig = PlotlyJS.make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in ln
                PlotlyJS.add_trace!(fig, 
                PlotlyJS.scatter(x=df.year, y=df[:,i],
                name=i));
            end
            if logy == true
                PlotlyJS.relayout!(fig,yaxis_type="log",
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            else
                PlotlyJS.relayout!(fig,
                height=600*fact,width=900*fact,
                title_text="Series of "*ti)
            end
        end
        display(fig)
    end

    function dfpjsold(df::DataFrame;)
        nrows=size(df)[2]-1 
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact,logy = 1,0
        if logy == true
            PlotlyJS.relayout!(fig,yaxis_type="log",
            height=600*fact,width=900*fact,
            title_text="Series of "*ti,
            xaxis=PlotlyJS.attr(    rangeslider_visible=true,rangeselector=PlotlyJS.attr(
            buttons=[
                PlotlyJS.attr(count=1, label="1m", step="month", stepmode="backward"),
                PlotlyJS.attr(count=6, label="6m", step="month", stepmode="backward"),
                PlotlyJS.attr(count=1, label="YTD", step="year", stepmode="todate"),
                PlotlyJS.attr(count=1, label="1y", step="year", stepmode="backward"),
                PlotlyJS.attr(step="all")
            ]    )))
        else
            PlotlyJS.relayout!(fig,
            height=600*fact,width=900*fact,
            title_text="Series of "*ti,
            xaxis=PlotlyJS.attr(rangeslider_visible=true,
        rangeselector=PlotlyJS.attr(        buttons=[
            PlotlyJS.attr(count=1, label="1m", step="month", stepmode="backward"),
                PlotlyJS.attr(count=6, label="6m", step="month", stepmode="backward"),
                PlotlyJS.attr(count=1, label="YTD", step="year", stepmode="todate"),
                PlotlyJS.attr(count=1, label="1y", step="year", stepmode="backward"),
                PlotlyJS.attr(step="all")
            ]    )))
        end
        display(fig)
    end

    function dfpjs(df::String;)
        df = waread(df)

        if names(df)[end]!="date"
            df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        end
        
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end

        s = (filter(x->!occursin(r"year|date",x),names(df)))
        #renamer - remove char _   
        for x in s
            newname=replace(x,"_"=>" ")
            rename!(df,Dict(x=>newname))
        end
        s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))

        ##make scores
        overall_pearson_r = cor(df[!,2], df[!,1])
        r2 = overall_pearson_r^2
        nse_score = nse(df)
        kge_score = kge(df)

        subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"

        nrows=size(df)[2]-1

        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact = .66
        PlotlyJS.relayout!(fig,
            template="seaborn",
            #template="simple_white",
            height=650*fact,
            width=1200*fact,
            title_text=ti,
            xaxis_rangeslider_visible=true,
            annotations=[attr(
                        text=subs,
                        #x=minimum(df[!,2]),
                        x=maximum(df.date),
                        xanchor="right",
                        yanchor="bottom",
                        xref="x",
                        yref="y",
                        showarrow=false,
                        bordercolor="#c7c7c7",
                        borderwidth=2,
                        borderpad=4,
                        bgcolor="#ff7f0e",
                        opacity=0.6
                    )],
            updatemenus=[
                Dict(
                    "type" => "buttons",
                    "direction" => "left",
                    "buttons" => [
                        Dict(
                            "args" => [Dict("yaxis.type" => "linear")],
                            "label" => "Linear Scale",
                            "method" => "relayout"
                        ),
                        Dict(
                            "args" => [Dict("yaxis.type" => "log")],
                            "label" => "Log Scale",
                            "method" => "relayout"
                        )
                    ],
                    "pad" => Dict("r" => 1, "t" => 10),
                    "showactive" => true,
                    "x" => 0.11,
                    "xanchor" => "left",
                    "y" => 1.1,
                    "yanchor" => "auto"
                ),
            ]
            )

        return fig
    end

    function dfpjs(df::Regex;)
        df = waread(df)
        if names(df)[end]!="date"
            df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        end
        
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end

        s = (filter(x->!occursin(r"year|date",x),names(df)))
        #renamer - remove char _   
        for x in s
            newname=replace(x,"_"=>" ")
            rename!(df,Dict(x=>newname))
        end
        s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))

        ##make scores
        overall_pearson_r = cor(df[!,2], df[!,1])
        r2 = overall_pearson_r^2
        nse_score = nse(df)
        kge_score = kge(df)

        subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"

        nrows=size(df)[2]-1

        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact = .66
        PlotlyJS.relayout!(fig,
            template="seaborn",
            #template="simple_white",
            height=650*fact,
            width=1200*fact,
            title_text=ti,
            xaxis_rangeslider_visible=true,
            annotations=[attr(
                        text=subs,
                        #x=minimum(df[!,2]),
                        x=maximum(df.date),
                        xanchor="right",
                        yanchor="bottom",
                        xref="x",
                        yref="y",
                        showarrow=false,
                        bordercolor="#c7c7c7",
                        borderwidth=2,
                        borderpad=4,
                        bgcolor="#ff7f0e",
                        opacity=0.6
                    )],
            updatemenus=[
                Dict(
                    "type" => "buttons",
                    "direction" => "left",
                    "buttons" => [
                        Dict(
                            "args" => [Dict("yaxis.type" => "linear")],
                            "label" => "Linear Scale",
                            "method" => "relayout"
                        ),
                        Dict(
                            "args" => [Dict("yaxis.type" => "log")],
                            "label" => "Log Scale",
                            "method" => "relayout"
                        )
                    ],
                    "pad" => Dict("r" => 1, "t" => 10),
                    "showactive" => true,
                    "x" => 0.11,
                    "xanchor" => "left",
                    "y" => 1.1,
                    "yanchor" => "auto"
                ),
            ]
            )

        return fig
    end

    function dfbarjs(df::Regex;)
        df = waread(df)
        nrows=size(df)[2]-1
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.bar(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact,logy = 0.66,0
        if logy == true
            PlotlyJS.relayout!(fig,
            template="seaborn",
            yaxis_type="log",
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        else
            PlotlyJS.relayout!(fig,
            template="seaborn",
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        end
        display(fig)
    end

    function kge1(simulations, evaluation)
        r = cor(simulations, evaluation)
        α = std(simulations) / std(evaluation)
        β = mean(simulations) / mean(evaluation)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    function tpjs(x::DataFrame)
        """
        theplot optimized to PlotlyJS
        """
        ndf = x
        nm=names(ndf)[2]     #obscol
        ##subset DF by value (all positive vals..)
        ndf = filter(nm => x -> x > 0, ndf)
        rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
        dropmissing!(ndf)
        overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
        r2 = overall_pearson_r^2
        #qpl(ndf)
        nse_score = nse(ndf)
        kge_score = kge(ndf)
        ti = try
            basename(last(collect(DataFrames.metadata(ndf)))[2])
        catch
        @warn "No basename in metadata!"
            raw""
        end 
        subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"
        p = Plots.plot(
        ndf[!, :Date], ndf[!, :Simulated], color=:red, 
        label="Modeled",
        title=ti, ylabel="[mm/day]", xlabel="modeled time", 
        yscale=:log10, 
        legend=:outerbottomleft);
        plot!(p, ndf[!, :Date], ndf[!, :Observed],
        line=:dash, color=:blue, 
        label="Observed")
        annotate!(
            last(ndf.Date), mean(ndf[!,2]),
            text("$subs", 10, :black, :right))
        return p
    end

    function tpjs(x::Regex)
        """
        theplot optimized to PlotlyJS
        """
        ndf = waread(x)
        dropmissing!(ndf)
        #nm=names(ndf)[2]     #obscol
        #ndf = filter(nm => x -> x > 0, ndf)
        rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
        
        overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
        r2 = overall_pearson_r^2
        nse_score = nse(ndf)
        kge_score = kge(ndf)
        ti = try
            basename(last(collect(DataFrames.metadata(ndf)))[2])
        catch
        @warn "No basename in metadata!"
            raw""
        end 
        subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"
        p = Plots.plot(
        ndf[!, :Date], ndf[!, :Simulated], color=:red, 
        label="Modeled",
        title=ti, ylabel="[mm/day]", xlabel="modeled time", 
        yscale=:log10, 
        legend=:outerbottomleft);
        plot!(p, ndf[!, :Date], ndf[!, :Observed],
        line=:dash, color=:blue, 
        label="Observed")
        annotate!(
            last(ndf.Date), mean(ndf[!,2]),
            text("$subs", 10, :black, :right))
        return p
    end

    function tpjs(x::AbstractString)
        """
        theplot optimized to PlotlyJS
        """
        df = DataFrame(CSV.File(x))
        nm=names(df)[end-1] #lastbefore column (qobs)
        ##subset DF by value (all positive vals..)
        df = filter(nm => x -> x > 0, df)
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        ndf = df[!,Not(1:4)]
        rename!(ndf, [:Simulated,:Observed,:Date]) #wie in gof3.r
        overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
        r2 = overall_pearson_r^2
        nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
        kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
        ti = first(split(basename(x),"_")) 
        #subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
        subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"

        p = Plots.plot(
        ndf[!, :Date], ndf[!, :Simulated], color=:red, label="Modeled",
        title=ti, ylabel="[mm/day]", xlabel="modeled time", 
        yscale=:log10, 
        legend=:outerbottomleft)

        #plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
        plot!(p, ndf[!, :Date], ndf[!, :Observed],line=:dash, color=:blue, 
        label="Observed")
        annotate!(
            #:bottomright,
            last(ndf.Date), mean(ndf[!,1]),
            #last(ndf.Date), 0.95*maximum(ndf[!,1]),
            text("$subs", 10, :black, :right))
            #nrow(ndf), 0.95*maximum(ndf.Observed),
            #nrow(ndf),minimum(ndf.Observed),
        return p
    end

    ftpjs = tpjs

    function plotlybaryr(df::DataFrame)
        p = PlotlyJS.plot(df, kind = "bar");
        #rename(df,replace(names(df),"_1"=>""))
        s = (filter(x->!occursin(r"year|date",x),names(df)))
        #renamer - remove chars   
        # for x in s
        #     newname=replace(x,"_1"=>"")
        #     rename!(df,Dict(x=>newname))
        # end
        
        s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
        
        for i in s;
            PlotlyJS.add_trace!(p, 
            PlotlyJS.bar(x=df.year, y=df[:,i],
            name=i)       );
        end
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        PlotlyJS.relayout!(p,
        template="seaborn",
        #height=600*1.5,width=900*1.5, #dirname(pwd()
        #title_text="Series of "*basename(path)
        title_text=ti)
        display(p)
    end

    baryrjs = plotlybaryr

    function dfpjs(df::DataFrame;)

        if names(df)[end]!="date"
            df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        end
        
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end

        s = (filter(x->!occursin(r"year|date",x),names(df)))
        #renamer - remove char _   
        for x in s
            newname=replace(x,"_"=>" ")
            rename!(df,Dict(x=>newname))
        end
        s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))

        ##make scores
        overall_pearson_r = cor(df[!,2], df[!,1])
        r2 = overall_pearson_r^2
        nse_score = nse(df)
        kge_score = kge(df)

        subs = "Pearson R²: $(round(r2, digits=2))<br>NSE: $(round(nse_score, digits=2))<br>KGE: $(round(kge_score, digits=2))"

        fig = PlotlyJS.make_subplots(shared_xaxes=true, shared_yaxes=true)

        for i in s
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:, i], name=i)
            )
        end

        fact = .66
        PlotlyJS.relayout!(fig,
            template="seaborn",
            #template="simple_white",
            height=650*fact,
            width=1200*fact,
            title_text=ti,
            xaxis_rangeslider_visible=true,
            annotations=[attr(
                        text=subs,
                        #x=minimum(df[!,2]),
                        x=maximum(df.date),
                        xanchor="right",
                        yanchor="bottom",
                        xref="x",
                        yref="y",
                        showarrow=false,
                        bordercolor="#c7c7c7",
                        borderwidth=2,
                        borderpad=4,
                        bgcolor="#ff7f0e",
                        opacity=0.6
                    )],
            updatemenus=[
                Dict(
                    "type" => "buttons",
                    "direction" => "left",
                    "buttons" => [
                        Dict(
                            "args" => [Dict("yaxis.type" => "linear")],
                            "label" => "Linear Scale",
                            "method" => "relayout"
                        ),
                        Dict(
                            "args" => [Dict("yaxis.type" => "log")],
                            "label" => "Log Scale",
                            "method" => "relayout"
                        )
                    ],
                    "pad" => Dict("r" => 1, "t" => 10),
                    "showactive" => true,
                    "x" => 0.11,
                    "xanchor" => "left",
                    "y" => 1.1,
                    "yanchor" => "auto"
                ),
            ]
            )

        return fig

    end

    function xdf(df::DataFrame)
        try
            nrows=size(df)[2]-1 
            fig = make_subplots(
                shared_xaxes=true, 
                shared_yaxes=true    
                );
            for i in 1:nrows;
                add_trace!(fig, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=names(df)[i]));
            end
            PlotlyJS.relayout!(fig,
            template="plotly_dark",
            yaxis_type="log")
            display(fig)
        catch e
            println("An error occurred: ", e)
        finally
            println("showing plot...")
            println(describe(df))
    end
    end


    function fdf(df::DataFrame)
        nrows=size(df)[2]-1 
        o = DataFrames.metadata(df)|>collect
        ti = basename(o[1][2])
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact = 0.7
        logy = true;
        if logy == true
            PlotlyJS.relayout!(fig,
            template="seaborn",
            yaxis_type="log",
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        else
            PlotlyJS.relayout!(fig,
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        end
        println(DataFrames.describe(df))
        println("showing plot...")
        display(fig)
    end


    function dfplot(df::DataFrame)
        nrows=size(df)[2]-1
        st=[]
        for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
        p = make_subplots(rows=nrows, cols=1, 
        shared_xaxes=true, 
        shared_yaxes=false,
        vertical_spacing=0.05,
        )
        for i in 1:nrows;
                add_trace!(p, 
                scatter(x=df.date, y=df[:,i],
                name=st[i]),   row=i,     col=1);
        end
        relayout!(p,height=600*1.5,width=900*1.5)
        return(p)
    end

    function dfplot(df::AbstractString)
        df=waread(df)
        nrows=size(df)[2]-1
        st=[]
        for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
        p = make_subplots(rows=nrows, cols=1, 
        shared_xaxes=true, 
        shared_yaxes=false,
        vertical_spacing=0.05,
        )
        for i in 1:nrows;
                add_trace!(p, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=st[i]),   row=i,     col=1);
        end
        PlotlyJS.relayout!(p,height=600,width=900)
        display(p)
    end

    plotdf=dfplot

    function pline(path::AbstractString)
        df = wa.waread(path)

        nrows=size(df)[2]-1
        st=[]
        for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
        p = make_subplots(rows=nrows, cols=1, 
        shared_xaxes=true, 
        shared_yaxes=false,
        vertical_spacing=0.05,
        #subplot_titles= st;
        )
        for i in 1:nrows;
                add_trace!(p, 
                PlotlyJS.scatter(x=df.date, y=df[:,i],
                name=st[i]),   row=i,     col=1);
        end
        #relayout!(p,height=600*2,width=900*2,title_text="Series of "*basename(path))
        PlotlyJS.relayout!(p,height=600*1.5,width=900*1.5,title_text="Series of "*basename(path))
        display(p)
    end

    function dfplotjs(df::DataFrame;logy::Bool,fact::Float64)
        nrows=size(df)[2]-1 
        #length(names(df))-1
        o = DataFrames.metadata(df)|>collect
        ti = basename(o[1][2])
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
        #rows=2, cols=2
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact = isnothing(fact) ? 1 : fact; #nice
        logy = isnothing(logy)==true ? logy==false : logy==true;
        if logy == true
            PlotlyJS.relayout!(fig,yaxis_type="log",
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        #elseif isnothing(log) 
        else
            PlotlyJS.relayout!(fig,
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        end
        display(fig)
    end

    function dfplotjs(df::AbstractString;logy::Bool,fact::Float64)
        df=waread(df)
        nrows=size(df)[2]-1
        #length(names(df))-1
        o = DataFrames.metadata(df)|>collect
        ti = basename(o[1][2])
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact = isnothing(fact) ? 1 : fact; #nice
        logy = isnothing(logy)==true ? logy==false : logy==true;
        if logy == true
            PlotlyJS.relayout!(fig,yaxis_type="log",
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        else
            PlotlyJS.relayout!(fig,
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        end
        display(fig)
    end

    function dfplotjs(filepath::AbstractString)
        dfplotjs(filepath;logy=false,fact=1.0)
    end

    function dflogjs(filepath::AbstractString)
        dfplotjs(filepath;logy=true,fact=1.0)
    end

    function lplotjs(x::Regex)
        df=waread(x)
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
        ti = raw"" 
        end
        #o = collect(DataFrames.metadata(df))[1][2] |>basename
        ln = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
        #@df df Plots.plot(:date,cols(ln),yaxis=:log,title=ti)
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        
        #for i in 1:size(df)[2]-1
        for i in ln
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(
                x=df.date, 
                #y=select(df,i),
                y=df[:,i],
                name=i))
                #name=names(df)[i]));
        end
        PlotlyJS.relayout!(fig,
                template="seaborn",        
                yaxis_type="log",
                #height=600*fact,width=900*fact,
                #title_text=ti*" (log)",
                title_text=ti,
                xaxis_rangeslider_visible=true,
                updatemenus=[
                    Dict(
                        "type" => "buttons",
                        "direction" => "left",
                        "buttons" => [
                            Dict(
                                "args" => [Dict("yaxis.type" => "linear")],
                                "label" => "Linear Scale",
                                "method" => "relayout"
                            ),
                            Dict(
                                "args" => [Dict("yaxis.type" => "log")],
                                "label" => "Log Scale",
                                "method" => "relayout"
                            )
                        ],
                        "pad" => Dict("r" => 1, "t" => 10),
                        "showactive" => true,
                        "x" => 0.11,
                        "xanchor" => "left",
                        "y" => 1.1,
                        "yanchor" => "auto"
                    )]
                
                
                )
        return(fig)
    end

    function lplotjs(x::DataFrame)
        df = x 
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
        ti = raw"" 
        end
        ln = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );

        for i in ln
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(
                x=df.date, 
                #y=select(df,i),
                y=df[:,i],
                name=i))
                #name=names(df)[i]));
        end
        PlotlyJS.relayout!(fig,
                template="seaborn",        
                yaxis_type="log",
                #height=600*fact,width=900*fact,
                #title_text=ti*" (log)",
                title_text=ti,
                xaxis_rangeslider_visible=true,
                updatemenus=[
                    Dict(
                        "type" => "buttons",
                        "direction" => "left",
                        "buttons" => [
                            Dict(
                                "args" => [Dict("yaxis.type" => "linear")],
                                "label" => "Linear Scale",
                                "method" => "relayout"
                            ),
                            Dict(
                                "args" => [Dict("yaxis.type" => "log")],
                                "label" => "Log Scale",
                                "method" => "relayout"
                            )
                        ],
                        "pad" => Dict("r" => 1, "t" => 10),
                        "showactive" => true,
                        "x" => 0.11,
                        "xanchor" => "left",
                        "y" => 1.1,
                        "yanchor" => "auto"
                    )]
                
                
                )
        return(fig)
    end

    function recursive_glob_prfx(rootdir=".", prefix="")
        results = []
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
        return results
    end

    #fwin=raw"C:\Users\Public\Documents\Python_Scripts\julia\func-win.jl"
    # fwin=raw"C:\Users\chs72fw\.julia\dev\WaSiM\src\WaSiM.jl"
    # include(fwin)

    function waread2(x::String)
        """
        Read the text file, preserve line 1 as header column
        Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
        This avoids reading the entire file into memory at once, 
        which can be more memory-efficient for large datasets.
        """
        ms = ["-9999", "lin", "log", "--"]
        df = CSV.File(x; delim="\t", header=1, normalizenames=true, missingstring=ms, types=Float64) |> DataFrame
        dropmissing!(df,1)
        dt2 = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
        select!(df, Not(1:4))
        df.date = dt2
        metadata!(df, "filename", x, style=:note)
        return df
    end

    function waread(x::String)
        """
        Fastest Reader. is also dfr.
        Read the text file, preserve line 1 as header column
        """
        ms = ["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame; delim="\t", header=1, missingstring=ms, normalizenames=true, types=Float64)
        df = dropmissing(df, 1)
        dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
        df.date = dt2
        df = select(df, Not(1:4))
        DataFrames.metadata!(df, "filename", x, style=:note)
        for x in names(df)
            if startswith(x,"_")
                newname=replace(x,"_"=>"C", count=1)
                rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    function waread(x::Regex)
        """
        Read the text file, preserve line 1 as header column
        """
        x = dfonly(x)|>first
        ms = ["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame; delim="\t", header=1, missingstring=ms, normalizenames=true, types=Float64)
        df = dropmissing(df, 1)
        dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
        df.date = dt2
        df = select(df, Not(1:4))
        metadata!(df, "filename", x, style=:note)
        #renamer
        for x in names(df)
            if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    function dfonly(x1::AbstractString)
        v = filter(file -> occursin(Regex(x1,"i"),file), readdir());
        z = v[broadcast(x->!endswith(x,r"nc|png|svg|jpg"),v)];
        return(z)
    end
    
    function dfonly(x1::Regex)
        z = filter(file -> occursin(x1,file), 
        readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
        return(z)
    end

    function kge2(simulated::Vector{Float64}, observed::Vector{Float64})
        r = cor(simulated, observed)
        α = std(simulated) / std(observed)
        β = mean(simulated) / mean(observed)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    function kge2(df::DataFrame)
        observed, simulated = df[:,6],df[:,5]
        r = cor(observed, simulated)
        α = std(simulated) / std(observed)
        β = mean(simulated) / mean(observed)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    function kge_read(path::AbstractString, ext::AbstractString)
        files = readdir(path)
        for file in files
            file_path = joinpath(path, file)
            if isfile(file_path) && endswith(file, ext) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|tif|jpg|png|svg",file))
                dd = CSV.read(file_path,DataFrame,missingstring="-9999",delim="\t")
                # observed  = dd[:,5]
                # simulated = dd[:,6]
                # kge_value = kge2(observed, simulated)
                kge_value = kge2(dd)
                println(replace("KGE value is $kge_value on $file_path", "\\"  => "/"))
            elseif isdir(file_path)
                dfs_in_subdir = kge_read(file_path, ext)
            end
        end
    end

    function nse(predictions::Vector{Float64}, targets::Vector{Float64})
        return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
    end

    function nse(df::DataFrame)
        simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
        return (1 - (sum((simulated .- observed).^2) / sum((observed .- mean(observed)).^2)))
    end

    function kge(df::DataFrame)
        simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
        r = cor(observed, simulated)
        α = std(simulated) / std(observed)
        β = mean(simulated) / mean(observed)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end


    #function dfpl(df::DataFrame;logy::Bool,fact::Float64)
    function dfljs(df::DataFrame;logy=true,fact=.66)
        nrows=size(df)[2]-1 
        o = DataFrames.metadata(df)|>collect
        ti = basename(o[1][2])
        fig = PlotlyJS.make_subplots(
            shared_xaxes=true, 
            shared_yaxes=true    
            );
        for i in 1:nrows;
            PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:,i],
            name=names(df)[i]));
        end
        fact = isnothing(fact) ? 1 : fact; #nice
        logy = isnothing(logy)==true ? logy==false : logy==true;
        if logy == true
            PlotlyJS.relayout!(fig,
            template = "seaborn",
            yaxis_type="log",
            height=600*fact,width=900*fact,
            title_text=o[1][2]*" (log)")
        else
            PlotlyJS.relayout!(fig,
            height=600*fact,width=900*fact,
            title_text="Series of "*ti)
        end
        return fig
    end

      


end #module