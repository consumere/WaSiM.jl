# 
function tbx(df::DataFrame)
        ti = try
            DataFrames.metadata(df) |> only |> last |> basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
    

        month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
        if any(x->occursin("year",x),names(df))
            #df = df[!,Not("year")]
            s = Symbol.(filter(x->!occursin(r"date|year|month|day",x),names(df)))
            #str = unique(df.year)
            @df df StatsPlots.violin(cols(s), linewidth=0.1)
            #xticks!(0.5:11.5, month_abbr)
            title!(ti)
        else
            str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            fnt = Plots.font(
                family="sans-serif",
                pointsize=8,
                valign=:bottom,
                rotation=0.0,
                color=:black
            )
    
            dx = monmean(df)
            y_offset = dx[!,2]
            anns = map(x->string.(round(x,digits=2)), y_offset)
            xans = map(x->Plots.text(x, fnt), anns)
            #col_annotations = (Vector(dx[!, 1]) .- 0.5, y_offset, xans) # x y val
            #col_annotations = (Vector(dx[!, 1]) .- 0.5, y_offset, xans) # x y val
            col_annotations = (Vector(dx[!, 1]) .+ 0.5, y_offset, xans) # x y val
    
            # Create a colormap from the `monmean` values
            #k=colorfunction(dx[!,2])
            #k = cf2(dx[!,2])
            #mat = reshape(k, (1,nrow(dx)))
            #colors = colormap.(dx[!,2])
    
            @df df StatsPlots.boxplot(
                str,
                cols(s),
                group = str, # group by month, needed for colormap
                linewidth=0.1,
                outliers = false,
                whisker_width = :match,
                #colors=map(x -> (x > 3 ? :green : :red), y_offset),
                #fillcolor = mat|>reverse,
                #fillcolor = Colors.gray.(1:12),
                #fillcolor = :grey,
                colors = Plots.colormap("RdBu",size(dx, 2),mid=0.2),
                annotations = col_annotations,
                legend=false,
            )
    
            xticks!(0.5:1.4:16, month_abbr)
            #xticks!(0.5:1.5:12.5, month_abbr)
            title!(ti)
        end
    end

    """
    moisture_plot_with_confidence(df, timestep=month; kwargs...)
    """
    