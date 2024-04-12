function flow_duration_curve(flow; ti = nothing, normal = false, 
    gust = true)
    g = zeros(20, 7)
    g[1, :] = [975.7, 577.26, 20.49, 3.7, 1.73, 1, 0.38]
    g[2, :] = [904.17, 534.08, 22.69, 4.42, 2.13, 1.26, 0.51]
    g[3, :] = [838.77, 511.37, 25.1, 5.27, 2.62, 1.58, 0.67]
    g[4, :] = [776.04, 480.48, 27.86, 6.33, 3.25, 2, 0.88]
    g[5, :] = [719.91, 452.42, 30.82, 7.54, 3.99, 2.51, 1.16]
    g[6, :] = [667.48, 425.82, 34.11, 9, 4.92, 3.16, 1.53]
    g[7, :] = [618.22, 400.44, 37.81, 10.77, 6.07, 3.98, 2.02]
    g[8, :] = [572.53, 376.64, 41.82, 12.86, 7.47, 5.01, 2.65]
    g[9, :] = [520, 350.65, 45.1, 15.2, 9.16, 6.3, 3.46]
    g[10, :] = [472.29, 326.46, 48.64, 17.98, 11.22, 7.94, 4.52]
    g[11, :] = [428.96, 303.93, 52.46, 21.25, 13.75, 10, 5.89]
    g[12, :] = [389.6, 282.96, 56.57, 25.13, 16.86, 12.57, 7.69]
    g[13, :] = [353.86, 263.44, 61.01, 29.71, 20.66, 15.83, 10.03]
    g[14, :] = [321.39, 245.26, 65.79, 35.12, 25.32, 19.93, 13.08]
    g[15, :] = [291.65, 228.19, 71, 41.58, 31.09, 25.13, 17.11]
    g[16, :] = [264.89, 212.45, 76.57, 49.16, 38.1, 31.64, 22.32]
    g[17, :] = [240.09, 197.49, 82.6, 58.08, 46.67, 39.81, 29.13]
    g[18, :] = [206.89, 176.99, 89.91, 67.82, 56.95, 50.13, 39]
    g[19, :] = [178.28, 158.62, 97.86, 79.21, 69.5, 63.12, 52.22]
    g[20, :] = [153.69, 142.2, 106.49, 92.46, 84.77, 79.43, 69.85]
    
    p = [0.02, 0.05, 0.5, 0.8, 0.9, 0.95, 0.99]
    rank = sortperm(flow, rev=true)
    #rank .= maximum(rank) .- rank

    exceedtime = 1.0 * (rank .- 1) / length(flow)
    q = sort(100 * flow / mean(flow), rev=true)
    
    yl = "Percent of mean discharge"
    xl = "Exceedance probability (%)"
    ylims = (minimum(g), maximum(q))
    
    if normal
        exceed_z = quantile(Normal(), exceedtime)
        p_z = quantile(Normal(), p)
        xlims = (-3, 3)
        
        plot(exceed_z, q, 
        line=:solid, linewidth=1.2, 
        color=:blue, 
        #plot(1:length(q), q, line=:solid, linewidth=2, color=:blue, 
        yscale=:log2, 
        legend=:outertopright,
        xaxis=false, ylim=ylims, xlim=xlims, xlabel=xl, ylabel=yl)
        title!(ti, loc=3, line=1, c=:black, alpha=0.7)
        
        probs = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999]
        z_vals = quantile(Normal(), probs)
        xticks!(z_vals, labels=probs)
        
        if gust
            for k in 1:20
                plot!(p_z, g[k, :], line=:solid, color=:gray50)
                annotate!(p_z[7], g[k, 7], text(k - 1, color=:gray50, halign=:right, valign=:center, pointsize=8))
            end
            #plot!(exceed_z, q, line=:solid, linewidth=2, 
            # legend=:outertopright,
            #     color=:blue)
            #legend!("Flow Duration Curve with Gustard's Type Curves", color=:gray50, framestyle=:none)
            annotate!(0.2, 2.5,text("Flow Duration Curve with Gustard's Type Curves"), color=:gray50, halign=:left, valign=:center, pointsize=12)
            annotate!(2.5, g[20, 7], text("permeable", color=:gray50, halign=:left, valign=:center, pointsize=8))
            annotate!(2.5, g[1, 7], text("impermeable", color=:gray50, halign=:left, valign=:center, pointsize=8))
        end
        vline!([0.5], linestyle=:dash, linecolor=:red)
        hline!([100], linestyle=:dash, linecolor=:red)
    else
        xlims = (0, 1)
        plot(exceedtime, q, 
            seriestype=:line, line=:solid, linewidth=2,
            #seriestype=:scatter,
            legend=:outertopright, color=:blue, 
        yscale=:log2, ylim=ylims, xlim=xlims, xlabel=xl, ylabel=yl)
        title!(ti, loc=3, line=1, c=:black, alpha=0.7)
        
        if gust
            for k in 1:19
                plot!(p, g[k, :], line=:solid, color=:gray50)
                annotate!(p[7], g[k, 7], text(k, color=:gray50, halign=:right, valign=:center, pointsize=8))
            end
            #plot!(exceedtime, q, line=:solid, linewidth=2, color=:blue,legend=false)
                #legend="Flow Duration Curve with Gustard's Type Curves")
            #legend!("Flow Duration Curve with Gustard's Type Curves", color=:gray50, framestyle=:none)
            annotate!(p[7], g[19, 7], text("permeable", color=:gray50, halign=:left, valign=:center, pointsize=8))
            annotate!(p[7], g[1, 7], text("impermeable", color=:gray50, halign=:left, valign=:center, pointsize=8))
        end
        vline!([0.5], linestyle=:dash, linecolor=:red)
        hline!([100], linestyle=:dash, linecolor=:red)
    end
end


fl = dfr(r"Wolf")
flow_duration_curve(fl.Wolfsmuenster;ti= "Wolfsmuenster",gust=false)
flow_duration_curve(fl.Wolfsmuenster;ti= "Wolfsmuenster",normal=true)
flow_duration_curve(fl.Wolfsmuenster;ti= "Wolfsmuenster")
flow_duration_curve(fl.Wolfsmuenster;ti= "Wolfsmuenster",normal=false,gust=false)

flow = fl.Wolfsmuenster
rank = sortperm(flow, rev=true)
rank .= maximum(rank) .- rank
#rank = permutedims(rank)|>vec
#rank|>scatter
scatter(rank,0:length(rank))


rank = sortperm(flow, rev=true)
rank .= maximum(rank) .- rank
#rank greater than 0 
rank = rank[rank.>0]

exceedtime = 1.0 * (rank .- 1) / length(rank)
q = sort(100 * flow / mean(flow), rev=false)
exceed_z = quantile(Normal(), exceedtime)
plot(1:length(q), q, line=:solid, linewidth=2, color=:blue, 
        yscale=:identity)

plot(reverse(q),1:length(q), line=:solid, linewidth=2, color=:blue, 
        yscale=:identity)
        

scatter(exceedtime, q)