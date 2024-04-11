# 
function luscatter(dfs::Vector=dfs,x::String="rs_evaporation")
        month_abbr = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]
        p1 = plot(title=x);
        for i in 1:size(dfs,1)
            value_str = findindf(dfs[i], x).Value
            #value_str = replace.(value_str, r" " => "")
            value_str = string.(value_str)
            #value_str = strip.(value_str)
            vs = split.(value_str)
            z0 = [parse.(Float64, z) for z in vs]
            cls = split(dfs[i][1, 1], r"\W")[2] 
            scatter!(month_abbr, first(z0), label=cls,
                yaxis=:identity) #log10
            #push!(plots, plot(month_abbr, only(z0), label=cls))
        end
        #title!("LAI in [m2/m2]")
        Plots.show(p1)
        return p1    
    end

    """
    heatmap of landusedata 
    dfs::Vector=dfs,x::String="LAI"
    """
    