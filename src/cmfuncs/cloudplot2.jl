# 
function cloudplot2(df::DataFrame)
    
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
        end
    
        a,b = cmk.ctov(df)
        colors = Makie.wong_colors()
        unique_values = unique(a)
        value_to_color = Dict(unique_values[i] => colors[i % length(colors) + 1] for i in 1:length(unique_values))
        selected_colors = [value_to_color[value] for value in a]
    
        rainclouds(a,b;
            xlabel = "",
            ylabel = " ", 
            title = ti,
            gap=0.0002,
            #cloud_width = maximum(b),
            #clouds = violin,
            #violin_limits = extrema(b).*.8,
            xticklabelalign = (:right, :center),
            color = selected_colors)
            #clouds=hist,
            #xticklabelrotation = π / 4,
    end
    
    """
    plots timeseries
    """
    