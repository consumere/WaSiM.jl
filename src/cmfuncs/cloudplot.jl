# 
function cloudplot(df::DataFrame)

        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
        end

        a,b = ctov(df)
        # colors = Makie.wong_colors()
        # Int.(size(unique(a)))[1]
        colors = Makie.wong_colors()
        # Get unique values
        unique_values = unique(a)
        # Map unique values to colors
        value_to_color = Dict(unique_values[i] => colors[i % length(colors) + 1] for i in 1:length(unique_values))
        # Map values in a to colors
        selected_colors = [value_to_color[value] for value in a]
   
        #category_labels = names(df[!,Not(:date)]),
        
        rainclouds(a,b;
            xlabel = "", #"Basins",
            #ylabel = " ", 
            
            orientation = :horizontal,
            title = ti,
            plot_boxplots = true, 
            cloud_width=2.5, 
            side = :right, violin_limits = extrema,
            #xticklabelrotation = π / 4,
            #xticklabelalign = (:right, :center),
            #clouds=hist,
            color = selected_colors)
            # ax3.xticklabelrotation = π / 4
            # ax3.xticklabelalign = (:right, :center)

        #cbs = names(df[!,Not(:date)])
        #CairoMakie.yticks!(;ytickrange=1:length(cbs),yticklabels=cbs)

    end

    """
    cloudplot2(df::DataFrame)
    with less options
    """
    