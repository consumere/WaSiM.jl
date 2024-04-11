# 
function tspblack(df::DataFrame)
        dt = df.date
        fig = Figure( #resolution=(600, 400), 
            fonts=(;regular = "consolas")            
            )
        
        tempo = string.(dt)
        #lentime = nrow(df)
        lentime = size(df,1)
        slice_dates = range(1, lentime, step=lentime ÷ 8)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        
        ax3 = Axis(fig[1, 1], 
            title = replace(tit,r"_|so"=>" "),
                xlabel="Date", 
                #yscale = :log10,
                ylabel=first(names(df)))
    
        #cols = propertynames(df)  
        cols = names(df)  
        filter!(x->!occursin(r"date|year",x),cols)

        # colors = Makie.wong_colors()
        # # Map unique values to colors
        # value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # # Map values in a to colors
        # selected_colors = [value_to_color[value] for value in cols]
        # colors = Makie.wong_colors()
        #v = 1:length(cols)
        #for v, col in enumerate(cols)
        for col in cols
            #vals = select(df,col)|>Matrix|>vec   
            vals = Vector(df[:, col])
            lines!(ax3, 1:lentime, vals; 
            color=:black,         #     color=colors[v], 
             linewidth=0.85)
        end
        
        #vals2 = select(df,2)|>Matrix|>vec
        #line1 = lines!(ax3, 1:lentime, vals; color=:black, linewidth=0.85)
        #line2 = lines!(ax3, 1:lentime, vals2; color=:red, linewidth=0.85)
        #
        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        return fig
    end

    """
    plots timeseries
    """
    