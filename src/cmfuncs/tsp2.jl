# 
function tsp2(df::DataFrame)
        dt = df.date
        tempo = string.(dt)
        lentime = size(df, 1)
        slice_dates = range(1, lentime, step=lentime ÷ 8)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        cols = filter(x -> !occursin(r"date|year", x), names(df))
    
        fig = Figure( #resolution=(800, 400), 
        fonts=(;regular = "consolas")            
        )
                     
        
        for (i, col) in enumerate(cols)
            ax = Axis(fig[1, i], 
                       title = replace(tit, r"_|so" => " "),
                       xlabel = raw"", #"Date",
                       ylabel = col)
            
            #vals = select(df, col) |> Matrix |> vec
            vals = Vector(df[:, col])
            
            lines!(ax, 1:lentime, vals; linewidth = 0.85)
            
            ax.xticks = (slice_dates, tempo[slice_dates])
            ax.xticklabelrotation = π / 4
            ax.xticklabelalign = (:right, :center)
        end
        
        return fig
    end

    """
    plots timeseries for each column separately in one figure
    """
    