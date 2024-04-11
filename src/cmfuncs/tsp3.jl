# 
function tsp3(df::DataFrame)
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
    
        fig = Figure( #resolution=(800, 800), 
        fonts=(;regular = "consolas")            
        )
                     
        
        for (i, col) in enumerate(cols)
            row, cl = fldmod1(i, 2)
            ax = Axis( fig[row, cl],
                       title = replace(tit, r"_|so" => " "),
                       xlabel = "Date",
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
    barplots yearsum timeseries
    baryrsum
    """
    