# 
function dfp(x::Union{Regex,String,DataFrame};leg::Bool=true)
        if isa(x,DataFrame)
            df = (x)
        else
            df = waread(x)
        end

        dt = df.date
        fig = Figure()
                # fonts=(;regular = "consolas")            
                # #size=(600, 400), 
                # )
        
        tempo = string.(dt)
        lentime = size(df, 1)
        #slice_dates = range(1, lentime, step=lentime ÷ 8)
        slice_dates = range(1, lentime, step=lentime ÷ 10)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
        
        ax3 = Axis(fig[1, 1]) 

        cols = names(df)
        filter!(x -> !occursin(r"date|year", x), cols)

        colors = Makie.wong_colors()
        # Map unique values to colors
        value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # Map values in a to colors
        selected_colors = [value_to_color[value] for value in cols]

        for (i, col) in enumerate(cols)
            vals = Vector(df[:, col])
            lines!(ax3, 1:lentime, vals; color = selected_colors[i], linewidth = 0.85, label=col)
        end
        
        if leg
            lt = replace(tit, r"_|so" => " ")
            legend = Legend(fig, ax3, lt, 
                framevisible = true,
                nbanks = 2)

            fig[1, 2] = legend
        end

        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        
        return fig
    end

    tsdf = dfp

    """
    density plot of dataframe monthly values
    """
    