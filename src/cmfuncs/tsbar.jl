# 
function tsbar(df::DataFrame)
        df = yrsum(df)
        #dt = df.date
        dt = df.year
        fig = Figure( #resolution=(600, 400), 
                fonts=(;regular = "consolas")            
                )
    
        
        tempo = string.(dt)
        lentime = size(df, 1)
        #slice_dates = range(1, lentime, step=lentime ÷ 8)
        slice_dates = range(1, lentime, step=1)
        tit = try
            DataFrames.metadata(df) |> only |> last |> basename
            catch
                @warn "no metadata in df"
                raw""
            end
                
        
        ax3 = Axis(fig[1, 1], 
                   title = replace(tit, r"_|so" => " "),
                   xlabel = raw" ",
                   ylabel = first(names(df)))
  
        cols = names(df)
        filter!(x -> !occursin(r"date|year", x), cols)
        
        colors = Makie.wong_colors()
        # Map unique values to colors
        value_to_color = Dict(cols[i] => colors[i % length(colors) + 1] for i in 1:length(cols))
        # Map values in a to colors
        selected_colors = [value_to_color[value] for value in cols]
        
        # dfm = DataFrames.stack(df)
        # dfm.variable = string.(dfm.variable)
        # fig = Figure()
        # ax = Axis(fig[1,1])
        # barplot!(ax,
        #     dfm.year,
        #     dfm.value,
        #     color=dfm.year,
        #     stack=dfm.year )
      


        for (i, col) in enumerate(cols)
            #vals = select(df, col) |> Matrix |> vec
            vals = Vector(df[:, col]) #more efficient
            #lines!(ax3, 1:lentime, vals; color = selected_colors[i], linewidth = 0.85)
            # barplot!(ax3, 1:lentime, vals; 
            # stack = df.year,
            # color = selected_colors[i], alpha = 0.85)
            barplot!(ax3, 1:lentime, vals,
            stack = df.year,
            color = selected_colors[i])
        end
        
        ax3.xticks = (slice_dates, tempo[slice_dates])
        ax3.xticklabelrotation = π / 4
        ax3.xticklabelalign = (:right, :center)
        
        return fig
    end

    