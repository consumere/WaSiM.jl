# 
function moisture_plot_with_confidence(df, timestep=month; kwargs...)
        # Extracting relevant columns
        dk = select(df, [:date, :tot_average])
        if timestep != :day
            dk[!, :agg_date] .= timestep.(dk[!,:date])
        else
            dk[!, :agg_date] .= dk[!,:date]
        end    

        # Grouping by the specified timestep and calculating the minimum and maximum values of tot_average
        su = DataFrames.combine(groupby(dk, :agg_date), 
            :tot_average .=> (minimum, maximum) .=> 
            [:min_moisture, :max_moisture])

            lab = uppercase(string(timestep))
        # Creating the plot with confidence bands
        plot(su.agg_date, [su.min_moisture su.max_moisture], ribbon=true, fillalpha=0.2,
            xlabel=string(lab), 
            ylabel="Total Average Moisture",
            title="Moisture Plot with Confidence Bands",
            label=["Min Moisture" "Max Moisture"]; kwargs...)
    end

    