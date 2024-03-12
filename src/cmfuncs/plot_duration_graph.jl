# 
function plot_duration_graph(log_file::Tuple{Vector{Any}, Vector{Any}})
        timestamps, durations = log_file
        # Sort data by timestamps
        sorted_indices = sortperm(timestamps)
        sorted_timestamps = timestamps[sorted_indices]
        sorted_durations = durations[sorted_indices]

        sorted_timestamps = map(x->
            replace(basename(x), "xml" => splitdir(dirname(x))[2]), 
                sorted_timestamps)
    
        # Create bar plot
        p = plot(
            sorted_timestamps,
            sorted_durations,
            seriestype = :bar,
            #xlabel = "XML Files",
            xlabel = "",
            ylabel = "Duration (minutes)",
            title = "Duration of WaSiM runs",
            legend = false,
            xrotation = -35,  # Rotate x-axis labels for better readability
            grid = true,
            right_margin = 10mm,
            left_margin = 5mm,
            top_margin = 5mm,
            bottom_margin = 15mm, 
            bar_width = 0.6);
        
    
        # Add annotations for duration time in minutes
        for (x, y) in zip(sorted_timestamps, sorted_durations)
            annotate!(x, y, Plots.text("$y min",7,:black,:bottom,
            rotation=0            
            ))
        end
    
        return p
    end

    #