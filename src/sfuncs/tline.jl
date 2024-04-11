# 
function tline(df::DataFrame, date_col::Symbol=:date;lab=true)
        # Get the date column and column names for trendlines
        date_data = df[!, date_col]
        trendline_cols = setdiff(names(df), [string.(date_col)])
    
        p = plot()
    
        # Initialize an empty string to store the trendline equations
        trendline_equations = ""
    
        for y_col in trendline_cols
            y_data = df[!, y_col]
    
            # Perform linear regression to get the slope and intercept
            X = hcat(ones(length(date_data)), 1:length(date_data))
            y = y_data
            β = X \ y  # Linear regression
    
            # Extract the intercept and slope
            intercept, slope = β[1], β[2]
    
            # Generate the trendline using the linear equation
            trendline(x) = intercept + slope * x
            # Add the trendline equation to the trendline_equations string
            #trendline_equations *= "$y_col: y = $(round(slope, digits=2))x + $(round(intercept, digits=2))\n"
            if lab
                trendline_equations = "$y_col\ny = $(round(slope; sigdigits=2))x + $(round(intercept, sigdigits=1))"
            else
                trendline_equations = false
            end
    
            # plot
            plot!(p ,date_data, 
                trendline.(1:length(date_data)), 
                label=trendline_equations,
                linewidth=2, linestyle=:dash)
    
            
        end
        return p
    end


    """
    adds trendlines to existing plot
    (df::DataFrame; date_col=:date, lab=false)
    """
    