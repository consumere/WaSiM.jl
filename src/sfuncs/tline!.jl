# 
function tline!(df::DataFrame; date_col=:date, lab=false)
        # Get the date column and column names for trendlines
        date_data = df[!, date_col]
        trendline_cols = setdiff(names(df), [string.(date_col)])
    
        # get the last viewed plot # yes.wrks.
        p = current()
        
        for y_col in trendline_cols
            y_data = df[!, y_col]
    
            # Perform linear regression to get the slope and intercept
            X = hcat(ones(length(date_data)), 1:length(date_data))
            y = y_data
            β = X \ y  # Linear regression
    
            # Extract the intercept and slope
            intercept, slope = β[1], β[2]
    
            # Generate the trendline 