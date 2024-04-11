# 
function monc_f(df::DataFrame; confidence_level=0.95)
        # Make a copy of the input DataFrame
        dmean = copy(df)
        

        
        # Add a 'month' column based on the 'date' column
        dmean[!, :month] = month.(dmean[!, :date])
        
        # Drop the 'date' column
        select!(dmean, Not(:date))
        # Extract the columns of interest
        columns = filter(x -> !occursin(r"date|month", x), names(dmean))
        # before aggregation.
        for col in columns
            values = dmean[!, col]
            mean_value = mean(values)
            
            z_score = quantile(Normal(), 0.5 + confidence_level / 2) # Use 0.975 for 95% confidence level
            margin_of_error = z_score * std(values) / sqrt(length(values))
            lower_bound = mean_value - margin_of_error
            upper_bound = mean_value + margin_of_error
            
            dmean[!, col * "_mean"] .= mean_value
            dmean[!, col * "_lower"] .= lower_bound
            dmean[!, col * "_upper"] .= upper_bound
        end

        # get new names
        y = filter(x -> !occursin(r"date|month", x), names(dmean))
        dmean = DataFrames.combine(
            groupby(dmean, :month), 
            y .=> mean .=> y);

    
        return dmean
    end

    """
    running momean and add confidence intervals
    df::DataFrame; confidence_level=0.95
    """
    