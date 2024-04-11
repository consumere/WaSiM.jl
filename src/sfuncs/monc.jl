# 
function monc(df::DataFrame; confidence_level=0.95)
        # Make a copy of the input DataFrame
        dmean = copy(df)
        
        # Extract the columns of interest
        columns = filter(x -> !occursin("date", x), names(dmean))
        
        # Add a 'month' column based on the 'date' column
        dmean[!, :month] = month.(dmean[!, :date])
        
        # Drop the 'date' column
        select!(dmean, Not(:date))

        dmean = DataFrames.combine(DataFrames.groupby(dmean, :month)) do group
            result = DataFrame(month = group.month)
            for col in columns
                values = group[!, col]
                mean_value = mean(values)
                
                z_score = quantile(Normal(), 0.5 + confidence_level / 2) # Use 0.975 for 95% confidence level
                margin_of_error = z_score * std(values) / sqrt(length(values))
                lower_bound = mean_value - margin_of_error
                upper_bound = mean_value + margin_of_error
                
                result[!, col * "_mean"] .= mean_value
                result[!, col * "_lower"] .= lower_bound
                result[!, col * "_upper"] .= upper_bound
            end
            return result
        end
    
        return dmean
    end

    """
    monthly mean plot with ribbon
    x::Union{String,Regex,DataFrame};col::Any=1,confidence_level::Number=0.95)
    dfrib(df;confidence_level=.999,col=:tot_average) 
    """
    