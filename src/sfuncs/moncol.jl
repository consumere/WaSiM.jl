# 
function moncol(x::DataFrame)
        dmean = copy(x)
        s = map(Symbol, filter(x -> !occursin("date", x), names(dmean)))
        dmean[!, :month] = month.(dmean[!,:date])
        select!(dmean, Not(:date))
    
        # Calculate monthly mean and confidence interval using combine and groupby
        dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
            s .=> (dmean -> (mean(dmean), 1.96 * std(dmean) / sqrt(length(dmean)))) 
            .=> s)
    
        # Create separate columns for mean, lower bound, and upper bound
        for col in s
            dmean[!, string(col)*"_mean"]  .= [first(only(x)) for x in eachrow(dmean[!, col])]
            dmean[!, string(col)*"_lower"] .= [first(only(x)) - last(only(x)) for x in eachrow(dmean[!, col])]
            dmean[!, string(col)*"_upper"] .= [first(only(x)) + last(only(x)) for x in eachrow(dmean[!, col])]
        end
    
        # # Drop the original columns representing the mean and standard deviation
        # select!(dmean, Not.(:month, s...))
    
        return dmean
    end

    """
    aggregate and add confidence intervals
    df::DataFrame; confidence_level=0.95
    """
    