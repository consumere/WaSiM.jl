# 
function rmm(x::DataFrame;fun=mean)
        # Ensure the date column is of type Date
        df = copy(x)
        # Add a year and month column
        df.year = year.(df.date)
        df.month = month.(df.date)

        # Group by year and month and calculate the mean
        aggdat = DataFrames.combine(groupby(
            df, [:year, :month]), 
            names(df, Not([:date, :year, :month])) .=> fun;
            renamecols=false)

        # Recombine year and month to a date column
        aggdat.date = Date.(aggdat.year, 
            aggdat.month)

        # Remove the year and month columns
        select!(aggdat, Not([:year, :month]))

        return aggdat
    end

    