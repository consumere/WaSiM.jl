# 
function mbx(df::DataFrame;fun=mean)
        df.Month = month.(df.date)
        str = [ @sprintf("%02i", x) for x in (df.Month) ];
        
        month_abbr = ["Jan", "Feb", "Mar", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]

        ln = Symbol.(filter(x->!occursin(r"date|year|month"i,x),names(df)))
        # values = means[:, ncol(means)]
        # colors = colorfunction(values)

        # ln = Symbol.(filter(x->!occursin(r"date|year|month|mean"i,x),names(df)))

        p = @df df StatsPlots.boxplot(
            str , cols(ln),
            fillalpha=0.75, 
            linewidth=0.25,
            # seriescolor=:heat,
            # color = colors,
            notch = true,
            whisker_width = :match,
            legend=false)
        xticks!(0.5:11.5 , month_abbr)
        means = DataFrames.combine(groupby(df,:Month), ln[1] => fun)
        
        for i in eachrow(means)
            m = i[2]
            annotate!(i.Month - 0.5, m, #+ 1 
            text(round(m; digits=2), 6, :center, :top))
        end
        return p
    end

    """
    correlation plots on dataframe
    size=(1200,800)
    has to be sim,obs
    """
    