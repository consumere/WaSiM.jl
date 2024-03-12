# 
function vibx(df::DataFrame)
        str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
        ln = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
        @df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
        @df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
        @df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
    end

    """
    vioplot wasim timeseries 
    """
    