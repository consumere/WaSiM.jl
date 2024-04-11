# 
function aplot(df::Regex)
        df = readf(df)
        df = copy(df)
        df[!,:year]=year.(df[!,:date]) ;
        s = Symbol.(filter(x->!occursin("date",x),names(df)))
        o = DataFrames.metadata(df)|>collect
        ti = "AndrewsPlot of "*basename(o[1][2])
        @df df andrewsplot(:year, cols(s), legend = :topleft,title=ti)
    end
    
    