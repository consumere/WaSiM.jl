# 
function qpl(x::AbstractString)
        df = readdf(x)
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
        end
        s = names(df)[1:2]
        t2 = string.(ti,"\n",s[1],"|",s[2],ti)
        StatsPlots.plot( 
        qqplot(df[!,1],df[!,2], qqline = :fit), 
        qqplot(Cauchy,df[!,2]), 
        qqnorm(df[!,2], qqline = :R),
        title = t2)
    end

    qqp=qpl

    