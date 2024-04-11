# 
function ovio(df::DataFrame)
        str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
        end    
        if (any(x->occursin("year",x),names(df)))
            #s = Symbol.(filter(x->!occursin("year",x),names(df)))
            df = df[!,Not("year")]
            s = Symbol.(filter(x->!occursin("date",x),names(df)))
            @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
            #@df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
            title!(ti)
        else    
        s = Symbol.(filter(x->!occursin("date",x),names(df)))
        @df df StatsPlots.violin(str,cols(s),linewidth=0.1)
        #@df df StatsPlots.boxplot!(str,cols(s),fillalpha=0.75, linewidth=0.25,legend=false);
        title!(ti)
        #	@df df StatsPlots.boxplot!(str,cols(ln),fillalpha=0.75, linewidth=0.25,legend=false);
        #	@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.75,marker=(:black,stroke(1)),legend=false)
        end
    end

    