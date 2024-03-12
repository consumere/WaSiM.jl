# 
function baryrsum(df::Regex)
        """
        automatically sums if only datecolumn is available
        """
        df = globdf(df)|>first|>readf
        v = map(
            (x->occursin(r"date", x) & !occursin(r"year", x)),
            (names(df))
            )

        if any(v)
            df = yrsum(df) 
        end
        
        s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
        ti = try
            z=DataFrames.metadata(df)|>only|>last|>basename
            basename(pwd())*" $z"
        catch
            @warn "No basename in metadata!"
            ti = "Series of "*basename(pwd())
        end
        @df df groupedbar(df.year,cols(s), 
        legend = :outertopright,
        xticks = df.year,
        xrotation = 45,
        xlabel = "", ylabel = "[mm]", title = ti)
    end

    