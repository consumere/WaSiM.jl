# 
function baryrmean(df::DataFrame;leg=:outertopright)
        """
        automatically sums if only datecolumn is available
        """
        v = try map(
            (x->occursin(r"date", x) & !occursin(r"year", x)),
            (names(df))
            )
        catch
            @error "No date available!"
            return
        end

        if any(v)
            df = yrmean(df)
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
        legend = leg,
        xticks = df.year,
        xrotation = 45,
        xlabel = "", ylabel = "[mm]", title = ti)
    end

    """
    kge barplot 
    """
    