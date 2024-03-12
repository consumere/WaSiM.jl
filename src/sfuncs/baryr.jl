# 
function baryr(df::DataFrame)
        s = filter(x -> !(occursin(r"year|date", x)), names(df))
        for x in s
            newname = replace(x, "_1" => "")
            rename!(df, Dict(x => newname))
        end
        s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
        # Create a grouped bar plot for each column    
        # StatsPlots.groupedbar(df.year,[df[!, col] for col in s], 
        # #group = s,#repeat(s, outer = size(df, 1)),
        # xlabel = "Year", ylabel = "Value", title = "Grouped Bar Plot")
        ti = try
            z=DataFrames.metadata(df)|>only|>last|>basename
            basename(pwd())*" $z"
        catch
            @warn "No basename in metadata!"
            ti = "Series of "*basename(pwd())
        end
        @df df groupedbar(df.year,cols(s), 
        #group = s,#repeat(s, outer = size(df, 1)),
        legend = :outertopright,
        xticks = df.year,
        xrotation = 45,
        xlabel = "", ylabel = "[mm]", title = ti)
    end

    