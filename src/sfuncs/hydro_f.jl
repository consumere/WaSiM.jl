# 
function hydro_f(df::DataFrame; leg = :outertopright, logy=false)
        ti = try
            DataFrames.metadata(df) |> only |> last |> basename
        catch
            @warn "No basename in metadata!"
            raw""
        end
        #df = selt(df,4)
        s = Symbol.(filter(x -> !occursin(r"date|year|month"i, x), names(df)))
        years = unique(Dates.year.(df[!,:date]))
        #Scale of the axis. 
        #Choose from [:identity, :ln, :log2, :log10, :asinh, :sqrt].
        if logy 
            ylog = :log
        else
            ylog = :identity
        end

            su = @rsubset df year.(:date)==years[1]
            #unique(monthname.(su.date))
            tm_ticks = round.(su.date, Month(1)) |> unique
    
            hp1 = @df su Plots.plot(:date,cols(s),
                yaxis = ylog,
                xticks=(tm_ticks, 
                    #Dates.format.(tm_ticks, "uu/yyyy")), 
                    Dates.format.(tm_ticks, "uu")), 
                xrot=45, xminorticks=true, 
                xlim=extrema(su.date),
                label = (years[1]),  #first
                title=ti,
                formatter=:plain,
                legend = leg
                );
    
    
            for yr in years[2:end]
                su = @rsubset df year.(:date)===yr
                #tm_ticks = round.(su.date, Month(1)) |> unique; 
                hp1 = Plots.plot!(
                    vec(Matrix(su[!,Not(:date)])),
                    label = (yr),
                    formatter=:plain)
            end
    
            return hp1
        
    end

    """
    hydrgraph monthly mean plot
    hydromon(x::Union{Regex,String,DataFrame}, col = 1, leg = :outertopright, logy = false)
    """
    