# 
function hydromon(x::Union{Regex,String,DataFrame}, col = 1, leg = :outertopright, logy = false)
        if isa(x,DataFrame)
            df = (x)
        else
            df = waread(x)
        end
        
        ti = try
            DataFrames.metadata(df) |> only |> last |> basename
        catch
            @warn "No basename in metadata!"
            raw""
        end

        if names(df)[col] == "date"
            col = col + 1
            @info "date column is at position 1, skipping to col+1..."
        end
        

        df = select(df,Cols(col,:date))

        @info ("Selection:",names(df))
        colname = names(df[!,Not(:date)])|>only
        s = Symbol.(filter(x -> !occursin(r"date|year|month"i, x), names(df)))
        years = unique(Dates.year.(df[!, :date]))
        ylog = logy ? :log : :identity
    
        begin
            su = @rsubset df year.(:date)==years[1]
            su = monmean(su)
            @df su Plots.plot(:month,cols(s),
                label = years[1], 
                title=ti,
                legendtitle=colname,
                yaxis=ylog,
                legend=leg)
            for yr in years[2:end]
                su = @rsubset df year.(:date)===yr
                su = monmean(su)
                @df su Plots.plot!(:month,
                    cols(s),
                    yaxis = ylog,
                    label = yr)
            end
            month_abbr = ["Jan", "Feb", "Mär", "Apr", 
                        "Mai", "Jun", "Jul", "Aug", "Sep", 
                            "Okt", "Nov", "Dez"];
            #xticks!(0.5:11.5 , month_abbr)
            xticks!(1:12, month_abbr)
            plot!()
        end
    end

    """
    hydrgraph plot, selection of first column by default
    x::Union{Regex,String,DataFrame}; col = 1,
        leg = :outertopright, logy = false)
    """
    