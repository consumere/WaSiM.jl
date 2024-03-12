# 
function dfrib(x::Union{String,Regex,DataFrame};
        col::Any=1,confidence_level::Number=0.95)
        if x isa String
            printstyled("reading $x\n",color=:light_red)
            df = waread2(x;silencewarnings=false)
            dropmissing!(df)    
        elseif x isa Regex
            x = first(dfonly(x))
            printstyled("reading $x\n",color=:light_red)
            df = waread2(x;silencewarnings=false)
            dropmissing!(df)
        else
            df = x #reorder_df(x)
            dropmissing!(df)
        end

        ti = try
            DataFrames.metadata(df) |> only |> last |> basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end       

        if !("date" in names(df))
            @error "no :date column found in $x !"
            return
        end

        if col isa Number && names(df)[col] == "date"
            if col + 1 >= ncol(df)
                col = col - 1
                @info "date column is at $col position, skipping to col-1..."
            elseif names(df)[col - 1] == "date"
                col = ncol(df)
                @info "date column is at $col-1 position, skipping to last col..."
            else
                col = col + 1
                @info "date column is at position $col, skipping to col + 1..."
            end
        end

        #colname = string(names(df)[1]) #Not(:date)
        #colname = string(names(df[!, Not(:date)]))
        ##has to be in this order...
        colname = only(names(select(df,Cols(col))))
        
        if length(ti)>1
            lab = ti[1:4]*colname
        else
            lab = colname
        end

        lab = replace(lab,r"_" => " ")
        ti = replace(ti,r"_" => " ")

        df = select(df,Cols(col,:date))
        @info "aggregation monthly mean of col: $col at confidence_level: $confidence_level"
        
        # if startswith(colname,"_")
        #     colname = "basin"*colname
        # end

        # #dmean = copy(df)
        # dmean = df
        # s = map(Symbol, filter(x -> !occursin("date", x), names(dmean)))
        # dmean[!, :month] = month.(dmean[!,:date])
        # select!(dmean, Not(:date))
        
        # # # Calculate monthly mean and confidence interval using combine and groupby
        # # dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
        # #     s .=> (dmean -> (mean(dmean), 1.96 * std(dmean) / sqrt(length(dmean)))) 
        # #     .=> s)

        # # Calculate monthly mean and narrower confidence interval using combine and groupby
        # dmean = DataFrames.combine(DataFrames.groupby(dmean, :month),
        #     s .=> (
        #         dmean -> 
        #         (mean(dmean), 1.1 * std(dmean) / sqrt(length(dmean)))) 
        #     .=> s)
    
    
        # # Create separate columns for mean, lower bound, and upper bound
        # for col in s
        #     dmean[!, string(col)*"_mean"]  .= [first(only(x)) for x in eachrow(dmean[!, col])]
        #     dmean[!, string(col)*"_lower"] .= [first(only(x)) - last(only(x)) for x in eachrow(dmean[!, col])]
        #     dmean[!, string(col)*"_upper"] .= [first(only(x)) + last(only(x)) for x in eachrow(dmean[!, col])]
        # end
        
        dmean = monc(df;confidence_level = confidence_level) #s.o. confidence_interval 95%
        mean_col = propertynames(select(dmean, Cols(r"mean")))[1]
        lower_col = propertynames(select(dmean, Cols(r"lower")))[1]
        upper_col = propertynames(select(dmean, Cols(r"upper")))[1]
        lower_bounds = dmean[!, lower_col]|>unique
        upper_bounds = dmean[!, upper_col]|>unique

        #plt = @df dmean plot(:month, Matrix(dmean[!,Cols(r"mean")]), 
            #ribbon=(dmean.Lai_lower, dmean.Lai_upper), 
            # ribbon=(Matrix(dmean[!,Cols(r"lower")]), 
            #         Matrix(dmean[!,Cols(r"upper")])), 
        plt = @df dmean plot(:month, cols(mean_col), 
                          ribbon=(lower_bounds, upper_bounds),
            labels = lab, #see above at ti   # "Lai_mean", 
            xlabel = "", #"Month", 
            ylabel = "",  #"Lai", 
            title = ti, #"Lai Mean with Ribbon",
            legend = :topleft,
            xticks = (1:12, 
            [ monthabbr(x) for x in unique(dmean.month) ]),
            xrotation = 35
            )    
        return plt
    end

    """
    reader for earth engine csv files
    """
    