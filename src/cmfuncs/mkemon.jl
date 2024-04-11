# 
function mkemon(x::Union{String,Regex,DataFrame};col::Any="tot_average",msk=false)
        if x isa String
            printstyled("reading $x\n",color=:light_red)
            df = dfr(x)
            dropmissing!(df)    
        elseif x isa Regex
            x = first(dfonly(x))
            printstyled("reading $x\n",color=:light_red)
            df = dfr(x)
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
        #check if "tot_average" in df
        # col = try 
        #     findfirst(x->occursin(r"tot_average",x),names(df))
        # catch
        #     @warn "No :tot_average column in dataframe - changed to first column!"
        #     first(names(df))
        # end

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
        
        try 
            select!(df, Cols(:date,col))
        catch
            @warn "No :date & $col column in dataframe!"
            return
        end


        if msk
            #only values > 0 
            df = try 
                DataFrames.subset(df,col => ByRow(>(0)))
            catch
                cn=propertynames(df[!,Not(Cols(r"date|month"i))])|>first
                #DataFrames.index(df,col)
                #DataFrames.subset(df,Int(col) => ByRow(>(0)))
                DataFrames.subset(df,cn => ByRow(>(0)))
            end
        end


        df[!, :month] = month.(df[!,:date]);
        grp = groupby(df, :month)

        months = ["Januar", "Februar", "März", "April",
            "Mai", "Juni", "Juli", "August", "September",
            "Oktober", "November", "Dezember"]
        
        f = Figure()
        Axis(f[1, 1], title = ti*" Basin:"*string(col),
            yticks = ((1:12) ./ 4,  reverse(months)))
        for i in 12:-1:1
            values = select(grp[i], Not(Cols(r"date|month"i)))|>Matrix|>vec
            d = density!(values, offset = i / 4,
                        color = :x, colormap = :thermal, 
                        #colorrange = (-5, 5),
                        strokewidth = 1, strokecolor = :black)
            #Apply an absolute translation to the Scene
            translate!(d, 0, 0, -0.1i) 
        end
        return f
    end

    # """
    # uses xarray to read netcdf and makie to plot.
    # using PyCall
    # @pyimport xarray as xr
    # """
    # 