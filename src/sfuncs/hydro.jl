# 
function hydro(x::Union{Regex,String,DataFrame}; col::Union{Int64,Regex,String} = 1,
        leg = :outertopright, logy = false)
        if isa(x,DataFrame)
        df = copy(x)
        else
        df = waread(x)
        end

        ti = try
            DataFrames.metadata(df)|>values|>only|>basename
        catch
            @warn "No basename in metadata!"
            raw""
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

        
        # if isa(col,Regex)
        #     @info "grep first match of $col ..."
        #     col = filter(x->occursin(col,x),names(df))[1]
        # end
        
        try
            select!(df,Cols(col,:date))
        catch
            @error "col not found!"
            return
        end

        if length(names(df)) < 2
            col = string.(col)
            @error "no $col column found!"
            #cn = names(x)
            #@info "p $nc"
            return
        end


        @info ("Selection:",names(df))

        colname = names(df[!,Not(:date)])|>only
        #ti = colname*" - "*ti

        s = filter(x -> !occursin(r"(?i)date|year|month", 
        string(x)), names(df))
        years = unique(year.(df.date))
        years_str = string.(years)

        ylog = logy ? :log : :identity
        
        mn = [ monthabbr(x) for x in unique(month.(df.date)) ]

        hp1 = plot(#xticks = mn,
               xrot = 45,
               xminorticks = false,
               xmajorticks = false,
               yaxis = ylog,
               #xlim = extrema(df.date),
               title = ti,
               legendtitle = colname,
               #formatter = :plain,
               legend = leg)
        #title!("Legende", legendtitle = "Titel")

        for yr in years #[2:end]
            su = filter(row -> year(row.date) == yr, df)
            hp1 = plot!(vec(Matrix(select(su, Not(:date)))),
                        label = yr) #,formatter = :plain)
        end
        lng = size(df,1) ./ length(years) ./ 12
        nln = size(df,1) ./ length(years)
        st = 15:lng:nln
        xticks!(st, mn)

        return hp1
    end

    cdinto = cdof

    """
    reads from stationdata, reprojects to EPSG
    #ProjString("+proj=longlat +datum=WGS84 +no_defs")
    """
    