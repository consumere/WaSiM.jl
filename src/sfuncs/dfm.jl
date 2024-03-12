# 
function dfm(x::Union{Regex, String, DataFrame}; 
        ann = true, 
        log = false, 
        title = true,
        leg = false,  
        fun = monmean, 
        mode=:line)
        if isa(x, DataFrame)
            df = x
        else
            df = waread(x)
        end
        
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
            ti = raw""
        end

        ln = Symbol.(filter(x -> !occursin(r"date|month|year", x), names(df)))
        
        if fun == false
            # Reorder the DataFrame
            dx = hcat(df[:, Cols(r"date")], df[!, Not(Cols(r"date"))])
            @info "No 