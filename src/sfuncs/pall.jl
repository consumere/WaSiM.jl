# reduces + merges by date
function pall(files::Vector{DataFrame};toyr=true,leg=:outertopright)
        #files = dfs
        bns = try
            broadcast(x->DataFrames.metadata(x)|>only|>last|>basename,files)
        catch
            @warn "No basenames in metadata!"
            raw""
            end
    
        try
            for i in 1:length(files)
                s = Symbol.(filter(x->!occursin(r"year|date",x),names(files[i])))
                nn = bns[i]|>x->split(x,".")|>first
                for x in s
                    #println(string(x)*"-"*nn)
                    newname = string(x)*"-"*nn
                    #rename!(i,s[x]=>s[x]*"-"*bns[x])
                    rename!(files[i],x=>newname)
                end
            end
        catch
            @warn "error in renaming!"
            @debug begin
                nms=map(x->names(x),files)
                "names are: $nms"
            end
            #@warn "error in renaming!"
        end
    
        if toyr==true
            dfs = broadcast(x->yrsum(x),files)
            df = reduce((left, right) -> 
            innerjoin(left, right, on = :year,
            makeunique=true
            ),dfs)
        else
            df = reduce((left, right) -> 
            innerjoin(left, right, on = :date,
            makeunique=true
            #renamecols = lowercase => uppercase
            ),files)
        end
        
        
        ti = try
            #DataFrames.metadata(df)|>only|>last|>basename
            z = map(x->split(x,".")|>first,bns)
            if length(z)>5
                z = z[1:5]
                "Merged Series"
            end
            #"Series of "*join(z," ")
            join(z," ")
        catch
            @warn "No basenames for title"
        ti = raw""
            end
    
        if length(ti)>30
            ti = ti[1:30]*"..."
        end
    
        if all_values_equal(df[!,Not(Cols(r"date|year|month|day"))])==true
            @error "all values are equal!"
            return
        end
    
        if (any(x->occursin("year",x),names(df)))
            s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
            p = @df df Plots.plot(:year,cols(s),legend = leg, title=ti)
        else    
            s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
            p = @df df Plots.plot(:date,cols(s),legend = leg, title=ti)
        end
        return p
    end

    """reduces + merges by date"""
    