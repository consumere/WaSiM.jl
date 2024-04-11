# 
function dfmo(x::Union{Regex,String,DataFrame};leg = :topright, fun=monmean,mode=:line)
        if isa(x,DataFrame)
            df = (x)
        else
            df = waread(x)
        end
            
        
        #df = waread(mm)
        #DataFrames.metadata(df)|>collect|>only|>last
        ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
            ti = raw""
        end
        
        if fun==false
            #reorder           
            dx = hcat(df[:,Cols(r"date")],df[!,Not(Cols(r"date"))])
            #printstyled("no 