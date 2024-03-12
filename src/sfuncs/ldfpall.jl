# 
function ldfpall(x::Regex)
        files = rglob(x)
        dfs = DataFrame[]
        for file in files
            if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
            file_path = file
        println("reading ",file_path,"...")
        p1 = readdf(file_path)
        ##renamer
        for x in 1:size(p1,2)-1
            rename!(p1,x=>basename(file_path)*names(p1)[x])
        end
        push!(dfs, p1)
            end
        end
        df = reduce((left, right) -> 
        innerjoin(left, right, on = :date,makeunique=true),dfs)
        ##to preserve column order and names (date at last position)
        #df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])  
        y = filter(x->!occursin("date",x), names(df))
        s = map(y -> Symbol(y),y)
        @df df Plots.plot(:date,
                cols(s), yaxis = :log,
                legend = :outertopright)
    end

    