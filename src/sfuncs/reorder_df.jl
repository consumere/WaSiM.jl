# 
function reorder_df(df::DataFrame,tofirst::Bool=false)
        if tofirst
            df = select(df, :date, :)
        else
            df = select(df, Not(:date), :)
        end
        #df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        return(df)
    end

    # 