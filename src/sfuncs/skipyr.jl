# with DataFrame input
function skipyr(df::DataFrame;fun=nothing)
        fst = year.(df.date)[1]
        data = filter(:date => x -> Dates.year(x) > fst, df)
        #return fun ? fun(data) : data #if boolean
        if isnothing(fun)
            return data
        else
            return fun(data)
        end
    end


    
    """with DataFrame input"""
    