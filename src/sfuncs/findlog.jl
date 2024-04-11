# 
function findlog(;lb=.4)
        v = qall(;recursive=true)
        #map(x->DataFrames.subset(x,3 => ByRow(<(1))),v)
        k = try 
            map(x->DataFrames.subset(x,3 => ByRow(>(lb))),v)
            catch
                @error "no df on lowerbound! "
                return
        end
        df = try 
            reduce(vcat,k) 
            catch
                @error "vcat failed! "
                return
        end
        
        df = try 
            DataFrames.subset(df,3 => ByRow(<(1.0)))
            catch
                println(first(k))
                @error "no df on upperbound! "
                return
        end
        
        
        return df
    end

    """
    skips first year and returns df
    optionally apply fun to df
    """
    