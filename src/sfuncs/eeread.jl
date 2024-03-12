# 
function eeread(filename)
        df = CSV.read(filename, DataFrame)
        df = unique(df, 1) # remove duplicate date rows
        rename!(df, 1 => :date)
        df.date .= [Date(d, dateformat"u d, y") for d in df.date]
        z = basename(filename)
        DataFrames.metadata!(df, "filename", z, style=:note)
        return df
    end
    
    """
    usage like: kernelplot("route.txt") but with linlog option
    reads with  qba() and names of df or (route.txt)
    lin: cols=Cols(2,4)
    log: cols=Cols(3,5)
    """
    