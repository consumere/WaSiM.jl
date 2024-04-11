# 
function fsoil(fn::String)
        soiltable = open(fn) do io
            a = readbetween(io, "soil_table", "substance_transport")
            return(join(a[3:end-1]," "))
        end
        mysoildict = build_soil_dictionary(soiltable)
        xdf = DataFrame()
        # Iterate over the key-value pairs in the dictionary
        for (key, df) in mysoildict
            # Add a column named "key" with the current key value to the DataFrame
            df.key = fill(key, size(df, 1))    
            # Append the current DataFrame to the combined DataFrame
            append!(xdf, df)
        end
        #propertynames(xdf)|>cb
        nms=[:horizon, :ksat, :theta_res, :theta_sat, :alpha, :Par_n, :thickness, :maxratio]
        for col in nms
            xdf[!, col] .= [parse.(Float64, split(string(x))) for x in xdf[!,col]]
        end
        xdf.sums = [sum(x) for x in xdf.thickness]
        return xdf
    end

    """
    [parse.(Float64, split(string(x))) for x in nm]
    """
    