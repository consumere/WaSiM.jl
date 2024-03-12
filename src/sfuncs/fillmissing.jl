# 
function fillmissing(df,fillval=-9999)
        for col in eachcol(df)
            col .= coalesce.(col, fillval)
        end
        return df
    end

    """
    build_soil_dictionary from soiltable
    see fsoil
    soiltable = open(fn) do io
        a = readbetween(io, "soil_table", "substance_transport")
        return(join(a[3:end-1]," ")) #rem first 2 and last lines
    end
    """
    