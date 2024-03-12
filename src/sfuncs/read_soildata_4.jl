# 
function read_soildata_4(filename::String)
        file = open(filename) do io
            a = readbetween(io, "soil_table", "special_output")
            return(a)
        end
        file = file[2:end-1]
        #data = Grep.grep(r"^[0-9].* {",file)
        data = Grep.grep(r"^[0-9]|^[ ;]|^[;]",file)
        filter!(x->length(x)>3,data)
        return data
    end
    
    import InteractiveUtils: clipboard
    