# 
function ctlook(infile::String)
        data = open(filename) do io
            readbetween(io, "soil_table", "substance_transport")
        end
        data = data[3:end]
        numbers = extract_numbers(data)
        # Remove comments and unnecessary characters
        data = broadcast(x -> replace(x,    
                r"^#.*" => "",
                r"^[[].*" => "",
                r"method" => "",
                r"MultipleHorizons" => "",
                r"}" => "",
                r" = " => "=",
                #r"[?*.{=;]" => "",      #problems
                r";" => "" ), data)
        data = Grep.grep(r"^(?! Evap.*$|^[0-9].*$)",data)
        data = strip.(data)
        filter!(x -> !isempty(x), data)
        sel = Grep.grep(r"^Name", data)
        # Split each string into "Name" and "Value"
        split_name_value = split.(sel, '=', limit=2)
        
        ks = split.(Grep.grep(r"^ksat", data), '=', limit=2)
        ths = split.(Grep.grep(r"^theta_sat", data), '=', limit=2)
        thr = split.(Grep.grep(r"^theta_res", data), '=', limit=2)
        thk = split.(Grep.grep(r"^thickness", data), '=', limit=2)
        parn = split.(Grep.grep(r"^Par_n", data), '=', limit=2)
        #ks = Dict(first(ks)[1] => getindex.(ks, 2))
        # Create a DataFrame
        df = DataFrame( id = numbers,
                        Name = getindex.(split_name_value, 2),
                        ksat = getindex.(ks, 2),
                        theta_sat = getindex.(ths, 2),
                        theta_res = getindex.(thr, 2),
                        Par_n = getindex.(parn, 2),
                        Thick = getindex.(thk, 2))
        return df
    end

    """
    reads controlfile
    select landuse table
    returns a Vector of DataFrames
    """
    