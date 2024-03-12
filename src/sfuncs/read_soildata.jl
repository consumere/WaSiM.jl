# 
function read_soildata(filename::String)
        data = open(filename) do io
            a = readbetween(io, "soil_table", "substance_transport")
            return(a)
        end
        data = broadcast(x -> replace(x,    
        r"^#.*" => "",
        r"^[[].*" => "",
        r"method" => "",
        r"MultipleHorizons" => "",
        r"}" => "",
        r" = " => "=" ), data)


        filter!(s -> !isempty(s), data)
        lines = data[2:end]
        result = []
        # iterate over each line
        for line in lines
            # skip empty lines
            if isempty(line)
                continue
            end
            
            # split the line into number and fields
            number, fields = split(line, " {", limit=2)
            
            # convert the number to an integer
            number = parse(Int, number)
            
            # split the fields into individual fields
            fields = split(fields, ';')
            fields = map(f->strip(f),fields)
            
            # initialize a dictionary to store the data for this line
            data = Dict{String, Any}()
            
            # store the number in the dictionary
            data["number"] = number
            
            # iterate over each field
            try
            for field in fields
                # check if the field contains the " = " substring
                if occursin("=", field)
                    # split the field into key and value
                    key, value = split(field, "=")
                    
                    # check if the key is "Name"
                    if key == "Name"
                        # keep the value as a string
                    else
                      
                        # check if the value is a number
                        if occursin(r"^-?\d+(\.\d+)?$", value)
                            # convert the value to a float
                            value = parse(Float64, value)
                        elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                            # convert the value to a float (scientific notation)
                            value = parse(Float64, value)
                        elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                            # convert the value to a float (scientific notation)
                            value = parse(Float64, value)
                        elseif occursin(r"^\d+ \d+", value)
                            # convert the value to an array of integers
                            value = parse.(Int, split(value))
                            #value = parse.(Float64, split(value))
                        elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)

                            # convert the value to an array of floats
                            value = parse.(Float64, split(value))
                                                    
                        end

                    end
                    
                    # store the key-value pair in the dictionary
                    data[key] = value
                end
            end
        catch
            @warn "could not parse $value"
            continue
        end
            
            # append the dictionary to the result array
            push!(result, data)
        end
        df = DataFrame(result)
        df.ksat = [parse.(Float64, split(string(x))) for x in df.ksat]
        return df
    end  

    