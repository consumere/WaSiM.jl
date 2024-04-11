# 
function build_soil_dictionary(soiltable::String)     
        entries = split(soiltable, "}")
        filter!(s -> !isempty(s), entries)
        dictionary = Dict{Int, DataFrame}()
        
        for entry in entries
            # Extract key and value
            key = parse(Int, match(r"(\d+)\s*\{", entry).captures[1]|>strip)
            println("check: $key")
            value = match(r"\{\s*(.*)", entry).captures[1]|>strip
            
            # Process value string and convert to DataFrame
            value = replace(value, "method = MultipleHorizons;  EvapMaxDepth = 0.15;" => "")
            value = replace(value, r";" => ",")
            value = replace(value, r";" => "")
            value = replace(value, r"\$" => "")
            value = replace(value, r"#" => "")
            value = replace(value, r"\s+" => " ")
            value = replace(value, r"^\s+|\s+$" => "")
            value = replace(value, r",$" => "")
            pairs = split(value, ",")
            # Create an empty dictionary to store the column names and values
            #col_dict = Dict{String, Vector{String}}()
    
            # Iterate over the key-value pairs and extract the column names and values
            filter!(s -> !isempty(s), pairs)
            z=[]
            for pair in pairs
                key2, value = split(strip(pair), " = ")
                push!(z,DataFrame(key2 => value))
            end
            # Create the DataFrame using the extracted column names and values
            df = hcat(z...)
            # Store key-value pair in the dictionary
            dictionary[key] = df
        end
        
        return dictionary
    end

    """
    soiltable reader wrapper func.
    """
    