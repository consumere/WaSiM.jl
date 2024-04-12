
filename=infile
npp(filename)

#function rsd(filename::String)
    data = open(filename) do io
        a = readbetween(io, "soil_table", "substance_transport")
        return(a[3:end-1]) # remove first 2 and last line
    end
    
    data = broadcast(x -> replace(x,    
        r"^#.*" => "",
        r"^[[].*" => "",
        r"method" => "",
        r"MultipleHorizons" => "",
        #r"}" => "",
        r" = " => "=" ), data)

    filter!(s -> !isempty(s), data)
#    da = hcat(data)|>string
    da = join(data, " ")

    key_value_pairs = split(da, ";")

    da = broadcast(x -> replace(x,    
        r"{=" => "",
        r";\n" => "; " ), )

    result = []
    # iterate over each line
    for line in data
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
#    return df
#end  



using DataFrames

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
        # for pair in pairs
        #     key2, value = split(strip(pair), " = ")
        #     push!(col_dict[key2], split(value))
        # end


        z=[]
        for pair in pairs
            key2, value = split(strip(pair), " = ")
            push!(z,DataFrame(key2 => value))
        end

        # Create the DataFrame using the extracted column names and values
        df = hcat(z...)
        #value = replace(value, r"(\w+)\s*=" => s"'\1'=")
        #value = replace(value, r"(\w+)\s*=" => s"'\1'=")
        #nm,dat = split(value, "=")
        #df = DataFrame(eval(Meta.parse("[$value]")))
        #df = DataFrame(value,:auto)
        
        # Store key-value pair in the dictionary
        dictionary[key] = df
    end
    
    return dictionary
end

filename = raw"D:\Wasim\sinn\control\s100-spin1_loc.ctl"
data = open(filename) do io
    a = readbetween(io, "soil_table", "substance_transport")
    return(a[3:end-1]) # remove first 2 and last line
end
soiltable = join(data, " ")
mysoildict = build_soil_dictionary(soiltable)

mysoildict[9050]


combined_df = DataFrame()
# Iterate over the key-value pairs in the dictionary
for (key, df) in mysoildict
    # Add a column named "key" with the current key value to the DataFrame
    df.key = fill(key, size(df, 1))
    
    # Append the current DataFrame to the combined DataFrame
    append!(combined_df, df)
end

combined_df
combined_df.key
combined_df[6,3]

sel = select(combined_df, [:ksat,:key])
sel.ksat = [parse.(Float64, split(string(x))) for x in sel.ksat]

p1 = plot()
for i in eachrow(sel[1:25,:])
    println(i[1]," -> ",i[2])
    scatter!(i[1],label=i[2])
end
p1

#@df sel[1:20,:] plot(:key, :ksat, group=:key, legend=:topleft)
#@df sel groupedbar(sel.key,cols(:ksat), legend = :outertopright)



filename = raw"D:\Wasim\regio\control\rcm200_r9-cl4.ctl"
data = open(filename) do io
    a = readbetween(io, "soil_table", "substance_transport")
    return(a[3:end-1]) # remove first 2 and last line
end
soiltable = join(data, " ")
mysoildict = build_soil_dictionary(soiltable)
xdf = DataFrame()
# Iterate over the key-value pairs in the dictionary
for (key, df) in mysoildict
    # Add a column named "key" with the current key value to the DataFrame
    df.key = fill(key, size(df, 1))    
    # Append the current DataFrame to the combined DataFrame
    append!(xdf, df)
end

function fparse(nm)
    [parse.(Float64, split(string(x))) for x in nm]
end
k = hcat(xdf.key,fparse(xdf.thickness))
k = DataFrame(k,:auto)
k.x2[1]|>sum
k.sums = [sum(x) for x in k.x2]
k

sel = select(xdf, [:thickness,:key])
sel.thickness = [parse.(Float64, split(string(x))) for x in sel.thickness]
sel.sums = [sum(x) for x in sel.thickness]
@df sel groupedbar(:sums, cols(:key), 
        group=:sums, 
        legend=:outertopleft)


#wrapper func:
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
    #names(xdf)|>println
    #propertynames(xdf)|>cb
    nms=[:horizon, :ksat, :theta_res, :theta_sat, :alpha, :Par_n, :thickness, :maxratio]
    for col in nms
        xdf[!, col] .= [parse.(Float64, split(string(x))) for x in xdf[!,col]]
    end
    # xdf.thickness .= [parse.(Float64, split(string(x))) for x in xdf.thickness]
    # xdf. .= [parse.(Float64, split(string(x))) for x in xdf.thickness]
    xdf.sums = [sum(x) for x in xdf.thickness]
    return xdf
end

filename = raw"D:\Wasim\regio\control\rcm200_x22-cl4.ctl"
okd = fsoil(filename)
