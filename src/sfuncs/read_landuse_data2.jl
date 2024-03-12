# 
function read_landuse_data2(filename::AbstractString)
        data = open(filename) do io
            readbetween(io, "landuse_table", "special_output")
        end
        data = broadcast(x -> replace(x,
        r"\\t" => "",        #strip comments
        r"^#.*" => "",        #strip comments
                            r";.*$" => "",         # match everything after semicolon until the end the line
                            r" = " => "=",      #strip spaces around "="
                            r";" => "" ), data)
                            #r"^[[].*" => "",      #strip module names
                            #r"method" => "",     #startflag
                            #r"}" => "",          #this is needed for end flag
                            
                            #r"[?*.{=;]" => "",      #problems


        class_dataframes = []
        io = open(filename, "r")

        for line in data
            if occursin(r"(?i)method", line)
                classmatch = strip(line[3:15])  # first 10 chars, stripping leading/trailing whitespaces
                println(classmatch)
                filtered_lines = readbetween(io, Regex(classmatch), r"}$")
                filtered_lines = replace.(filtered_lines, r";.*" => "")
                filtered_lines = filter(x -> !isempty(x), filtered_lines)
                # Add a header
                header = ["Parameter=Value"]
                filtered_lines = vcat(header, filtered_lines)
                # Create a CSV File from the filtered lines
                csv_file = CSV.File(IOBuffer(join(filtered_lines, '\n')), delim='=', quotechar=' ', ignorerepeated=true)
                # Convert the CSV File to a DataFrame
                class_df = DataFrame(csv_file)
                push!(class_dataframes, class_df)
            end
        end
        close(io)
        return class_dataframes
    end

    """
    scatterplot of landusedata
    dfs::Vector=dfs,x::String="rs_evaporation"
    """
    