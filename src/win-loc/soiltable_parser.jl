#I can try to write a Reader function in Julia for this kind of data. Here is one possible way:


    # define a struct to store the data for each horizon
    struct Horizon
        method::String
        ThicknessScaling::Int
        horizon::Vector{Int}
        Name::Vector{String}
        ksat::Vector{Float64}
        k_recession::Vector{Float64}
        theta_sat::Vector{Float64}
        theta_res::Vector{Float64}
        alpha::Vector{Float64}
        Par_n::Vector{Float64}
        Par_tau::Vector{Float64}
        thickness::Vector{Float64}
        layers::Vector{Int}
    end
    
    # define a Reader function that takes a file name as input and returns a vector of Horizon structs
    function Reader(filename::String)
        # open the file for reading
        open(filename, "r") do file
            # initialize an empty vector to store the horizons
            horizons = Horizon[]
            # loop through each line of the file
            for line in eachline(file)
                # split the line by spaces and remove the curly braces
                tokens = filter(x -> x != "{" && x != "}", split(line," "))
                # get the index of the first equal sign
                i = findfirst(==("="), tokens)
                # get the method name from the first token
                method = tokens[1][1:i-1]
                #method = tokens[1:i-1]
                # get the ThicknessScaling value from the second token
                ThicknessScaling = parse(Int, tokens[2][i+1:end])
                # initialize empty vectors to store the values for each field
                horizon = Int[]
                Name = String[]
                ksat = Float64[]
                k_recession = Float64[]
                theta_sat = Float64[]
                theta_res = Float64[]
                alpha = Float64[]
                Par_n = Float64[]
                Par_tau = Float64[]
                thickness = Float64[]
                layers = Int[]
                # loop through the remaining tokens and parse the values according to the field name
                for token in tokens[3:end]
                    # get the index of the equal sign
                    i = findfirst(==("="), token)
                    # get the field name and value from the token
                    field = token[1:i-1]
                    value = token[i+1:end]
                    # parse the value according to the field name and append it to the corresponding vector
                    if field == "horizon"
                        append!(horizon, parse.(Int, split(value)))
                    elseif field == "Name"
                        append!(Name, split(value))
                    elseif field == "ksat"
                        append!(ksat, parse.(Float64, split(value)))
                    elseif field == "k_recession"
                        append!(k_recession, parse.(Float64, split(value)))
                    elseif field == "theta_sat"
                        append!(theta_sat, parse.(Float64, split(value)))
                    elseif field == "theta_res"
                        append!(theta_res, parse.(Float64, split(value)))
                    elseif field == "alpha"
                        append!(alpha, parse.(Float64, split(value)))
                    elseif field == "Par_n"
                        append!(Par_n, parse.(Float64, split(value)))
                    elseif field == "Par_tau"
                        append!(Par_tau, parse.(Float64, split(value)))
                    elseif field == "thickness"
                        append!(thickness, parse.(Float64, split(value)))
                    elseif field == "layers"
                        append!(layers, parse.(Int, split(value)))
                    else
                        error("Unknown field: $field")
                    end
                end
                # create a Horizon struct from the parsed values and push it to the horizons vector
                horizon = Horizon(method, ThicknessScaling, horizon, Name, ksat,
                                  k_recession, theta_sat, theta_res, alpha,
                                  Par_n, Par_tau, thickness, layers)
                push!(horizons, horizon)
            end
            # return the horizons vector
            return horizons
        end
    end
    
    # test the Reader function with a sample file
    pt="/mnt/d/temp/saale/control/thulba_soiltable.txt"

    horizons = Reader(pt)

    
    #This script will read the data from a file named "sample.txt" and return a vector of Horizon structs. You can access the fields of each struct by using dot notation. For example:
    
    
    horizons[1].method # returns "MultipleHorizons"
    horizons[2].ksat # returns [1.33969641503796e-6 1.33969641503796e-6 1.33969641503796e-6 1.83888680967466


  

    using CSV
    using DataFrames
    
    # Define a struct to store the data for each horizon
    struct Horizon
        method::String
        ThicknessScaling::Int
        horizon::Vector{Int}
        Name::Vector{String}
        ksat::Vector{Float64}
        k_recession::Vector{Float64}
        theta_sat::Vector{Float64}
        theta_res::Vector{Float64}
        alpha::Vector{Float64}
        Par_n::Vector{Float64}
        Par_tau::Vector{Float64}
        thickness::Vector{Float64}
        layers::Vector{Int}
    end
    
    # Define a function to parse a line of data into a Horizon object
function parse_line(line::String)
        # Split the line by spaces and remove the curly braces
        tokens = filter(x -> x != "{" && x != "}", split(line))
        
        # Extract the method and ThicknessScaling values
        method = tokens[1]
        ThicknessScaling = parse(Int, tokens[3])
        
        # Initialize empty vectors to store the rest of the values
        horizon = Int[]
        Name = String[]
        ksat = Float64[]
        k_recession = Float64[]
        theta_sat = Float64[]
        theta_res = Float64[]
        alpha = Float64[]
        Par_n = Float64[]
        Par_tau = Float64[]
        thickness = Float64[]
        layers = Int[]
        
        # Loop through the tokens and append the values to the corresponding vectors
        for i in 4:length(tokens)
            token = tokens[i]
            if token == "horizon"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(horizon, parse(Int, tokens[i]))
                end
            elseif token == "Name"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(Name, tokens[i])
                end
            elseif token == "ksat"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(ksat, parse(Float64, tokens[i]))
                end
            elseif token == "k_recession"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(k_recession, parse(Float64, tokens[i]))
                end
            elseif token == "theta_sat"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(theta_sat, parse(Float64, tokens[i]))
                end
            elseif token == "theta_res"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(theta_res, parse(Float64, tokens[i]))
                end
            elseif token == "alpha"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(alpha, parse(Float64, tokens[i]))
                end
            elseif token == "Par_n"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(Par_n, parse(Float64, tokens[i]))
                end
            elseif token == "Par_tau"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(Par_tau, parse(Float64, tokens[i]))
                end
            elseif token == "thickness"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(thickness, parse(Float64, tokens[i]))
                end
            elseif token == "layers"
                i += 1 # Skip the equal sign
                while tokens[i+1] != ";"
                    i += 1
                    push!(layers, parse(Int, tokens[i]))
                end
            else 
                continue # Skip any other token that is not a field name or value 
            end 
        end
        
        # Return a Horizon object with the extracted values 
        return Horizon(method, ThicknessScaling, horizon, Name, ksat,
                       k_recession, theta_sat, theta_res,
                       alpha,Par_n, Par_tau, thickness, layers)
end


function read_file(filename::String)
        # Read the file as a CSV file with space delimiter and no header
        df = CSV.File(filename, delim=' ', header=false) |> DataFrame
        
        # Initialize an empty vector to store the Horizon objects
        horizons = Horizon[]
        
        # Loop through each row of the DataFrame and parse it into a Horizon object
        for row in eachrow(df)
            # Get the line of data as a string
            line = row[1]
            
            # Parse the line into a Horizon object
            horizon = parse_line(line)
            
            # Append the horizon to the vector
            push!(horizons, horizon)
        end
        
        # Return a DataFrame with one column of Horizon objects
        return DataFrame(Horizon=horizons)
end

zz = read_file(pt)

filename=pt

line = df[1,:]
line = string(line)
parse_line(line)


using CSV, DataFrames

function parse_soiltable(pt::AbstractString)
data = []
open(pt, "r") do file
    for line in eachline(file)
        line = strip(line)
        if !isempty(line)
            parts = split(line, ';')
            values = Dict{String, String}()
            for part in parts
                key_value = split(strip(part), '=')
                if length(key_value) >= 2
                    key = strip(key_value[1])
                    value = strip(key_value[2])
                    values[key] = value
                end
            end
            push!(data, values)
        end
    end
end
return(data)
end




pt="/mnt/d/temp/saale/control/thulba_soiltable.txt"
data = parse_soiltable(pt)
dfs=[]
for i in 1:size(data)[1]
    push!(dfs,DataFrame(data[i]))
end
#dx = map(x->x, dfs)
df = DataFrame(data[4])
permutedims(df,1,split(names(df)[1]," ")[1]) #or

df = permutedims(df,1,replace(names(df)[1], r".{.*"=>""))
rename!(df,2=>"data")

#map(x->parse(Float64,x),df.data[Not(1)])



replace(names(df)[1],r".{method "," ")



println(df)
