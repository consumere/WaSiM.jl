# 
function ctl3()
        # Loop through the current directory and its subdirectories
        matches::Vector{Any} = []
        for (root, dirs, files) in walkdir(".")
        # Loop through each file name
        for file in files
            # If the file name ends with .xml
            if endswith(file, ".xml")
            # Join the root and file name to get the full path
            path = joinpath(root, file)
            # Open the file for reading
            open(path) do f
                # Loop through each line of the file
                for line in eachline(f)
                # If the line contains 'compiling symbols in control file '
                if occursin("compiling symbols in control file ", line)
                    # Split the line by whitespace and get the fields from index 9 to 15
                    fields = split(line)[8:end] #," "
                    # Join the fields by space and print them
                    println(join(fields, " "))
                    out=join(fields, " ")
                    push!(matches,out)
                end
                end
            end
            end
        end
        end
        fl = last(matches)
        fl = split(fl)|>last
        fl = split(fl,"\"")|>first
        if !occursin("Wasim",fl)
            fl = replace(fl,"control"=>"D:/Wasim/regio/control")
        end
        
        return(string(fl))
    end

    """
    checks state of wasim routing table
    isroute(ctl3()) 
    isroute(filename::AbstractString;match=r"timeoffset")
    """
    