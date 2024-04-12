#import Grep:grep

"""
smarter way to find control file in xmls
"""
function ctl()
    # Print a message
    println("looks for control file in all xmls")
    # Loop through the current directory and its subdirectories
    matches::Vector{Any} = []
    if (any(x->isdir(x),readdir()))
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
                        out = join(fields, " ")
                        ol=out|>split|>last
                        println("check: $ol")
                        push!(matches,out)
                    end
                    end
                end
                end
            end
        end
        return(matches)
    else
        for file in (filter(x->isfile(x) && endswith(x, ".xml"),readdir()))
            # Join the root and file name to get the full path
            root = "."
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
                    out = join(fields, " ")
                    ol=out|>split|>last
                    println("check: $ol")
                    push!(matches,out)
                end
                end
            end
            end
        end
        return(matches)  
end

"""
windows path to wsl
"""
function towsl(file_path::Union{String, Nothing})
    if isnothing(file_path)
        #return error("File path cannot be nothing")
        drive = lowercase(pwd()[1])
        rst = replace(pwd()[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    else
        drive = lowercase(file_path[1])
        rst = replace(file_path[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    end
    return wsl_path
end


ctf = ctl()
fls = []
for i in ctf
    #tmp = i|>split|>last
    tmp = replace(last(split(i)),r"ctl.*"=>"ctl")
    
    if (!occursin("regio",tmp) && occursin("regio",pwd()))
        tmp = replace(tmp,
    #    r"ctl.*"=>"ctl",
        "control"=>"D:/Wasim/regio/control")
    elseif (!occursin("brend",tmp) && occursin("brend",pwd()))
        tmp = replace(tmp,
    #    r"ctl.*"=>"ctl",
        "control"=>"D:/Wasim/Tanalys/DEM/brend_fab/control")
    elseif (!occursin("temp",tmp) && occursin("saale",pwd()))
        tmp = replace(tmp,
    #    r"ctl.*"=>"ctl",
        "control"=>"D:/temp/saale/control")
    end

    if Sys.isunix()
        tmp = towsl(string(tmp))
    end

    if isfile(tmp)
        printstyled("found: $tmp\n", color=:yellow)
        push!(fls,tmp)
    else
        @warn "File not found: $tmp"
    end
        
end


function isspin(filename::AbstractString;match=r"1 lin")
    START = false

    function process_line(line::AbstractString)
        m = strip(line)
        if occursin(match, line)
            return m, true
        end
        return m, false
    end

    for (i, line) in enumerate(eachline(filename))
        # if i > 50 && if occursin(r"^\[SpinUp]", line)
        #     m, done = process_line(line)
        #     printstyled("$m\n", color=:yellow)
        #     if done
        #         break
        #     end
        # end
        
        if i > 100 && occursin(r"^\[SpinUp]", line)
            START = true
        end
        if START
            m, done = process_line(line)
            printstyled("$m\n", color=:yellow)
            if done
                break
            end
        end


    end
end

fls = fls|>unique
map(x->vcat(println("-> ",replace(x, "\\"=> "/")),
    isspin(x;
    match=r"\$time")),fls)