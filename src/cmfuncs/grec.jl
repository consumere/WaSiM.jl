# 
function grec(dir_path::String, file_ending::String, match::String)
        """
        grep recursive and looks for matches in files
        """
        for (looproot, dirs, filenames) in walkdir(dir_path)
            for file in filenames
                fx = joinpath(looproot, file)
                if isfile(fx) && occursin(Regex(file_ending*"\$","i"), fx)
                    #printstyled("check $file\n",color=:yellow)
                    #fx = joinpath(looproot, file)
                    prevline = ""
                    for (i, line) in enumerate(eachline(fx))
                        if findfirst(match, line) !== nothing
                            printstyled("File: $fx\n",color=:red)
                            println("$(prevline)")
                            printstyled("$line\n",color=:green)
                            nextline = readline(fx)
                            println("$(nextline)")
                        end
                        prevline = line
                    end
                end
            end
        end
    end

    