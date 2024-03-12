# 
function ctlg(dir_path::String, match::String;file_ending="ctl")
        """
        dir_path::String, match::String;file_ending="ctl"
        """
        for file in readdir(dir_path)
            if occursin(file_ending, file)
                fx = joinpath(dir_path, file)
                prevline = ""
                for (i, line) in enumerate(eachline(fx))
                    if findfirst(match, line) !== nothing
                        printstyled("File: $file\n",color=:red)
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

    macro sdjl() fx="C:/Users/Public/Documents/Python_Scripts/julia/sd2.jl";include(fx);end
    #@sdjl
    macro sf(s) glob(s);end
    macro listdir() ls();end

    