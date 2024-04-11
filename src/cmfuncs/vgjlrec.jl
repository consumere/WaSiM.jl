# 
function vgjlrec(snippet::AbstractString)
        """
        recursive grep
        """
        owd="C:/Users/Public/Documents/Python_Scripts/julia"
        for (root, dirs, files) in walkdir(owd)
            for file in files 
                if (endswith(file, ".jl"))
                    pt=(joinpath(root, file))
                    open(pt) do f
                        counter = 0 # Zähler initialisieren
                        for line in eachline(f)
                            counter += 1 # Zähler erhöhen
                            if contains(line,snippet)
                                printstyled("$counter:\t",color=:light_red) 
                                printstyled("$pt:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                                printstyled("$line\n",color=:green,bold=true) 
                            end
                        end
                    end
                end
            end
        end
    end

    