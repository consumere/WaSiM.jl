# 
function vgctl(snippet::AbstractString)
        """
        hint. vgctl("set \$TS")
        """
        owd=pwd()
        nwd="D:/Wasim/regio/control/"
        nwd2="D:/temp/saale/control/"
        nwd3="D:/Wasim/Tanalys/DEM/brend_fab/control/"
        nwd4="D:/Wasim/regio/control/"
        nwd5="D:/Wasim/streu/control/"
        cd(nwd)
        #println("greps from *ctl from  \n$nwd and \n$nwd2...")
        println("greps from *ctl from  \n$nwd \n$nwd2 \n$nwd3 \n$nwd4 \n$nwd5...")
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd2)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd2",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd3)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd3",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd4)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd4",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(nwd5)
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter: $nwd5",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        end
        cd(owd)
    end

    