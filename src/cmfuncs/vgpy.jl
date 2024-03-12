# 
function vgpy(snippet::AbstractString)
        owd="C:/Users/Public/Documents/Python_Scripts"
        files = filter(file -> endswith(file, ".py"), readdir(owd,join=true))
        for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled(rpad("$file:",50),color=:light_magenta)
                        #printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled(lpad("$line\n",30),color=:green,bold=true)
                        #underline = true 
                    end
                end
            end
        end
    end

    