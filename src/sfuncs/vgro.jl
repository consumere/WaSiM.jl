# 
function vgro(snippet::AbstractString)
    owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "D:/Fernerkundungsdaten/Klassifikation/R-Sessions"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions"
    end
    
    cd(script_dir)

    files = filter(file -> endswith(file, ".R"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if Base.contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(owd)
end

