# 
function vgpy(snippet::AbstractString)
    #owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "C:/Users/Public/Documents/Python_Scripts"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/c/Users/Public/Documents/Python_Scripts"
    end
    
    files = filter(file -> endswith(file, ".py"), readdir(script_dir,join=true))
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if Base.contains(line,snippet)
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

"""
Usage: vgctl("set \$TS")
"""
