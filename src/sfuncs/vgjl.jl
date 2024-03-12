# 
function vgjl(snippet::AbstractString)
    owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "C:/Users/Public/Documents/Python_Scripts/julia"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
    end
    
    cd(script_dir)
    
    # files = []
    files = filter(file -> endswith(file, ".jl"), readdir())
    fwin = filter(file -> endswith(file, ".jl"), 
        readdir(script_dir*"/win", join=true))
    files = vcat(files,fwin)
    for file in files
        open(file) do f
            counter = 0 # Initialize the counter
            for line in eachline(f)
                counter += 1 # Increment the counter
                if Base.contains(line, snippet)
                    printstyled("$counter:\t", color=:light_red)
                    printstyled("$file:\t", color=:light_magenta, underline=true, blink=false, bold=true)
                    printstyled("$line\n", color=:green, bold=true)
                end
            end
        end
    end
    
    cd(owd)
end

