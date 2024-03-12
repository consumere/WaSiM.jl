# 
function vgctl(snippet::AbstractString)
    owd = pwd()
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform

    if platform == "windows"
        nwd = "D:/Wasim/regio/control/"
        nwd2 = "D:/temp/saale/control/"
        nwd3 = "D:/Wasim/Tanalys/DEM/brend_fab/control/"
        nwd4 = "D:/Wasim/Goldbach/revision/control/"
        nwd5 = "D:/Wasim/streu/control/"
    else
        # Modify these paths for your WSL setup
        nwd = "/mnt/d/Wasim/regio/control/"
        nwd2 = "/mnt/d/temp/saale/control/"
        nwd3 = "/mnt/d/Wasim/Tanalys/DEM/brend_fab/control/"
        nwd4 = "/mnt/d/Wasim/Goldbach/revision/control/"
        nwd5 = "/mnt/d/Wasim/streu/control/"
    end

    paths = [nwd, nwd2, nwd3, nwd4, nwd5]

    for path in paths
        cd(path)
        println("Searching in directory: $path")
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0
                for line in eachline(f)
                    counter += 1
                    if Base.contains(line, snippet)
                        printstyled("$counter: $path", color=:light_red)
                        printstyled("$file:\t", color=:light_magenta, underline=true, blink=false, bold=true)
                        printstyled("$line\n", color=:green, bold=true)
                    end
                end
            end
        end
    end

    cd(owd)
end

