# 
function vgr(regex, file_ending)
    rootdir=pwd()
    println("starting on: $rootdir...\n searching for >> $regex << with file ending >> $file_ending <<\n")
    files = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
            #if (occursin(Regex(prefix,"i"),filename))
            if (endswith(filename, file_ending))
                push!(files, joinpath(looproot, filename)) 
            end
        end
    end
    #files = filter(file -> endswith(file, file_ending), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if occursin(Regex(regex,"i"), line)
                    println("$file: $counter:\t $line")
                end
            end
        end
    end
end

