#pt=/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/vg.jl
#julia --startup-file=no $pt "kol2" ctl
#julia --startup-file=no $pt "alpha = 0.3838 0.3838 1.1319 1.1319 1.1319 ;Par_n = 1.500993 " v4.ctl

"""
recursive grep
"""
function vgrec(snippet::AbstractString, file_ending::AbstractString)

    owd=pwd()
    file_ending = lowercase(file_ending)
    
    
    
    for (root, dirs, files) in walkdir(owd)
        for (index, file) in enumerate(files)
            if (endswith(file, file_ending))           
                #printstyled("check $index: $file\n",color=:light_blue)
                printstyled("$index...",color=:light_blue)
                pt=(joinpath(root, file))
                open(pt) do f
                    counter = 0 # Zähler initialisieren
                    for line in eachline(f)
                        counter += 1 # Zähler erhöhen
                        #if contains(line,snippet)
                        if occursin(Regex(snippet,"i"),line)
                            printstyled("\n$counter:\t",color=:light_red) 
                            printstyled("$pt:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                            printstyled("$line\n",color=:green,bold=true) 
                        end
                    end
                end
            end
        end
    end
end

if length(ARGS) < 2
    println("Usage: vg <pattern> <file_ending>")
    exit(0) #indicating that the program completed successfully.
end

#vg(ARGS...)
vgrec(ARGS...)


########junk###########
"""
gpt grep
"""
function vg(pattern::AbstractString, file_ending::AbstractString)
    if length(ARGS) < 2
        println("Usage: vg <pattern> <file_ending>")
        return
    end

    regex = Regex(pattern)
    file_ending = lowercase(file_ending)

    current_dir = pwd()
    walker = readdir(current_dir)
    #ck = 0

    for file in walker
        if isfile(file) && endswith(lowercase(file), ".$file_ending")
            #give a check counter
            # ck =+ 1 
            # printstyled("check $ck: $file\n",color=:light_blue)
            open(file) do f
                for (line_number, line) in enumerate(eachline(f))
                    if match(regex, line) !== nothing
                        #printstyled(stdout, ("$(abspath(file)):$(line_number): $line \n"), color=:light_green)
                        printstyled(stdout, ("$(file):$(line_number): $line \n"), color=:light_green)
                        #printstyled("$(abspath(file)):$(line_number): $line ")
                        #printstyled(stdout*"\n", color=:light_green)
                        
                    end
                end
            end
        end
    end
end