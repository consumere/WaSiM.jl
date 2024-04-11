# 
function vgrepl(snippet::AbstractString)
        """
        greps from repl_history
        """
        #file = raw"/home/ubu/.julia/logs/repl_history.jl"
        file = raw"C:\Users\chs72fw\.julia\logs\repl_history.jl"
        #files = filter(file -> endswith(file, ".py"), readdir())
        #for file in files
            open(file) do f
                counter = 0 # Zähler initialisieren
                for line in eachline(f)
                    counter += 1 # Zähler erhöhen
                    if contains(line,snippet)
                        printstyled("$counter:\t",color=:light_red) 
                        printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                        printstyled("$line\n",color=:green,bold=true) 
                    end
                end
            end
        #end
        #cd(owd)
    end

    
    
    """
    date to last/first position
    df = select(df, Not(:date), :)
    """
    