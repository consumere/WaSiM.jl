
"""
filters internal WaSiM stats of routed discharge files
works recursively
"""
function qga(;rootdir=".", prefix="qgk")
    files = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                push!(files, joinpath(looproot, filename)) 
            end
        end
    end
    
    for z in files
        println(raw"file:	",basename(z),"...")
        m = filter(line -> occursin(r"^[LIN. R]|^[LOG. R]|^CO",line), readlines(open(z)))
        str = []
        for l in m
            x = replace(l, r"\s+" => "\t")
            x = replace(x, ".\t" => " ")
            #x = replace(x, "\t0" => "")
            x = replace(x, "R-SQUARE\t" => "R-SQUARE:  ")
            x = replace(x, "LOG\t" => "LOG: ")
            x = replace(x, "LIN\t" => "LIN: ")
            x = replace(x, "\t" => " |")
            x = replace(x, r"( 0 | )" => " ")
            #x = replace(x, r"(\t|0 \| )" => "\t|\t")
            #x = replace(x, ": " => "\n")
            #x = replace(x, " | " => "\n")
            #x = replace(x, "\t+" => " | ")
            println(strip(x))
            #z = split(x, " | ")
            # z = (split(x, "|"))|>collect
            #z = adjoint(z)
            #println(z)
            # #println(lpad.(z, 15))
            # #string.(z)|>permutedims|>println
            #joined_string = join(string.(z), "\n")
            #adjoint(joined_string)|>println
            # push!(str,joined_string)
        end
        #println(str)
        # for l in m
        #     replace!(l, r"\s+" => "\t")
        #     replace!(x, ".\t" => " ")
        # end
        # mt = transpose(m)Union{}[ERROR: LoadError: MethodError: no method matching transpose(::String)      
        # println(mt)
    end
    return nothing
end

qga()

#julia --startup-file=no --optimize=0 --threads 8  "C:\Users\Public\Documents\Python_Scripts\julia\qga.jl"