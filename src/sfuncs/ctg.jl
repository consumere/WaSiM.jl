# 
function ctg(match::String ; dir=".", file_ending=".ctl")
    for file in readdir(dir)
        if occursin(file_ending, file)
            fx = joinpath(dir, file)
            prevline = ""
            for (i, line) in enumerate(eachline(fx))
                if findfirst(match, line) !== nothing
                    printstyled("File: $file\n",color=:red)
                    println("$(prevline)")
                    printstyled("$line\n",color=:green)
                    nextline = readline(fx)
                    println("$(nextline)")
                end
                prevline = line
            end
        end
    end
end

