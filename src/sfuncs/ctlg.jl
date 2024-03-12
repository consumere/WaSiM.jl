# 
function ctlg(dir_path::String, match::String;file_ending="ctl")
    for file in readdir(dir_path)
        #if ( occursin(file_ending, file) &&
        if ( endswith(file,file_ending) &&
            !occursin(r"tar$|tex$|pl$|sh$|csv|html|xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg", file))
            fx = joinpath(dir_path, file)
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

