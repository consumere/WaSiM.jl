# 
function grep_with_context(pattern, filename, context)
        lines = readlines(filename)   
        for (i, line) in enumerate(lines)
            if occursin(pattern, line)
                    #println("$(filename):$(i): $(line)")
                    for j in max(1, i - context):min(length(lines), i + context)
                        #printstyled("$(filename):$j: $(lines[j])\n",color=:red)
                        printstyled("$(filename):$j:",color=:red)
                        printstyled("$(lines[j])\n",color=:green)
                    end
                    println("-" ^ 80)
            end
        end
    end
        
    #grep_with_context("routing_model", file_paths[1], 2)

    """
    pts = readdir(dirname(infile);join=true)
    filter!(x->occursin(r"ctl",x),pts)
    grep_files("routing_model", pts, 2)
    """
    