# 
function process_line(line::AbstractString)
            m = strip(line)
            if occursin(match, line)
                return m, true
            end
            return m, false
        end

        # if isnothing(filename)
        #     filename = ctl3()
        # end
    
        for (i, line) in enumerate(eachline(filename))
            if i > 100 && occursin(r"^\[routing_model]", line)
                START = true
            end
            if START
                m, done = process_line(line)
                printstyled("$m\n", color=:green)
                if done
                    break
                end
            end
        end
    end

    """
    lpro projects from 25832 to 4326
    """
    