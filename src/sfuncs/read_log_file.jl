# 
function read_log_file(log_file::AbstractString)
        timestamps = []
        durations = []
    
        open(log_file) do io
            while !eof(io)
                line = readline(io)
                if startswith(line, "XML File: ")
                    filename = split(line, "XML File: ")[2]
                    push!(timestamps, filename)
                elseif startswith(line, "Duration: ")
                    duration_str = split(line, "Duration: ")[2]
                    duration_parts = split(duration_str, " hours and ")
                    hours = parse(Int, duration_parts[1])
                    minutes = parse(Int, split(duration_parts[2], " minutes")[1])
                    push!(durations, hours * 60 + minutes)
                end
            end
        end
    
        return timestamps, durations
    end
    
    