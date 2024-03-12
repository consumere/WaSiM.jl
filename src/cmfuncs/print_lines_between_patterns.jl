# 
function print_lines_between_patterns(filename::AbstractString, start_pattern::AbstractString, end_pattern::AbstractString)
        in_range = false
        
        for line in eachline(filename)
            if occursin(start_pattern, line)
                in_range = true
            end
            
            if in_range
                println(line)
            end
            
            if occursin(end_pattern, line)
                in_range = false
            end
        end
    end

    