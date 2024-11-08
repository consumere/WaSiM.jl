# 
function readbetween(io::IO, start::Regex, stop::Regex)
    output = Vector{String}()
    while !eof(io)
        line = readline(io)
        if occursin(start,line)
            push!(output, line)
            while !eof(io)
                line = readline(io)
                if occursin(stop,line)
                    skip(io, length(line))
                    break
                end
                push!(output, line)
            end
            break
        end
    end
    return output
end