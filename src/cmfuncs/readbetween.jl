# 
function readbetween(io::IO, start::String, stop::String)
        output = Vector{String}()
        while !eof(io)
            line = readline(io)
            if contains(line, start)
                push!(output, line)
                while !eof(io)
                    line = readline(io)
                    if contains(line, stop)
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


    