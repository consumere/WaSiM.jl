# 
function tff2(x::Vector{String})
        for filename in x
            b = Dict{String, Float64}()
            m = Dict{String, Float64}()
            h = Dict{String, Float64}()
            cnte = Dict{String, Int64}()
            ncols = 0

            open(filename) do file
                first_line = true
                

                for (i, line) in enumerate(eachline(file))
                    fields = split(line)
                    if first_line
                        println("filename: $filename")
                        ncols = length(fields)
                        println("no of fields: $ncols")
                        println("year\t$(join(fields[5:end], "\t"))")
                        first_line = false
                        continue
                    end

                    if match(r"^\d{4}$", fields[1]) != nothing
                        cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
                        b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
                        m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
                        h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
                    end
                end
            end

            if ncols <= 5
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols == 6
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols >= 7
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", m[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", m[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            end
        end
    end

    