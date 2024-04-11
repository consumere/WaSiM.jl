# 
function jldfnm(x::Vector{String})
        """
        with names
        """
        result = Vector{DataFrame}()

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
                        ncols = length(fields)
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

            df = DataFrame(
                Year = String[],
                B = Float64[],
                M = Float64[],
                H = Float64[],
                B_Mean = Float64[],
                M_Mean = Float64[],
                H_Mean = Float64[],
                Counts = Int64[],
                Name = String[]
            )
            metadata!(df, "filename", filename, style=:note);
            for key in sort(collect(keys(b)))
                push!(df, (
                    Year = key,
                    B = b[key],
                    M = m[key],
                    H = h[key],
                    B_Mean = b[key] / cnte[key],
                    M_Mean = m[key] / cnte[key],
                    H_Mean = h[key] / cnte[key],
                    Counts = cnte[key],
                    Name = filename
                ))
            end
            push!(result, df)
        end

        return result
    end

    