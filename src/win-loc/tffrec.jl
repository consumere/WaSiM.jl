if length(ARGS) == 0
	println("need args! <file>...")
    exit()
end
using Printf

function tff3(x::Vector{String})
    for filename in x
        if (
            (!occursin(r"html|sh|txt|xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",
                filename))
            )
        
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
                        #println("filename: $filename")
                        printstyled("filename: $filename\n",color=:yellow)
                        ncols = length(fields)
                        datcols = ncols - 4
                        println("no of fields: $datcols")
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

            if ncols <= 4
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols == 5
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            elseif ncols >= 6
                for key in sort(collect(keys(b)))
                    println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", m[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", m[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
                end
            end
        end
    end
end

fullpath = pwd()
nms::Vector{String} = []
#ARGS="qg"

#xm = parse.(string,ARGS[1])
xm = ARGS[1]
println("looking for match of *$xm* ...")

for (subdir, _, filenames) in walkdir(fullpath)
    for filename in filenames
        if (occursin(Regex(xm,"i"),filename) && (!occursin(r"sh|yrly|nc|png|svg",filename)) )
        pt = joinpath(subdir, filename)
        #sizes[fullpath] = stat(fullpath).size
        push!(nms,pt)
        end
    end
end
#println(    nms    )
#println(  typeof(nms)    )

#Vector{String}
# # Neuen Vektor initialisieren
# vector_string = Vector{String}(undef, length(nms))

# # Elemente des ursprünglichen Vektors in Strings umwandeln
# for i in 1:length(nms)
#     vector_string[i] = string(nms[i])
# end


ti= try 
    tff3(nms)
catch 
    wd = pwd();
    @error "
    no recursive match for $ARGS on $wd !
    ...exiting now"; 
    exit()
end

#in wsl:
#julia --startup-file=no $(wslpath "C:\Users\Public\Documents\Python_Scripts\julia\tffrec.jl") rg*

