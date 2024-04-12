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


function tff4(x::Vector{String})
    for filename in x
        if !occursin(r"html|sh|txt|xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg", filename)
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
                        printstyled("filename: $filename\n", color=:yellow)
                        ncols = length(fields)
                        datcols = ncols - 4
                        println("no of fields: $datcols")
                        println("year\t$(join(fields[5:end], "\t"))")
                        first_line = false
                        continue
                    end

                    if match(r"^\d{4}$", fields[1]) != nothing
                        cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
                        
                        for i in 5:length(fields)
                            field_value = parse(Float64, replace(fields[i], r"-9999" => "0"))
                            if i == 5
                                b[fields[1]] = get(b, fields[1], 0.0) + field_value
                            elseif i == length(fields) - 1
                                m[fields[1]] = get(m, fields[1], 0.0) + field_value
                            elseif i == length(fields)
                                h[fields[1]] = get(h, fields[1], 0.0) + field_value
                            end
                        end
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


ti= try 
    #tff3(ARGS)
    tff4(ARGS)
catch 
    wd = pwd();
    @error "
    no match for $ARGS on $wd !
    this is non-recursive... exiting now"; 
    exit()
end

#in wsl:
#julia --startup-file=no $(wslpath "C:\Users\Public\Documents\Python_Scripts\julia\tff.jl") rg*


# function tff2(x::Vector{String})
#     for filename in x
#         b = Dict{String, Float64}()
#         m = Dict{String, Float64}()
#         h = Dict{String, Float64}()
#         cnte = Dict{String, Int64}()
#         ncols = 0

#         open(filename) do file
#             first_line = true
            

#             for (i, line) in enumerate(eachline(file))
#                 fields = split(line)
#                 if first_line
#                     #println("filename: $filename")
#                     printstyled("filename: $filename\n",color=:yellow)
#                     ncols = length(fields)
#                     println("no of fields: $ncols")
#                     println("year\t$(join(fields[5:end], "\t"))")
#                     first_line = false
#                     continue
#                 end

#                 if match(r"^\d{4}$", fields[1]) != nothing
#                     cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
#                     b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
#                     m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
#                     h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
#                 end
#             end
#         end

#         if ncols <= 5
#             for key in sort(collect(keys(b)))
#                 println("$key\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
#             end
#         elseif ncols == 6
#             for key in sort(collect(keys(b)))
#                 println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
#             end
#         elseif ncols >= 7
#             for key in sort(collect(keys(b)))
#                 println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", m[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", m[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
#             end
#         end
#     end
# end

# function rglob(prefix::AbstractString)
#     rootdir="."
#     results::Vector{String} = []
#     for (looproot, dirs, filenames) in walkdir(rootdir)
#         for filename in filenames
#             if (occursin(Regex(prefix,"i"),filename)
#                 && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",
#                     filename))
#                 )
#                 push!(results, joinpath(looproot, filename)) 
#             end
#         end
#     end
#     return results
# end

#printstyled("with txt grep!\n",color=:green)
# printstyled(typeof(ARGS[1]),"\t",
# ARGS,"!\n",color=:green)

#prefix=ARGS[1]

#filenames = ARGS

#vec = rglob(parse.(string,ARGS))
#vec = rglob(string(ARGS))
#vec = rglob(string(ARGS[1]))
# vec = rglob(ARGS[1])
# tff2(vec)


# ts::Vector{String} = []
# path=pwd()
# for file in ARGS
#     #if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",file))
#     if isfile(file) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",file))
#        file_path = joinpath(path, file)
#        printstyled("eval on ",file_path," ...\n",color=:green)
#       #p1 = readf(file_path) #new
#       push!(ts, file)
#     end
# end

# #geht auch.
# #ts = filter(x -> !occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg",x), ARGS)

# tff2(ts)
