# 
function ctlx()
    matches::Vector{Any} = []
    for file in readdir(".")
        if endswith(file, ".xml")
            path = joinpath(pwd(), file)
            open(path) do f
                for line in eachline(f)
                    if occursin("compiling symbols in control file ", line)
                        fields = split(line)[8:end]
                        println(join(fields, " "))
                        out = join(fields, " ")
                        push!(matches, out)
                    end
                    if occursin("looking for starting date in ", line)
                        fields = split(line)[9:end]
                        str = replace(join(fields, " "),r"\"/>.*" => " -> ")
                        printstyled("discharge file is: $str\n",color=:light_green)
                    end
                end
            end
        end
    end
    if !isempty(matches)
        fl = first(matches)
        fl = split(fl) |> last
        fl = split(fl, "\"") |> first
        if (!occursin("regio",fl) && occursin("regio",pwd()))
            fl = replace(fl,"control"=>"D:/Wasim/regio/control")
        elseif (!occursin("brend",fl) && occursin("brend",pwd()))
            fl = replace(fl,"control"=>"D:/Wasim/Tanalys/DEM/brend_fab/control")            
        elseif (!occursin("temp",fl) && occursin("saale",pwd()))
            fl = replace(fl,"control"=>"D:/temp/saale/control")
        end
        return string(fl)
    else
        @warn "no control file or xml found!..."
        return ""
    end
end

