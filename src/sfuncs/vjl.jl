# 
function vjl(regex)
    # greps jl from julia folder
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia";
    file_ending=".jl"
    files = filter(file -> endswith(file, file_ending), readdir(pt,join=true))
    for file in files
        open(file) do f
            counter = 0
            for line in eachline(f)
                counter += 1
                if occursin(Regex(regex,"i"), line)
                    println("$file: $counter:\t $line")
                end
            end
        end
    end
end

