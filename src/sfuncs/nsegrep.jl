# 
function nsegrep()
    path = pwd()
    files = glob(r"_output.txt|_outputjl") #non-recurse
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"mNSE", output)
        if !isempty(match)
            fn = first(split(file, "_qout"))
            for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:", 45), lpad("$line\n", 10), color = :green)
            end
        end
    end
end

