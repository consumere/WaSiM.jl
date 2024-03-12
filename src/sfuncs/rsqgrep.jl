# 
function rsqgrep()
    path = pwd()
    files = glob(r"_output.txt|_outputjl")
    @printf("Searching for R² values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"R2.*[0-9].[3-9]", output)
        if !isempty(match)
            fn = first(split(file, "_qout"))
            for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                line = strip(line)
                line = join(split(line), " ")
                printstyled(rpad("$fn:", 30), lpad("$line\n", 10), color = :green)
            end
        end
    end
end
#rsqgrep()
