# 
function kgegrepr()
    path = pwd()
    files = rglob(r"_output.txt|_outputjl") #recurse
    @printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file,'\t', String)
        match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
        if !isempty(match)
            fn = first(split(file,"_qout"))
            for line in match
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
            end
        end
    end
end

