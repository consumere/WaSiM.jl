# 
function nsegrepr()
        path = pwd()
        files = rglob(r"_output.txt|_outputjl") #recurse
        @printf("Searching for NSE values > 0.3 in files matching pattern %s\n", path)
        for file in filter(file -> endswith(file, "_output.txt"), files)
            output = DelimitedFiles.readdlm(file,'\t', String)
            match = Grep.grep(r"NSE.*[0-9].[3-9]",output)
            if !isempty(match)
                fn = first(split(file,"_qout"))
                for line in match
                    line = strip(line) 
                    line = join(split(line), " ")  ##remove inner whitespaces
                    printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
                end
            end
        end
    end

    """
    --- reader with fewer constrains ---
    no |> dropmissing 
    df[!,Cols(r"^Col|date")] |>dfp  
    """
    