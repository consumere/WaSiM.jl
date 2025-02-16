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
 
    ##very nice meta programming... 
    macro vv(s) vgjl(s);end
    #@vv "unter"
    macro vpy(s) vgpy(s);end
    #@vpy "climate"
    macro vr(s) vgr(s);end
    #@vr "climate"
    macro vct(s) vgctl(s);end
    #@vct "das ist"
    macro rg(s) rglob(s);end
    macro glb(s) glob(s);end
    macro gl(s) glob(s)|>first;end

    macro ncrm() ncrem="C:/Users/Public/Documents/Python_Scripts/julia/ncremover.jl";include(ncrem);end
    macro rasrm() remer="C:/Users/Public/Documents/Python_Scripts/julia/raster_remover.jl";include(remer);end
    macro nco(s) nconly(s);end
    
    macro wajs() pt="C:\\Users\\Public\\Documents\\Python_Scripts\\julia\\wajs.jl";include(pt);end

    