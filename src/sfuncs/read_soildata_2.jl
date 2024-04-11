# 
function read_soildata_2(filename::String)
        #data = readlines(filename)
        # data = open(fl) do io
        #     a = readbetween(io,"{","}")
        #     return(a)
        # end
        data = open(filename) do io
            a = readbetween(io, "soil_table", "special_output")
            return(a)
        end
        
        #output = DelimitedFiles.readdlm(filename,';', String)
        data = Grep.grep(r"^th.*|^[0-9]",data)
        data = broadcast(x -> replace(x,    
                r"^#.*" => "",
                r"[[]].*" => "",
                r"[{].*" => "",
                r"method" => "",
                r"MultipleHorizons" => "",
                r"}" => ""), data)
        bks = broadcast(x -> split(x, " ",limit=3), data)
        nos = broadcast(x -> if length(x)==2 x[1] end, bks)
        nos = filter(x -> !isnothing(x),nos)
        nos = broadcast(x -> parse.(Int, x),nos)
        
        #tck = broadcast(x -> if (length(x)==3  & startswith(x[1],"thickness")) x[3] end, bks)
        tck = broadcast(x -> if (startswith(x[1],"thickness")) x[3] end, bks)
        filter!(x -> !isnothing(x),tck)
        tck = broadcast(x -> split(x),tck)
        
        #foreach(x -> parse.(Int, x[1]),tck)
        
        #collect(tck)
        
        Z=[]
        for m in (tck)
            i = [parse(Float64, x) for x in m]
            push!(Z,i)
        end
        
        #[cumsum(arr) for arr in Z]
        
        #zd = DataFrame(permutedims(Z),:auto)
        zd = DataFrame((Z),:auto)
        for x in 1:length(nos)
            newname = Symbol.(nos[x])
            rename!(zd,Dict(names(zd)[x]=>newname))
        end
        # col_sums = sum.(eachcol(zd))
        # hcat(col_sums*.1, nos)
       
    
        return zd
    end

    """
     soil types
     returns a DataFrame
     needs rework
    """
    