# 
function klog(df::Union{String,DataFrame};lin=false)
        if lin 
            cols=Cols(2,4)
        else
            cols=Cols(3,5)
        end

        if df isa String
            ofl = df #should be "route.txt"
            df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
            rename!(df,1=>"sim",2=>"obs",3=>"name")
            df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
            df.name=map(x->replace(x,r"_>.*" => ""),df.name)
        end
        #dsu = qbb()|>last #recursive
        dsu = qba()
        M = parse.(Float64,dsu[!,cols])|>y->subset(y,1=> ByRow(>(0.)))|>Matrix
        mns = names(dsu[!,cols])
        B = kde(M)
        Plots.plot(B,legend=false,xlabel=mns[1],ylabel=mns[2])
        dsu.basin  .= df.name
        name_mapping = df    
        msk = @rsubset(parse.(Float64,dsu[!,cols]) .> 0)
        nmn = dsu[msk[!,1],:]
        nmn.basin=map(x->replace(x,r"_" => " "),nmn.basin)
        Plots.annotate!([(M[i,1], M[i,2], 
            Plots.text(
            nmn.basin[i], 
            8, :left, 
            halign=:center, 
            rotation=-15.0)) for i in 1:size(nmn, 1)])
        Plots.plot!(legend=false)
    end


    """
    pycall 