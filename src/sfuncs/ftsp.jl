# 
function ftsp(x::AbstractString)
        nc = NCDatasets.NCDataset(x);
        #nc.attrib
        dict = nc|>Dict   
        mykeys = keys(dict)
        println(string.(mykeys))
        time = nc["time"][:]
        v = filter(x->!occursin(r"time|lon|lat|x|y|spatial_ref",x),string.(mykeys))|>first
        xm = nc[v]|>size|>first
        xm = Int(round(median(1:xm);digits=0))
        ym = nc[v]|>size|>second
        ym = Int(round(median(1:ym);digits=0))
        plot(time, nc[v][xm,ym,:],
        label=nc[v].attrib["units"],
        title=nc[v].attrib["long_name"])        
    end

    