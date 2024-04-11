# 
function nctodfo(x::AbstractString)
        nc = NCDatasets.NCDataset(x);
        #nc.attrib
        dict = nc|>Dict   
        mykeys = keys(dict)
        #println(string.(mykeys))
        v = filter(x->!occursin(r"time|lon|lat|x|y",x),string.(mykeys))|>first
        time = nc["time"][:]
        datetime_vector = coalesce.(time, missing)
        #df = hcat(nc[v][end,end,:],datetime_vector)
        xm = nc[v]|>size|>first
        xm = Int(round(median(1:xm);digits=0))
        ym = nc[v]|>size|>second
        ym = Int(round(median(1:ym);digits=0))
        df = DataFrame(
                v=>nc[v][xm,ym,:],      #x, y, indices
                "date"=>datetime_vector)
        # df = DataFrame(
        #         v=>nc[v][end,end,:],      #x, y, indices
        #         "date"=>datetime_vector)
        DataFrames.metadata!(df, "filename", x, style=:note);        
        #df.date = Date.(string.(df.x2),"yyyy-mm-ddTHH:MM:SS") #not needed
            # plot(time, nc[v][end,end,:],
                # label=nc[v].attrib["units"],
        # title=nc[v].attrib["long_name"])        
    end

    