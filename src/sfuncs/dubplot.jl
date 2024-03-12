# 
function dubplot(pt::Union{String,DataFrame})
        
        if pt isa String
            df = waread2(pt)
        else
            df = pt
        end
               
        
        #v = nonunique(df,:date)
        yl = try
            collect(DataFrames.metadata(df))[1][2]|>dirname|>splitpath|>last    
        catch
            @info "No basename in metadata!"
            raw""
        end
        
        
        dubs = df[findall(x -> count(==(x), df.date) > 1, df.date),:]
        grp = groupby(dubs,:date)
        v = [size(yx,1) for yx in grp]|>unique
        println("unique counts: $(v)")
        if length(v) == 0
            @info("Warning: no unique count!")
            return
        else
            println("unique counts: $(v)")
        end
        
        StatsPlots.plot()
        
        # i = 0;
        for group in grp
            # cntz = size(group,1)
            # i =+ 1
            @df group StatsPlots.violin!(
                cols(propertynames(dubs)[1]), 
                # annotations = Plots.text("no: $cntz"),
                # annotations = (i,
                # maximum(group[!,1]),
                # "no: $cntz", :top),
                ylabel=yl, 
                xlabel="Group", 
                legend=false)
        end
        for i in 1:size(grp, 1)
            Plots.annotate!(i,
            maximum(grp[i][!,1]),    
            (size(grp[i],1),9,:center,:top,:black))
        end
        title!("Bandwith of duplicate values")
    end

    """
    prints pretty_table of df
    and copies to clipboard
    kwars are passed to pretty_table
    tblcb(kd;backend = Val(:latex))
    """
    