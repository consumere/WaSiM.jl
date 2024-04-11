# 
function nseval()
        """
        nse barplot with values > 0
        """
        ds = kge_df3()
        ds.name=map(x->replace(x,r"-qoutjl.*" =>"","_" => " "),ds.name)
        dfi = filter(row -> row.NSE .> 0, ds)
        ann = map(x->string.(round(x;sigdigits=2)),dfi.NSE)
        Plots.bar(dfi.name, dfi.NSE, 
        xlabel = "Name", ylabel = "NSE", legend = false,
        title = splitpath(pwd())|>last, xrotation = 45, 
        #fmt = :png, size = (800, 600),
        fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
        annotations = (dfi.name,dfi.NSE, ann, :top),
        left_margin = 10mm,
        bottom_margin = 15mm, 
        bar_width = 0.6)        
    end

    """
    nse barplot with all values
    """
    