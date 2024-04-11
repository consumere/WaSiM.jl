# 
function nsx(dfi::DataFrame)
        """
        nse barplot with all values
        """
        dfi.name=map(x->replace(x,r"-qoutjl.*" => ""),dfi.name)
        ann = map(x->string.(round(x;sigdigits=3)),dfi.NSE)
        Plots.bar(dfi.name, dfi.NSE, xlabel = "Name", ylabel = "NSE", legend = false, 
            title = "x", xrotation = -25, fmt = :png, size = (800, 600), 
            fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
            annotations = (dfi.name,dfi.NSE, ann, :top),
            xtickfont = font(6),    
            bar_width = 0.77)        
    end

    macro sdjl() fx=src_path*"/sd2.jl";include(fx);end
    #@sdjl
    macro sf(s) glob(s);end
    macro listdir() ls();end

    