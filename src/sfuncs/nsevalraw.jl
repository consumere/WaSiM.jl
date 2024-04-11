# 
function nsevalraw()
        ds = kge_df3()
        ds.name=map(x->replace(x,r"-qoutjl.*" =>"","_" => " "),ds.name)
        ann = map(x->string.(round(x;sigdigits=2)),ds.NSE)
        p1 = Plots.bar(ds.name, ds.NSE, 
            xlabel = "Name", ylabel = "NSE", legend = false, 
            title = splitpath(pwd())|>last, 
            xrotation = 45,
            ylims = (extrema(ds.NSE) .+ [-1.0, 0.5]), 
            #fmt = :png, 
            size = (800, 600), 
            fillcolor = ifelse.(ds.NSE .> 0, "cornflowerblue", "coral2"),
            xaxis = "",
            left_margin = 10mm,
            bottom_margin = 15mm, 
            bar_width = 0.6);

        for i in 1:length(ds.name)
            Plots.annotate!(ds.name[i],ds.NSE[i],
            (ann[i],11,
                :center,:top,:black))
            
        end
        return p1
        
        
        #ann = map(x->string.(round(x;sigdigits=2)),dfi.NSE)
        #fontfamily="Computer Modern",
        # dfi = ds
        # Plots.bar(
        #     dfi.name, 
        #     dfi.NSE, 
        #     xlabel = "Name", ylabel = "NSE", 
        #     legend = false, 
        #     title = splitpath(pwd())|>last, xrotation = 45, 
        #     #fmt = :png, size = (800, 600), 
        #     fillcolor = ifelse.(dfi.NSE .> 0, "cornflowerblue", "coral2"),
        #     annotations = (dfi.name,dfi.NSE, ann, :top),
        #     bar_width = 0.6)        
    end
   

    """
    calculates water balance externally
    using Images
    img=load("waba-jl.png")
    """
    