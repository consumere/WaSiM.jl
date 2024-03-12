# 
function dsbar(ds::DataFrame)
        ds.name=map(x->replace(x,r"-qoutjl*" => ""),ds.name)
        ann = map(x->string.(round(x;sigdigits=1)),ds.KGE)
        bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
            title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
            fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
            annotations = (ds.name,ds.KGE, ann, :top),
            bar_width = 0.6)
    end

    