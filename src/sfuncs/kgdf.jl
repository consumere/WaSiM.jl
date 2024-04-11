# 
function kgdf(ds::Vector{Any})
        
    #     ds.name=map(x->replace(x,r"-qoutjl.*" =>"","_" => " "),ds.name)
    #     ann = map(x->string.(round(x;sigdigits=2)),ds.KGE)
    #     p1 = Plots.bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
    #     title = splitpath(pwd())|>last, 
    #     xrotation = 45, 
    #     fmt = :png, 
    #     size = (800, 600), 
    #     fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
    #     #annotations = (ds.name,ds.KGE, ann, :top),
    #     xaxis = "",
    #     left_margin = 10mm,
    #     bottom_margin = 15mm, 
    #     bar_width = 0.6);

    #     for i in 1:length(ds.name)
    #         Plots.annotate!(ds.name[i],ds.KGE[i],(ann[i],11,
    #             :center,:top,:black))
    #         #println(ann[i]*" added")
    #     end
    #     return p1
    # end


    