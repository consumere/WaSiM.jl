# 
function corrbar(a::Vector{Float64}, b::Vector{Float64})
        
        df_A = a
        df_B = b
        #ti = "$a vs $b"
        ti = raw""

        #cor(df_A, df_B)^2
    
        #Calculate correlations and replace NaN with 0
        correlations = Vector{Float64}(undef, size(df_A,1))
        for i in 1:size(df_A, 1)
            correlations[i] = cor(df_A, df_B)^2
        end
        replace!(correlations, NaN => 0)
    
        p0 = Plots.bar(1:size(df_A, 1), correlations,
            legend = false,
            title = ti,
            fillcolor = ifelse.(correlations .> 0.35, 
                "cornflowerblue", "coral2"),
            xticks = (1:size(df_A, 1), propertynames(df_A)),
            xrotation = 45,
            xlabel = "",
            ylabel = "Correlation R²",
            left_margin = 10mm,
            bottom_margin = 2mm);
    
        ann = map(x->string.(round(x;sigdigits=2)),correlations)
    
        for i in 1:size(df_A, 1)
            Plots.annotate!(i,correlations[i],
            (ann[i],9,:center,:top,:black))
            # println("R² "*ann[i]*" of Basin 
            # "*(df_A)[i]*" added")
        end
    
       return p0
    end
    
    
    """
    selects first dfcol (for so)
    dt annotations
    cmplot(;temp=r"^so_temper",prec=r"^so_prec")
    see also:
    climateplot(r"^tem",r"^pre")
    col = subbasin of interest
    climateplot(r"^temp",r"pre";col="tot_average")
    ws. prc and temp tauschen un opacity einstellen....
    #