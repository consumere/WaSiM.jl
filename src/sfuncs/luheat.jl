# 
function luheat(dfs::Vector=dfs,x::String="LAI";colors=:grays)
        month_abbr = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
        ti = findindf(dfs[1],Regex(x,"i"))[1,1]|>first|>strip
        #ti = findindf(dfs[1],Regex(x,"i")).Parameter|>first|>strip
        a = map(k->findindf(k,Regex(x,"i")),dfs)
        a = vcat(a...)|>DataFrame
        z = vcat(map(x->split(replace(
            string(x[1, 1]),
            r"\s+"=>" ",
            r"{"=>"",
            r"method"=>"")),
            dfs)...)
        a.nm = filter(x->length(x)>2,z)
        a = a[:, sort(names(a),rev=true)] #sort columns
        a.nm .= replace.(a.nm, r"_" => " ")
        begin 
            value_str = string.(a.Value)
            vs = split.(value_str)
            z0 = [parse.(Float64, z) for z in vs]
            m = hcat(z0...)'
            #myc = colormap("Grays",10;logscale=true)
            #myc = colormap("Grays",Int(round(findmax(m)[1]));logscale=false)
            #heatmap(m,c=:grays)
            heatmap(m,c=colors)
            #heatmap(m,c=myc[5:end-1])
            xticks!(1:12,month_abbr)
            yticks!(1:size(m,1),a.nm)
            
            title!("$ti")

            #add annotation for each field
            qv = quantile(m)[4] #|>round
            for i in 1:size(m,1)
                for j in 1:size(m,2)
                    #tcolor = m[i,j] < qv ? :black : :white
                    tcolor = m[i,j] < qv ? :white : :black
                    annotate!(j,i,text(string(
                        round(m[i, j]; digits=2)
                        ),7,tcolor,rotation=-45))
                end
            end
            return plot!()
        end
    end
    
    

    """
    runs perl wsl \$bsp x x 4 5 or 
    \$f1 \$f2 \$c1 \$c2 ∇
    pers(x::Union{String,Regex}=r"qoutjl";c1=4,c2=5,f2=nothing,tofile=false)
    """
    