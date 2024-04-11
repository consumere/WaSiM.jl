# 
function heat(x::Union{String,Regex,DataFrame};ann=true,kw...)
        if x isa String
            printstyled("reading $x\n",color=:light_red)
            df = waread2(x;silencewarnings=false)
            dropmissing!(df)    
        elseif x isa Regex
            x = first(dfonly(x))
            printstyled("reading $x\n",color=:light_red)
            df = waread2(x;silencewarnings=false)
            dropmissing!(df)
        else
            df = copy(x)
            dropmissing!(df)
        end

        if ncol(df)<2
            @error "only one column available!"
            return
        end

        if ncol(df)>20
            @warn "more than 20 columns, heatmap will be unreadable!"
            return
        end
        #md = cor(Matrix(df[:, Not(:date)])).^2
        #md = Matrix(select(df, Not(Cols(r"date|year|month|day"))))
        md = cor(Matrix(select(df, Not(Cols(r"date|year|month|day"))))).^2
        md = convert(Matrix{Union{Float64, Missing}}, md)
        #md = cor(md)^2
        replace!(md, 1.0 => missing)
        replace!(md, 0.0 => missing)

        p = heatmap(md, 
            #c=:balance, 
            #c=:lightrainbow, grep(r"gr",cs)
            #c=:grays,
            #c = Plots.colormap("RdBu",nrow(df);logscale=true),
            c = Plots.colormap("RdBu",size(md, 1),mid=0.2),
            linewidth=0.8, 
            cbar=true,
            colorbar_title = "RSQ",
            margins=5mm,
            xticks =false,
            yticks =false,
            grid = :false,
            # legend = :outertopright,
            # legendtext = string(names(df[:, Not(:date)])),#xlabel="", #ylabel="", 
            title="Correlation Heatmap";kw...)
        if ann
            for i in 1:size(md, 1)
                for j in 1:size(md, 2)
                    value = round(md[i, j], digits=2)
                    color = ismissing(value) ? :transparent : :black
                    #color = value==1.0 ? :transparent : :black
                    Plots.annotate!(j, i, 
                        Plots.text(string(value), 9, color, :center, 
                        halign=:center, rotation=-35.0))
                    # annotate!(j - 0.5, i,
                    #     text(value, 
                    #         8, color, :center))
                end
            end
        end

        nm = names(df[!,Not(Cols(r"date|year|month|day"))])
        nm = replace.(nm,
            r"_" => " ",
            r"qoutjl|qout" => "",
            #"-qoutjl"=>"",
            r"#" => "")
        Plots.xticks!(1:length(nm), nm) #reverse(nm)
        Plots.yticks!(1:length(nm), nm;rotation=-35)

        return(p)
    end

    """
    returns the position of each element in a vector
    """
    