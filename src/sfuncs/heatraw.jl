# 
function heatraw(x::Union{String,Regex,DataFrame};selcol::Union{String,Regex}=r"dr",ann=true,kw...)
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

        df = select(df, Cols(selcol))

        if ncol(df)<2
            @error "only one column available!"
            return
        end

        if ncol(df)>30
            @warn "more than 30 columns, heatmap will be unreadable!"
            println(names(df))
            return
        end

        md = cor(Matrix(df)).^2
        md = convert(Matrix{Union{Float64, Missing}}, md)
        replace!(md, 1.0 => missing)
        replace!(md, 0.0 => missing)
        p1 = heatmap(md, 
            c = Plots.colormap("RdBu",size(md, 1),mid=0.33),
            linewidth = 0.9, 
            cbar = true,
            # cbar=true,
            #colorbar_title = "RSQ",
            #colorbar_scale = :log10,
            #colorbar_title_location	= :top,
            #@show plotattr(:Series)
            #https://docs.juliaplots.org/latest/generated/attributes_subplot/
            #colorbar_titlepadding=5,
            colorbar_titlefontsize=12,
            #colorbar_titlefontvalign=:bottom,
            margins=6mm,
            xticks =false,
            yticks =false,
            grid = :false; kw...);
            # # legend = :outertopright,
            # # legendtext = string(names(df[:, Not(:date)])),#xlabel="", #ylabel="", 
            # title="Correlation Heatmap"
        if ann
            for i in 1:size(md, 1)
                for j in 1:size(md, 2)
                    value = round(md[i, j], digits=2)
                    color = ismissing(value) ? :transparent : :black
                    Plots.annotate!(j, i, 
                        Plots.text(string(value), 8, 
                        "Computer Modern",
                        color, 
                        :center, 
                        halign=:center, rotation=-35.0))
                end
            end
        end
        nm = names(df)
        nm = replace.(nm,
            r"_" => " ",r"NaN" => " ",
            r"qoutjl|qout" => "",  r"#" => "")
        Plots.xticks!(1:length(nm), nm) #reverse(nm)
        Plots.yticks!(1:length(nm), nm;rotation=-35)
        return p1
    end

    """
    represents how many standard deviations a given data point is from the mean of its dataset.
    """
    