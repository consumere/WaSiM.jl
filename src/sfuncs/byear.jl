# 
function byear(x::Union{String,Regex,DataFrame};)
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
            #df = copy(x)
            df = reorder_df(x)
            dropmissing!(df)
        end
        
        df[!, :year] = year.(df[!,:date]);
        nm = collect(DataFrames.metadata(df))[1][2]|>basename;
        nm = replace.(nm,
            r"_" => " ",
            r"-qoutjl|qout" => "",
            r"#" => "")
        grouped_df = groupby(df, :year)
        
        DataFrames.combine(grouped_df) do group
            #for i in names(group)
            #getproperty(group,))
            #make shure to select only 2 cols without date.
            df = select(group, Not(:date))
            simulated, observed = vec(Matrix(df[!,Cols(1)])),vec(Matrix(df[!,Cols(2)]))
            #simulated, observed = vec(Matrix(df[!,1])),vec(Matrix(df[!,2]))
            local kge = kge2(simulated, observed)
            local nse = nse2(simulated, observed)
            #local nse = (1 - (sum((simulated .- observed).^2) / sum((observed .- mean(simulated)).^2)))
            local ve = vef(simulated, observed)
            #grouping key :year is returned as first column
            dout = DataFrame(kge=kge,nse=nse,ve=ve,nm=nm)
            return dout
        end
    end

    """
    lat lon reverse
    """
    