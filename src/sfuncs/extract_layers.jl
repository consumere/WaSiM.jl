# 
function extract_layers(df::DataFrame, colkey::String="ksat")
        num_layers = size(select(df,Cols(colkey))[1,1],1)
        layer_columns = [Symbol("layer$i") for i in 1:num_layers]
        layer_values = Vector{Vector{Float64}}(undef, num_layers)   
        for i in 1:num_layers
            layer_values[i] = [ r[1][i] for r in eachrow(select(df, colkey)) ]
        end    
        df1 = DataFrame(layer_columns .=> layer_values)
        dout = rename(hcat(df.key,df1), 1 => :key)
        return dout    
    end

    """
    a reader using filtermask and looks_like_number
    kw passed to CSV.read
    """
    