# 
function ctov(df::DataFrame)
        df = df[!,Not(Cols(r"date|year|time"))]|>dropmissing   
        nr = DataFrames.nrow(df)
        
        col_names = String[]
        # for str in names(df)
        #     a = [str for i in 1:nr]
        #     push!(col_names, a)
        # end
    
        for str in names(df)
            a = [str for i in 1:nr]
            append!(col_names, a)
        end
    
        col_vectors = Float64[]
    
        for col in eachcol(df)
            append!(col_vectors, col)
        end
        return col_names, col_vectors
    end

    #https://docs.makie.org/stable/examples/plotting_functions/rainclouds/


    """
    cloudplot(df::DataFrame)
    #https://docs.makie.org/stable/examples/plotting_functions/rainclouds/
    """
    