# 
function vecdf(x::DataFrame, col::Any)
        m = try
                propertynames(df)[findfirst(x->occursin(col,x),names(df))]
            catch
                @error "no match found"
                return
            end
            @info "$m found!"
        return getproperty(x, m)
    end
    """
    extract_layers(df::DataFrame, colkey::String="ksat")
    form fsoil(infile)
    or build_soil_dictionary(infile) -> df
    """
    