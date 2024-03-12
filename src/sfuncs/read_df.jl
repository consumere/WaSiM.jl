# 
function read_df(s::Union{String,Regex})
    """
    reads a DataFrame from a file w dlm and tryparse subsetting
    """
    if s isa Regex
        s = dfonly(s)|>first
    end
    data, colnames = DelimitedFiles.readdlm(s, '\t', String, '\n', header=true)
    df = DataFrame(Matrix{Any}(data), :auto)
    rename!(df,Symbol.(colnames)|>vec)
    df = df[map(x->!isnothing(x),tryparse.(Int,df[!,1])),:]
    for i in 1:size(df, 2)
        df[!, i] = map(x -> parse(Float64, x), df[!, i])
    end
    for i in 5:size(df,2)
        df[!,i]=replace(df[!,i],-9999.0 => missing)
    end 
    for i in 5:size(df,2)
        replace!(df[!,i],-9999.0 => missing)
    end
    
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df = df[:,Not(1:4)]
    metadata!(df, "filename", s, style=:note);
    return df
end

