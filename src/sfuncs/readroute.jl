# 
function readroute(x::Union{String,AbstractString})
        df = CSV.read(x,DataFrame,header=false,
            skipto=8,delim="\t",footerskip=1,lazystrings=false)
        rename!(df,1=>"sim",2=>"obs",3=>"name")
        df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
        df.name=map(x->replace(x,r"_>.*" => ""),df.name)
        sort!(df, :sim)
        return df
    end

    """
    rename_columns!(df::DataFrame, name_mapping::DataFrame)
    inplace
    """
    