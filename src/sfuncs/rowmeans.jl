# 
function rowmeans(df::DataFrame)
        return DataFrames.combine(df, names(df) .=> mean)
    end
    
    #all to one mean: df[!,:mean] = mean.(eachrow(df[!,Not(:year)]))

    """
    addname(indf::DataFrame,nmdf::DataFrame)
    indf: input dataframe
    nmdf: name dataframe

    indf = dfonly(r"gwst")|>first|>waread

    ofl="route.txt"
    begin
        df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
        rename!(df,1=>"sim",2=>"obs",3=>"name")
        df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
        df.name=map(x->replace(x,r"_>.*" => ""),df.name)
        sort!(df, :sim)
    end
    
    addname(indf,df)
    select(indf,Not(Cols(r"^[0-9]")))|>dfp

    """
    