# 
function addname(indf::DataFrame,nmdf::DataFrame)
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
        
        addname(indf,nmdf::df)
        select(indf,Not(Cols(r"^[0-9]")))|>dfp
    
        """
        rename!(indf,map(x->replace(x,r"^C" => ""),names(indf)))
        for row in eachrow(nmdf)
            old_name = string(row.sim)
            new_name = string(row.name)
            if old_name in names(indf)
                println("renaming $old_name")
                rename!(indf, (old_name) => (new_name*"_"*old_name))
            end
        end
        return indf
    end

    