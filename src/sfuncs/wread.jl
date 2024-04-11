# 
function wread(x::Regex;skip=3)
    """
    Read wasim ts with DelimitedFiles.readdlm, skipto line 3 
    no header column
    """
    rgx = glob(x)|>first
    println("loading $rgx ...")
    df = DelimitedFiles.readdlm(rgx, '\t', Float64, '\n';
        header=false,skipstart=skip)
    df = DataFrame(df,:auto)
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
    df=df[:,Not(1:4)]
    metadata!(df, "filename", x, style=:note);
end

