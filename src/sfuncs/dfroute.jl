# 
function dfroute(;ofl="route.txt")
    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"_>.*" => "",r"-" => "_"),df[:,3])
    sort!(df, :sim)
    return df    
end

"""
non-recursively search for control file in current directory
"""
