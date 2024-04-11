# 
function fparse(nm)
        [parse.(Float64, split(string(x))) for x in nm]
    end

    """
    Convert DataFrame Column to a Vector
    returns only first match, see tovec for multiple matches
    or:
    m=Symbol.(filter(x->occursin(r"k",x),names(df)))
    map(x->getproperty(df, x),m)
    DataFrame(map(x->getproperty(df, x),m),:auto)
    """
    