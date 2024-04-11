# 
function dfon(x1::AbstractString)
    z = filter(file -> occursin(Regex(x1,"i"),file), 
    readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
    return(z)
end

# x=nothing
# typeof(x)
# isa(x,Any) #true
#Nothing,

