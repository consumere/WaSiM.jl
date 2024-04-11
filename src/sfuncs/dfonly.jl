# 
function dfonly(x1::Regex)
    z = filter(file -> occursin(x1,file), 
    readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
    return(z)
end

