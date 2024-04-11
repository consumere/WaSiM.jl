# 
function nconly(Any)
        v = readdir();
        z = v[broadcast(x->endswith(x,"nc"),v)];
        return(z)
    end

    