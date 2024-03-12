# 
function findalls(x::Vector)
        map(e -> e => findall(==(e), x), unique(x))
    end

    