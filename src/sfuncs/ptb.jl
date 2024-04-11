# 
function ptb(x::DataFrame)
    strs = sprint(io -> pretty_table(io, x, header=uppercasefirst.(names(x)), backend = Val(:text)))
    clipboard(strs)
end

