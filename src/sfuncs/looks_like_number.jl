# 
function looks_like_number(str::AbstractString)
    try
        parse(Float64, str)
        return true
    catch
        return false
    end
end

