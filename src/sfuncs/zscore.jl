# 
function zscore(x::Vector{Float64})
        μ = mean(x)
        σ = std(x)
        z = (x .- μ) ./ σ
        return z
    end