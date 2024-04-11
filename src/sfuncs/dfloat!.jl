# 
function dfloat!(df::DataFrame)
    for colname in names(df)
        if eltype(df[!, colname]) != Date && eltype(df[!, colname]) == String
            df[!, colname] .= tryparse.(Float64, df[!, colname])
        end
    end
end

#endof func declaration

println("you are here: ")
printstyled(pwd()*"\n",color=:green)

## if gr fails, try to plot GR.histogram(randn(1000))
# import GR
# 