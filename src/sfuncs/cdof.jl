# 
function cdof(x::Union{String,DataFrame})
    if x isa DataFrame
        try
            d = collect(DataFrames.metadata(x))[1][2]
            cd(dirname(d))
        catch
            @error "no basename in present!"
            return nothing
        end
    else
        try
            cd(dirname(x))
        catch
            @error "cd errored!"
            return nothing
        end
    end
    d=pwd()
    println("current dir: $d")
end

#k=raw"C:/Users/Public/Documents/Python_Scripts/julia/smallfuncs.jl"
k=src_path*"/smallfuncs.jl"
println("script loc is $k")
#homedir()|>cd
#home() ##necessary for ssup to work

# @vv "surface"
# x=raw"C:/Users/chs72fw/Documents/Promotionsstudium/Dropbox/brendfab/m3/Soil_Temperature_Stack.nc"
# r = Raster(x)
# r = readras(x)
# plot(r[t=4:5])
# surface(r[t=5],camera = (0,-75),legend=false)

