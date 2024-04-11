# 
function pyread_meteo(s::AbstractString;hdr=0)
    #     pd = pyimport("pandas")
    #     ddf = pd.read_csv(s, delim_whitespace=true, 
    #         header=hdr,
    #         na_values=-9999,
    #         low_memory=false,
    #         verbose=true)
    #     #ddf.filename=basename(s)
    #     ddf = wa.pydf(ddf)
        
    #     #dt = Date.(map(x -> join(x, "-"), eachrow(ddf[:, 1:3])))
    #     dt = select(ddf,1:3)
    #     Z = []
    #     for i in eachcol(dt)
    #         col = tryparse.(Int, i)
    #         push!(Z,col)
    #     end
    #     dt = DataFrame(Z,:auto)
        
    #     msk = tryparse.(Int, ddf.YY)
    #     ddf = ddf[findall(!isnothing,msk),:]
    #     dt = dt[findall(!isnothing, dt.x1),:]
    #     dt = Date.(map(k -> join(k, "-"), eachrow(dt[:, 1:3])))
    #     nd = hcat(ddf[!, 5:end], dt)
       
    #     for col in names(nd)[(eltype.(eachcol(nd)) .<: String)]
    #         nd[!, col] .= tryparse.(Float64, col)
    #     end
    #     rename!(nd, ncol(nd) =>"date")
    #     metadata!(nd, "filename", s, style=:note);
    #     return nd
    # end

    # 