# 
function nctodf(str::String)
        #x = @gl "tsoil"
        rr = Raster(str)
        # dims(rr)|>collect
        # rr.dims|>collect
        # rr[Dim{:t}(Rasters.Between(2,end))] |>contourf
        
        mt = split(string.(rr[:,:,1].refdims)[1],",")|>first
        #occursin("Dim{:t}",)
    
        #if last(dims(rr))=="Dim{:t}"
        if mt=="Dim{:t"
            ti = Rasters.lookup(rr, Dim{:t})
        else
            ti = try
                Rasters.lookup(rr, Ti)
            catch
                @error "no time dimension found!"
                return
            end 
        end
        
        #map(x->mean(x),rr.dims)
        #dims(rr, (Dim{:t})) isa Vector{Float64}
        
        
        #ag = Rasters.aggregate(Rasters.Center(), rr, (Y(20), X(20));)
        #plot(ag)
        #x,y,z = map(x->round(x ./2;digits=0),size(rr))
        x,y = map(x->round(x ./2;digits=0),size(rr)[1:2])
        #x,y,z = map(x->parse.(Int,x ./2),size(rr))
        #length(rr[1,:,:])
        
        df = DataFrame(rr[X=Int(x),Y=Int(y)]',:auto)|>permutedims
        df = hcat(df, parent(ti),makeunique=true)
        #rename!(df,1=>Rasters._maybename(rr),2=>"date")
        rename!(df,1=>Rasters._maybename(rr),2=>"layer")
    
        DataFrames.metadata!(df, "filename", str, style=:note);        
        return df
    
    end

    """
    simulated, observed=df[:,5],df[:,6]
    """
    