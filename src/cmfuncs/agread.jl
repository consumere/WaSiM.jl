# 
function agread(x::String)
        dataset = ArchGDAL.read(x)
        # Get the first band
        band = ArchGDAL.getband(dataset, 1)
        # A = permutedims(band, (2, 1))
        # A = A.a
        # A = reverse(A, dims=1)
        A = reverse(band, dims=2)
        B = replace(A, -9999.0 => missing)
        xm = .!all(ismissing, B, dims=2)
        ym = .!all(ismissing, B, dims=1)
        A = B[xm[:], ym[:]]
        k = Raster(A,(X,Y))
        println(extrema(k))
        return k
    end  

    """
    streamplot
    streamplot(f::function, xinterval, yinterval; color = norm, kwargs...)
    """                     
    