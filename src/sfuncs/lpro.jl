# 
function lpro(x::Union{String,DataFrame})
        if x isa String

            if endswith(x,"shp")
                fl = GeoDataFrames.read(x)
                @info "reprojecting polygon..."
                geom=fl.geometry|>first
                out = ArchGDAL.reproject(geom,
                    GeoFormatTypes.EPSG(25832),
                    GeoFormatTypes.EPSG(4326))
                od = DataFrame(geometry=out, 
                    name=propertynames(fl)[end])
               return od 
            else
                fl = CSV.read(x,DataFrame;limit=4)
            end

        else
            fl = x
        end
        xc = fl[2,5:end]|>collect
        yc = fl[3,5:end]|>collect
        pts = ArchGDAL.IGeometry[]
        for i in 1:length(xc)
            pt = ArchGDAL.createpoint([xc[i],yc[i]])
            pt = ArchGDAL.reproject(pt,
            GeoFormatTypes.EPSG(25832),
            GeoFormatTypes.EPSG(4326))
            ##reverse cords.
            ptc = tuple(GeoInterface.coordinates(pt)...)
            prev = ptc[2],ptc[1]
            pt = ArchGDAL.createpoint(prev)
            push!(pts,pt)
        end
        od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
    
        return od    
    end

    """
    only for wasim output files
    """
    