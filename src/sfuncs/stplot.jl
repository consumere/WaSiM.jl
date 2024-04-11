# 
function stplot(fn::String)
        fl = CSV.read(fn,DataFrame;limit=4)
        xc = fl[2,5:end]|>collect
        yc = fl[3,5:end]|>collect
        pts = ArchGDAL.IGeometry[]
        for i in 1:length(xc)
            pt = ArchGDAL.createpoint([xc[i],yc[i]])
            pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
            push!(pts,pt)
        end
        od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
        p = plot(od.geometry);
        for (i, pt) in enumerate(od.geometry)
            #x = od.xc[i]
            x = ArchGDAL.getx(od.geometry[i], 0)
            #y = od.yc[i]
            y = ArchGDAL.gety(od.geometry[i], 0)
            name = od.name[i]
            annotate!(x, y, text(name, 8, :black, :bottom, :left))
        end
        return plot!(p)
    end

    """
    runs inside REPL
    """
    