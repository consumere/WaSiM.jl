# 
function addplot(specfile::String)
        fl = CSV.read(specfile,DataFrame;limit=4)
        xc = fl[2,5:end]|>collect
        yc = fl[3,5:end]|>collect
        pts = ArchGDAL.IGeometry[]
        for i in 1:length(xc)
            pt = ArchGDAL.createpoint([xc[i],yc[i]])
            #pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
            push!(pts,pt)
        end
        od = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
        Plots.plot!(od.geometry;c=:black,shape = :x) #:utriangle)
        for (i,pt) in enumerate(od.geometry)
            x = ArchGDAL.getx(od.geometry[i], 0)
            y = ArchGDAL.gety(od.geometry[i], 0)
            name = od.name[i]
            #Plots.text(name,:consolas, 8, :black, :bottom, :left)
            #Plots.annotate!(x, y, Plots.text(name, 8, :black, :bottom, :left))
            Plots.annotate!(x, y, Plots.text(name, 8, :black, :bottom, :left))
        end
        Plots.plot!()
    end

    """
    plots a df and applies a function
        default: monmean
        monsum, yrsum, yrmean
        mode can be :bar, :scatter, :line, :steppre, :steppost,
        :hist, :box

        dfm(s;fun=yrsum,mode=:scatter,leg=false)
        dfm(s;fun=monmean,mode=:box,leg=false)

    """
    