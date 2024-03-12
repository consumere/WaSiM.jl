# 
function stplot!(fn::Union{String,DataFrame})
        if isa(fn,DataFrame)
            fl = fn
            xc = fl.xc
            yc = fl.yc
            pts = ArchGDAL.IGeometry[]
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
            #trim chars to length 10
            len = 10
            mns = map(y->((length(y)> len)  ? y[1:len] : y),string.(fl.name))
            od = DataFrame(geometry=pts, name=mns, xc=xc, yc=yc)
        else
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
            
        end
        
        Plots.plot!(od.geometry;c=:black,shape = :x)
        for (i, pt) in enumerate(od.geometry)
            x = ArchGDAL.getx(od.geometry[i], 0)
            y = ArchGDAL.gety(od.geometry[i], 0)
            name = od.name[i]
            annotate!(x, y, text(name, 8, :black, :bottom, :left))
        end
        Plots.plot!()
    end

    