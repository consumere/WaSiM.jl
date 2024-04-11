# 
function stp(fn::String;repro=false)
        fl = CSV.read(fn,DataFrame;limit=4)
        xc = fl[2,5:end]|>collect
        yc = fl[3,5:end]|>collect
        pts = ArchGDAL.IGeometry[]
        if repro #projects to 4326
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
        else
            for i in 1:length(xc)
                pt = ArchGDAL.createpoint([xc[i],yc[i]])
                #pt = AG.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
                push!(pts,pt)
            end
        end
        nd = DataFrame(geometry=pts, name=propertynames(fl)[5:end], xc=xc, yc=yc)
        return nd
    end

    stread=stp

    """
    reads from stationdata, reprojects and plots
    """
    