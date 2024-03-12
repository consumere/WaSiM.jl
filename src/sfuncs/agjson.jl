# 
function agjson(jsonfile::AbstractString)
        """
        reads json and transforms points...
        but AG transformation is wrong :(


        ERROR: type UnionAll has no field read
        works only in WSL
        """
        #const AG = ArchGDAL

        geojson_file=jsonfile
        jsonbytes = read(geojson_file) # read the geojson file as bytes
        fc = GeoJSON.read(jsonbytes)
        pts=[]
        for geom in fc.geometry
            xc = [(x) for x in geom]|>first|>first
            yc = [(x) for x in geom]|>first|>last
            pt = ArchGDAL.createpoint(xc,yc)
            pt = ArchGDAL.reproject(pt,EPSG(25832),ProjString("+proj=longlat +datum=WGS84 +no_defs"))
            push!(pts,pt)
        end
        # df = [] ##geht auch
        # for x in pts
        #     x1 = AG.getpoint(x,0)
        #     tmp = DataFrame(x=[x1[1]], y=[x1[2]])
        #     push!(df,tmp)
        # end
        # df = reduce(vcat, df) 
        #Plots.plot(df.x, df.y, seriestype=:scatter)
        df = DataFrame( 
            "x" => [AG.getpoint(x,0)[1] for x in pts],
            "y" => [AG.getpoint(x,0)[2] for x in pts] )
        return(df)
    end

    # works only in current REPL -> see pwc()
    # import InteractiveUtils.clipboard as cb
    # wslpath()|>cb

    