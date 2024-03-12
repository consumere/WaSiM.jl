# 
function reverse_coords(polygon)
        # Extract the coordinates from the polygon
        # crds = GeoInterface.coordinates.(polygon)
        # # Reverse each pair of coordinates
        # reversed_crds = [[[(lon, lat) for (lat, lon) in ring] for ring in polygon] for polygon in crds]
        
        # coords = GeoInterface.coordinates.(gd.geometry)
        coords = GeoInterface.coordinates(polygon)
        # Flatten the list of points
        #points = reduce(vcat, coords[1])
        ptc = coords[1]
        #ptc = GeoInterface.coordinates.(points)
        ptc_tuples = [tuple(pt...) for pt in ptc]
        prev = [(Y,X) for (X,Y) in ptc_tuples]
        #
        #rpt = ArchGDAL.createpoint.(prev)
        #reversed_polygon = ArchGDAL.createpolygon(coords)
        reversed_polygon = ArchGDAL.createpolygon(Vector(prev))
        
        # Create a new polygon with the reversed coordinates
        #reversed_polygon = ArchGDAL.createpolygon.(reversed_crds)
        #reversed_polygon = ArchGDAL.createpolygon.(prev)

        return reversed_polygon
    end

    """
    plots timeseries like cmk, but with a legend
    https://docs.makie.org/stable/reference/blocks/legend/
    """
    