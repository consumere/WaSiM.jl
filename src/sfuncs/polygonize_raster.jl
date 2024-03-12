# 
function polygonize_raster(input_raster_path::String, output_shapefile_path::String;epsg=25832)
        gdal = pyimport("osgeo.gdal")
        ogr = pyimport("osgeo.ogr")
        osr = pyimport("osgeo.osr")

        # Open the raster dataset
        dataset = gdal.Open(input_raster_path)

        # Get the first band
        band = dataset.GetRasterBand(1)

        # Get the "ESRI Shapefile" driver
        driver = gdal.GetDriverByName("ESRI Shapefile")

        # Create a new shapefile dataset
        out_ds = driver.Create(output_shapefile_path, 0, 0, 0, gdal.GDT_Unknown)

        # Create a spatial reference object
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)

        # Create a new layer
        layer = out_ds.CreateLayer("polygonized", srs, ogr.wkbPolygon)

        # Polygonize the raster
        try
            gdal.Polygonize(band, py"None", layer, -1, [], callback=py"None")
            @info "new shapefile created at: $output_shapefile_path ..."
        catch
            @error "gdal.Polygonize failed ..."
        end
        
        # Close the dataset to write it to the disk
        out_ds = py"None"
    end



    """
    reads controlfile
    uses Grep.grep to select lines
    returns a DataFrame
    """
    