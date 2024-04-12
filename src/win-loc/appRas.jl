using Dash
#using Rasters
using Plots
#using PlotlyJS
using Base64
#using Dates
#import NCDatasets
using NCDatasets
#plotlyjs()

app = Dash.dash(;update_title="Loading...")

app.layout = html_div() do
    [
        dcc_upload(
            id = "upload-data",
            children = [
                html_div("Drag and Drop or "),
                html_a("Select Files")
            ],
            style = Dict(
                "width" => "100%",
                "height" => "60px",
                "lineHeight" => "60px",
                "borderWidth" => "1px",
                "borderStyle" => "dashed",
                "borderRadius" => "5px",
                "textAlign" => "center",
                "margin" => "10px"
            ),
            multiple = true
        ),
        html_div(id = "output-graph")
    ]
end

function parse_contents(contents, filename)
    # Read the contents of the uploaded file
    content_type, content_string = split(contents, ',')
    # Decode the file content
    decoded = base64decode(content_string)
    #xr = read(Raster(IOBuffer(decoded);crs=EPSG(25832),missingval=0))
    #M = read(IOBuffer(decoded)) #::Vector{UInt8}
    #M = decoded
    #println(first(M))
    
    # # Save the decoded content to a temporary file
    tmpfile = tempname()
    tmpfile = tmpfile*".nc"
    buf = IOBuffer(decoded)
    vec = take!(buf) # convert IOBuffer to Vector{UInt8}
    write(tmpfile, vec)
    nc = NCDataset(tmpfile)

    #M = raw"D:\Wasim\regio\out\c8\ei__rcm_8100.sum.nc"
    #xr = read(Raster(tmpfile;crs=EPSG(25832),missingval=0))
    #Plots.plot(Raster(M))

    # buf = IOBuffer(decoded)
    
    # vec = take!(buf) # convert IOBuffer to Vector{UInt8}
    # nc = NCDataset(vec)
    #str = String(vec) # convert Vector{UInt8} to String
    #str=raw"D:\Wasim\regio\out\v9\ei__rcm.2017.nc"
    #nc = NCDataset(str)


    dict = nc|>Dict
    mykeys = keys(dict)
    #println(string.(mykeys))
    #v = filter(x->!occursin(r"time|lon|lat|x|y|spatial_ref",x),string.(mykeys))|>first
    #time = nc["time"][:]
    # datetime_vector = coalesce.(time, missing)
    # df = DataFrame(
    #         v=>nc[v][50,50,:],      #x, y, indices
    #         "date"=>datetime_vector)
    # DataFrames.metadata!(df, "filename", x, style=:note);    

    # lon = nc["x"][:]
    # latl = nc["y"][:]
    # #time = nc["time"][:]   
    # time = nc["t"][:]
    #
    varname = filter(x->!occursin(r"^t$|^y|^x",x)    
                ,string.(mykeys))|>first
    
    #rvar = nc[varname][:]

#    rvar[:,:,1]'|>Plots.plot


# fig = Plots.heatmap(lon, latl, rvar[:,:,1]', 
#             xlabel = "Longitude", ylabel = "Latitude", 
#             plot_title = string(varname));

    fig = Plots.heatmap(nc[varname][:,:,end]',
            xlabel = "Longitude", 
            ylabel = "Latitude", 
            #xticks = lon, yticks = latl,
            legend_title = string(varname),
            xflip=true)

# M = nc[varname][:,:,end]
# perm = (2, 1);
# mx = permutedims(M,perm)
# Plots.heatmap(mx,
#             xlabel = "Longitude", ylabel = "Latitude", 
#             legend_title = string(varname))


    #xr = read(Raster(String(decoded);crs=EPSG(25832),missingval=0))
    #xr = Raster(IOBuffer(decoded))

    #xr = xr[Dim{Rasters.name(xr.dims)[end]}(Rasters.Where(x -> x <= 6))] 
        #|>Plots.plot
    #fig = Plots.plot(xr;c=cgrad(:thermal),size=(1200*.8, 800*.8))
    #fig = plot(xr)
    #fig = Plots.plot(xr[t=1];c=cgrad(:thermal),size=(1200*.8, 800*.8))
    return fig
end


callback!(
    app,
    Output("output-graph", "children"),
    [Input("upload-data", "contents")],
    [State("upload-data", "filename")]
) do contents, filenames
    if contents !== nothing
        graphs = []
        for (content, filename) in zip(contents, filenames)
            graph = html_div([
                html_h4(filename),
                dcc_graph(
                    id = filename,
                    figure = parse_contents(content, filename)
                )
            ])
            push!(graphs, graph)
        end
        return graphs
    end
end

run_server(app, "0.0.0.0", debug=true)


#run_server(app, "127.0.0.1", 8050)

# function rplot(reg::Regex)
#     """
#     rplot(x::Regex)
#     reads first match of regex wasim ncs
#     """
#     file = filter(x -> occursin(reg,x), readdir(pwd()))
#     println("subsetting first nc of $file...")
#     file = filter(x -> endswith(x,".nc"), file)|>first
#     xr = read(Raster(file;crs=EPSG(25832),missingval=0))
# 	Plots.plot(xr[t=1];c=cgrad(:thermal),size=(1200*.8, 800*.8))   
# end





function ncplotjl(temp_file)
    dx = Dataset(temp_file)
    #dx = Dataset(x)
    dict = dx|>Dict
    mykeys = keys(dict)
    varname = filter(x->!occursin(r"^t$|^y|^x",x)    
                ,string.(mykeys))|>first
    # Get X and Y coordinates
    x_coords = dx["x"][:]
    y_coords = dx["y"][:]
    # Transpose the raster data
    Z = permutedims(dx[varname][:,:,end], [2, 1])
 
    fig = PlotlyJS.make_subplots(
        shared_xaxes=true, 
        shared_yaxes=true    
        );


    PlotlyJS.add_trace!(fig, 
                        PlotlyJS.heatmap(
                                        x = x_coords,
                                        y = y_coords,
                                        z = Z,
                                        colorscale = "YlOrRd"))

    fact=0.7
    PlotlyJS.relayout!(fig,
                        template="seaborn",
                        height=600*fact,width=900*fact,
                        title_text = "",
    xaxis_title = "X",
    yaxis_title = "Y"     )

    return fig
end

str=raw"D:\Wasim\streu\out\c1\albestr.2013.nc"
ncplotjl(str)