using Dash
using PlotlyJS
#using Plots
# using Rasters
# import NCDatasets
using NCDatasets
using Base64

app = Dash.dash(;update_title="Loading...")

app.layout = html_div([
    dcc_upload(
        id = "upload-data",
        children = html_div([
            "Drag and Drop or ",
            html_a("Select Files")
        ]),
        style = Dict(
            "width" => "80%",
            "height" => "60px",
            "lineHeight" => "60px",
            "borderWidth" => "2px",
            "borderStyle" => "dashed",
            "borderRadius" => "5px",
            "textAlign" => "center",
            "margin" => "10px"
        ),
        multiple = true
    ),
    html_div(id = "output-graph")
])

function parse_contents(contents, filename)
    # Read the contents of the uploaded file
    content_type, content_string = split(contents, ',')
    # Decode the file content
    decoded = base64decode(content_string)
    # Write the decoded content to a temporary file
    #temp_file = tempname(tempdir = tempdir(), prefix = "", suffix = ".nc")
    temp_file = tempname()*".nc"
    write(temp_file, decoded)
    # #close(temp_file)
    
    # #temp_file = raw"D:\Wasim\Tanalys\DEM\brendpest\out_v4\AnnualTemperature.nc"
    # # Read the NetCDF file using NCDatasets
    # dx = Dataset(temp_file)
    # #dx = Dataset(x)
    # dict = dx|>Dict
    # mykeys = keys(dict)
    # varname = filter(x->!occursin(r"^t$|^y|^x",x)    
    #             ,string.(mykeys))|>first
    # # Get X and Y coordinates
    # x_coords = dx["x"][:]
    # y_coords = dx["y"][:]
    # # Transpose the raster data
    # Z = permutedims(dx[varname][:,:,end], [2, 1])
    # #Z = coalesce.(Z, missing)
    # #M = Matrix(x_coords,y_coords,Z)
    
    # fig = PlotlyJS.plot(
    #     [
    #         PlotlyJS.heatmap(
    #             x = x_coords,
    #             y = y_coords,
    #             z = Z,
    #             colorscale = "YlOrRd"
    #         )])
    # if isfile(temp_file)
    #     r=read(Raster(temp_file,missingval=0,mappedcrs=EPSG(25832)));
    #     if (r.dims[end]|>length == 1)
    #         @warn("only one layer available...")
    #         #fig=plot(r;)
    #         # xlabel="",
    #         # ylabel="",
    #         #c=cgrad(:thermal),
    #         # c=cgrad(:matter),
    #         # size=(1200, 800));
    #     else
    #         #@warn("no nc file parsed...")
    #         @warn("subsetting first layer...")
    #         ee = Int(r.dims[3][end])
    #         r = r[t=2:ee];    #subset till end
    #         # p=Plots.plot(rn;
    #         # xlabel="",
    #         # ylabel="",
    #     end
    # end
    # x_coords = Rasters.dims(r)[1]|>collect
    # y_coords = Rasters.dims(r)[2]|>collect
    # # # Transpose the raster data
    # #Z = permutedims(r[:,:,end], [2, 1])|>collect
    # Z = r[:,:,end]|>collect
    
    # fig = PlotlyJS.make_subplots(
    #     shared_xaxes=true, 
    #     shared_yaxes=true    
    #     );
# #    M = r|>collect

#     #zv = reshape(Z, size(Z)[1]*size(Z)[2])


#     # PlotlyJS.add_trace!(fig, 
#     #     PlotlyJS.heatmap(
#     #                     x = 1:length(x_coords),
#     #                     y = 1:length(y_coords),
#     #                     z = zv,
#     #                     colorscale = "YlOrRd")
#     #                     )

#     PlotlyJS.add_trace!(fig, 
#                         PlotlyJS.heatmap(
#                                         x = x_coords,
#                                         y = y_coords,
#                                         z = Z,
#                                         colorscale = "YlOrRd"))

#     fact=0.7
#     PlotlyJS.relayout!(fig,
#                         template="seaborn",
#                         height=600*fact,width=900*fact,
#                         title_text = "",
#     xaxis_title = "X",
#     yaxis_title = "Y"     )
    
#     # Stacktrace:
#     # [1] unsafe_string
#     #   @ .\strings\string.jl:85 [inlined]
#     # fig = Plots.plot(r;
#     #             # xlabel="",
#     #             # ylabel="",
#     #             #c=cgrad(:thermal),
#     #             c=cgrad(:matter),
#     #             size=(1200, 800));
    
#     return fig
# #    rm(temp_file;force=true)

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



    fig = ncplotjl(temp_file)
    
    return fig
end

callback!(
    app,
    Output("output-graph", "children"),
    [Input("upload-data", "contents")]
) do contents, filenames, value
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

run_server(app, "0.0.0.0", 8052, debug = true)
