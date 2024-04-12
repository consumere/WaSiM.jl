using Dash
using CSV
using DataFrames
using PlotlyJS
using Base64
using Dates

app = Dash.dash()

app.layout = html_div() do
    [
    html_h1("WaSiM timeseries data"),    
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
    ms = ["-9999.0", "-9999", "lin", "log", "--"]
    # df = CSV.read(IOBuffer(decoded), DataFrame;
    #     delim="\t",header=1,
    #     silencewarnings=true,
    #     missingstring=ms,
    #     types=Float64)
    df = CSV.File(IOBuffer(decoded); delim="\t", header=1, normalizenames=true, 
        missingstring=ms, types=Float64) |> DataFrame
    
    dropmissing!(df,1) #better than above

    # for i in 5:size(df,2)
    #     df[!,i]=replace(df[!,i],-9999.0 => missing)
    # end 
    # for i in 5:size(df,2)
    #     replace!(df[!,i],-9999.0 => missing)
    # end
    # map to int for dates
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df=df[:,Not(1:4)]
    
    fig = PlotlyJS.make_subplots(
        shared_xaxes=true, 
        shared_yaxes=true    
        );
    #ti = " " #placeholder
    ti = filename
    
    nrows=size(df)[2]-1
    
    for i in 1:nrows;
        PlotlyJS.add_trace!(fig, 
        #PlotlyJS.bar(
        PlotlyJS.scatter(   
            x=df.date, y=df[:,i],
            name=names(df)[i]));
    end
    
    fact, logy = 1.05, false
    if logy == true
        PlotlyJS.relayout!(fig,
        template="seaborn",
        yaxis_type="log",
        height=600*fact,width=900*fact,
        #title_text="Series of "*ti)
        title_text=ti)
    else
        PlotlyJS.relayout!(fig,
        template="seaborn",
        height=600*fact,width=900*fact,
        #title_text="Series of "*ti)
        title_text=ti)
    end
   
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
                #html_h4(filename),
                dcc_graph(
                    id = filename,
                    figure = parse_contents(content, filename) #, filename
                )
            ])
            push!(graphs, graph)
        end
        return graphs
    end
end

#run_server(app, "127.0.0.1", 8050)
run_server(app, debug=true)