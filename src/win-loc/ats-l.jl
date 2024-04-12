using Dash
using CSV
using DataFrames
using PlotlyJS
using Base64
using Dates
using Statistics


app = Dash.dash(
    external_stylesheets=[
        "https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css",
        "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"
    ],
    update_title="Loading..."
)

app.layout = html_div() do
    [
        html_h4("WaSiM Timeseries Analyzer"),    
        dcc_upload(
            id = "upload-data",
            children = [
                html_div("Drag and Drop"),
                html_a("Select Files")
            ],
            style = Dict(
                "width" => "99%",
                "height" => "60px",
                "lineHeight" => "60px",
                "borderWidth" => "1.5px",
                "borderStyle" => "dashed",
                "borderRadius" => "5px",
                "textAlign" => "center",
                "margin" => "10px"
            ),
            multiple = true
        ),
        dcc_loading(
            id = "loading",
            type = "circle",
            children = [
                html_div(id = "output-graph")
            ]
        )
    ]
end

function parse_contents(contents, filename)
    # Read the contents of the uploaded file
    content_type, content_string = split(contents, ',')

    decoded = base64decode(content_string)
    ms = ["-9999.0", "-9999", "lin", "log", "--"]
    df = CSV.File(IOBuffer(decoded); 
        header=1, 
        normalizenames=true,
        missingstring=ms, types=Float64) |> DataFrame
    dropmissing!(df, 1)
    for i in 1:3
        df[!, i] = map(x -> Int(x), df[!, i])
    end
    df.date = Date.(string.(df[!, 1], "-", df[!, 2], "-", df[!, 3]), "yyyy-mm-dd")
    df = df[:, Not(1:4)]
    
    printstyled("generating graphs...\n", color=:blue)
    
    s = Symbol.(filter(x -> !occursin(r"year|date", x), names(df)))
    
    fig = PlotlyJS.make_subplots(shared_xaxes=true, shared_yaxes=true)

    for i in s
        PlotlyJS.add_trace!(fig, 
            PlotlyJS.scatter(x=df.date, y=df[:, i], name=i)
        )
    end

    ti = filename
    fact = .88
    PlotlyJS.relayout!(fig,
        template="seaborn",
        height=650*fact,
        width=1200*fact,
        title_text="",
        xaxis_rangeslider_visible=false,
        updatemenus=[
            Dict(
                "type" => "buttons",
                "direction" => "left",
                "buttons" => [
                    Dict(
                        "args" => [Dict("yaxis.type" => "linear")],
                        "label" => "Linear Scale",
                        "method" => "relayout"
                    ),
                    Dict(
                        "args" => [Dict("yaxis.type" => "log")],
                        "label" => "Log Scale",
                        "method" => "relayout"
                    )
                ],
                "pad" => Dict("r" => 1, "t" => 10),
                "showactive" => true,
                "x" => 0.11,
                "xanchor" => "left",
                "y" => 1.1,
                "yanchor" => "auto"
            ),
        ]
    )
    
    fig_div = dcc_graph(id = filename, figure = fig)
    return fig_div
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
                dcc_loading(
                    id = "loading",
                    type = "circle",
                    children = [
                        parse_contents(content, filename)
                    ]
                )
            ])
            push!(graphs, graph)
        end
        return graphs
    end
end

run_server(app, "127.0.0.1", 8054, debug = true)
