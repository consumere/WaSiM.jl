#using HTTP
using CSV
using DataFrames
using Dash
#using JSON3

function readdf(x::String)
    # Decode the base64 content
    decoded_content = base64decode(x)

    # Create a buffer from the decoded content
    buffer = IOBuffer(decoded_content)

    # Read the buffer as a CSV file
    #df = CSV.read(buffer, DataFrame)
    ms=["-9999","lin","log"]
    df = CSV.read(buffer, DataFrame,
        missingstring=ms,
        types=Float64,
        delim="\t",
        silencewarnings=false,
        normalizenames=true,
        drop=(i, nm) -> i == 4)
    dropmissing!(df, 1)
    df.YY = map(x -> Int(x), df.YY)
    df.MM = map(x -> Int(x), df.MM)
    df.DD = map(x -> Int(x), df.DD)
    df.date = Date.(string.(df.YY, "-", df.MM, "-", df.DD), "yyyy-mm-dd")
    df = df[:, Not(1:3)]
    DataFrames.metadata!(df, "filename", x, style=:note)
    return df
end

function dfpjs(x)
    df = readdf(x)
    nrows = size(df)[2] - 1
    ti = try
        DataFrames.metadata(df) |> only |> last |> basename
    catch
        @warn "No basename in metadata!"
        ti = raw""
    end
    fig = PlotlyJS.make_subplots(
        shared_xaxes=true,
        shared_yaxes=true
    )
    for i in 1:nrows
        PlotlyJS.add_trace!(fig,
            PlotlyJS.scatter(x=df.date, y=df[:, i],
                name=names(df)[i]))
    end
    fact, logy = 1, 0
    if logy == true
        PlotlyJS.relayout!(fig, yaxis_type="log",
            height=600*fact, width=900*fact,
            title_text="Series of "*ti)
    else
        PlotlyJS.relayout!(fig,
            height=600*fact, width=900*fact,
            title_text="Series of "*ti)
    end
    return fig
end
        
app = Dash.app("app")

app.layout = html_div() do
    html_h1("Upload and Plot Files"),
    dcc_upload(id = "upload-data"),
    html_div(id = "output-graph")
end

app.callback(
    output("output-graph", "children"),
    [input("upload-data", "contents"), input("upload-data", "filename")]
) do contents, filenames
    if contents !== nothing
        graphs = []
        for (content, filename) in zip(contents, filenames)
            # Process the content and generate the plot
            plot_data = dfpjs(content)

            # Create a graph component with the plot data
            graph = dcc_graph(id = "graph", figure = plot_data)

            # Add the graph to the list of graphs
            push!(graphs, graph)
        end

        return graphs
    end
end

run_server(app, "127.0.0.1", 8050)

function readdf_old(x::String)
    ms=["-9999","lin","log"]
    df = CSV.read(x, DataFrame,
        missingstring=ms,
        types=Float64,
        delim="\t",
        silencewarnings=false,
        normalizenames=true,
        drop=(i, nm) -> i == 4)
    dropmissing!(df, 1)
    df.YY = map(x -> Int(x), df.YY)
    df.MM = map(x -> Int(x), df.MM)
    df.DD = map(x -> Int(x), df.DD)
    df.date = Date.(string.(df.YY, "-", df.MM, "-", df.DD), "yyyy-mm-dd")
    df = df[:, Not(1:3)]
    DataFrames.metadata!(df, "filename", x, style=:note)
    return df
end