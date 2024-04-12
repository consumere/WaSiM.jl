using Dash
using CSV
using DataFrames
using PlotlyJS
using Base64
using Dates
using Statistics

###doesnt work.

app = Dash.dash(;
external_scripts=["https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML"],
update_title="Loading...")

app.layout = html_div() do
    [
        html_h1("WaSiM Timeseries Data"),
        dcc_upload(
            id = "upload-data",
            children = [
                html_div("Drag and Drop or "),
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
        html_div(id = "output-graph")
        #,
        #html_button("Aggregate to Year Sums", id = "aggregate-button", n_clicks = 0)
    ]
end

function parse_contents(contents, filename)
    # Existing code for parsing contents
    content_type, content_string = split(contents, ',')

    # Decode the file content
    decoded = base64decode(content_string)
    ms = ["-9999.0", "-9999", "lin", "log", "--"]
    df = CSV.File(IOBuffer(decoded); delim="\t", header=1, normalizenames=true,
        missingstring=ms, types=Float64) |> DataFrame
    dropmissing!(df, 1)
    for i in 1:3
        df[!, i] = map(x -> Int(x), df[!, i])
    end
    df.date = Date.(string.(df[!, 1], "-", df[!, 2], "-", df[!, 3]), "yyyy-mm-dd")
    df = df[:, Not(1:4)]

    function yrsum(x::DataFrame)
        df = copy(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        return(df_yearsum)
    end

    df_yrsum = yrsum(df)





    fig = PlotlyJS.make_subplots(shared_xaxes=true, shared_yaxes=true)
    # Existing code for creating the figure
    tcols = size(df)[2] - 1

    for i in 1:tcols
        PlotlyJS.add_trace!(fig, 
        PlotlyJS.scatter(x=df.date, y=df[:, i], name=names(df)[i])
        )
    end



    ti = filename
    fact = 1.1
    #fact = .8

#    Add dropdown
    PlotlyJS.relayout!(fig,
    template="seaborn",
    #template="simple_white",
    #template="plotly_dark",
    height=650*fact,
    width=1200*fact,
    title_text=filename,
    #xperiod = first(df.date),
    #xperiodalignment = "start",
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
                ),
                Dict(
                    "args" => [Dict("yaxis.type" => "linear"), 
                    Dict("y" => df_yrsum[:, 2]),
                    Dict("x" => df_yrsum[:,:year])],
                    "label" => "Year Sums",
                    "method" => "update" ,
                    # "args2" => [Dict("visible" => [true, true] + [false] * (tcols - 1))],
                    "args3" => [Dict("title" => "Year Sums")]
                )
            ],
            "pad" => Dict("r" => 1, "t" => 10),
            "showactive" => true,
            "x" => 0.11,
            #"x" => 5.11,
            "xanchor" => "left",
            #"xanchor" => "auto",
            "y" => 1.1,
            #"yanchor" => "top"
            "yanchor" => "auto"
        ),
    ]
    )
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

run_server(app, debug=true)

