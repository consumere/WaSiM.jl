using Dash
using CSV
using DataFrames
using PlotlyJS
using Base64
using Dates


app = Dash.dash(
    external_stylesheets=[
        "https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css",
        "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"
    ],
    update_title="Loading..."
)
#methods(dash)
#methods(html_div)

app.layout = html_div() do
    [
    html_h3("WaSiM Timeseries Data"),    
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
    ]
end

function parse_contents(contents, filename)
    # Read the contents of the uploaded file
    content_type, content_string = split(contents, ',')

    # Decode the file content
    decoded = base64decode(content_string)

    ms = ["-9999.0", "-9999", "lin", "log", "--"]

    df = CSV.File(IOBuffer(decoded); delim="\t", header=1, normalizenames=true,
        missingstring=ms, types=Float64) |> DataFrame

    # df = CSV.read(IOBuffer(decoded), DataFrame;
    #     delim="\t",header=1,
    #     silencewarnings=true,
    #     missingstring=ms,
    #     types=Float64)
    
    
    dropmissing!(df, 1)

    for i in 1:3
        df[!, i] = map(x -> Int(x), df[!, i])
    end

    df.date = Date.(string.(df[!, 1], "-", df[!, 2], "-", df[!, 3]), "yyyy-mm-dd")
    #df.date = Date.(string.(df[!, 1], "\t", df[!, 2], "\t", df[!, 3]), "yyyy\tMM\tDD")

    df = df[:, Not(1:4)]
    fig = PlotlyJS.make_subplots(shared_xaxes=true, shared_yaxes=true)

    # tcols = size(df)[2] - 1
    
    # for i in 1:tcols
    #     PlotlyJS.add_trace!(fig, 
    #     PlotlyJS.scatter(x=df.date, y=df[:, i], name=names(df)[i])
    #     )
    # end
    s = (filter(x->!occursin(r"year|date",x),names(df)))
    #renamer - remove char _   
    for x in s
        newname=replace(x,"_"=>" ")
        rename!(df,Dict(x=>newname))
    end
    s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
    for i in s
        PlotlyJS.add_trace!(fig, 
        PlotlyJS.scatter(x=df.date, y=df[:, i], name=i)
        )
    end
    
    ti = filename
    fact = 1.1
    #fact = .8

    PlotlyJS.relayout!(fig,
    template="seaborn",
    #template="simple_white",
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

#parse_y
function parse_y(contents, filename)
    # Read the contents of the uploaded file
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

    df = yrsum(df)
    
    fig = PlotlyJS.plot(df, kind = "bar");
    #rename(df,replace(names(df),"_1"=>""))
    s = (filter(x->!occursin(r"year|date",x),names(df)))
    #renamer - remove chars   
    for x in s
        newname=replace(x,"_"=>" ")
        #println(newname)
        rename!(df,Dict(x=>newname))
    end
    s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
    
    for i in s;
        PlotlyJS.add_trace!(fig, 
        PlotlyJS.bar(x=df.year, y=df[:,i],
        name=i)       );
    end
    #filename="lala.dka"
    ti = split(filename,".")|>first
    #typeof(ti)
    fact = 1.1
    PlotlyJS.relayout!(fig,
    template="seaborn",
    #template="simple_white",
    #template="plotly_dark",
    height=650*fact,
    width=1200*fact,
    title_text="yearly sums of "*ti,
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
                dcc_graph(
                    id = filename,
                    figure = parse_contents(content, filename) #, filename
                )
            ])
            push!(graphs, graph)
            yrs = html_div([
            dcc_graph(
                id = filename,
                figure = parse_y(content, filename)
                )
            ])
            push!(graphs, yrs)
        end
        return graphs
    end
end


run_server(app, "0.0.0.0", 8052, debug = true)
#run_server(app, debug=true)
#run_server(app, "127.0.0.1", 8050)



# callback!(
#     app,
#     Output("output-graph", "children"),
#     Output("yrsum-graph", "children"),
#     [Input("upload-data", "contents")],
#     [State("upload-data", "filename")]
# ) do contents, filenames
#     if contents !== nothing
#         graphs = []
#         for (content, filename) in zip(contents, filenames)
#             graph = html_div([
#                 #html_h4(filename),
#                 dcc_graph(
#                     id = filename,
#                     figure = parse_contents(content, filename) #, filename
#                 )
#             ])
#             push!(graphs, graph)
#             yrs = html_div([
#                 #html_h4(filename), #yrsum
#                 dcc_graph(
#                     id = filename,
#                     figure = parse_y(content, filename)
#                 )
#             ])
#             push!(graphs, yrs)
#         end
#         return graphs
#     end
# end