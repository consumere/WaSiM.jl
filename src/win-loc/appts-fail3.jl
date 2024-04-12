using Dash
using PlotlyJS
using Dates
using DataFrames
using DataFramesMeta
using Base64
using CSV

app = Dash.dash()

function parse_contents(contents, filename)
    # Read the contents of the uploaded file
    content_type, content_string = split(contents, ',')
    
    # Decode the file content
    decoded = base64decode(content_string)
    df = CSV.read(IOBuffer(decoded), DataFrame; header=false)
    
    ms = ["-9999.0", "-9999", "lin", "log", "--"]
    df = filter(row -> !(occursin(r"[A-z]|-", row[1])), df)
    df = DataFrame(replace!(df, ms .=> missing))
    df = DataFrame(Meta.parse.(string.(df)), DataFrame(df))
    df = select!(df, Not(1:3))
    df.date = [Date(row[1], row[2], row[3]) for row in eachrow(df)]
    df[!, :date] = DateTime.(df.date)
    df = df[:, 4:end]
    df.date = Date(df.date)
    
    fig = PlotlyJS.Figure()
    ti = " " #placeholder
    for i in 1:size(df, 2)-1
        PlotlyJS.add_trace!(fig, Scatter(x=df.date, y=df[!, i], name=df.columns[i]))
    end
    
    fact, logy = 0.70, false
    if logy
        PlotlyJS.update_layout!(fig, yaxis_type="log", template="simple_white", title=ti,
                                xaxis_title="modeled time", yaxis_title="[unit/day]",
                                font_size=14, legend=dict(x=1.02, y=0.5))
    else
        PlotlyJS.update_layout!(fig, template="simple_white", title=ti,
                                xaxis_title="modeled time", yaxis_title="[unit/day]",
                                font_size=14, legend=dict(x=1.02, y=0.5))
    end
    
    return fig
end

function update_graph(contents, filenames)
    if contents !== nothing
        graphs = []
        for (content, filename) in zip(contents, filenames)
            graph = html_div(
                html_h4(filename),
                dcc.Graph(
                    id=filename,
                    figure=parse_contents(content, filename)
                )
            )
            push!(graphs, graph)
        end
        return graphs
    end
end

callback!(app, "output-graph", ("upload-data", "contents"), ("upload-data", "filename")) do contents, filenames
    return update_graph(contents, filenames)
end

#run_server(app, debug=true)

run_server(app, "127.0.0.1", 8050)
