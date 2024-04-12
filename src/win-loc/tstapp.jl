#######das geht 
# using Dash
# using DataFrames, CSV, PlotlyJS, RDatasets

# iris = RDatasets.dataset("datasets", "iris")

# p1 = Plot(iris, x=:SepalLength, y=:SepalWidth, mode="markers", marker_size=8, group=:Species)

# app = dash()

# app.layout = html_div() do
#     html_h4("Iris Sepal Length vs Sepal Width"),
#     dcc_graph(
#         id = "example-graph-3",
#         figure = p1,
#     )
# end
# run_server(app, "0.0.0.0", debug=true)

using Dash
# using DashHtmlComponents
# using DashCoreComponents
using DataFrames
using CSV
#using Plots
using StatsPlots
using Base64
using Dates


app = dash()

app.layout = html_div() do
    html_h1("WaSiM timeseries data"),
    dcc_upload(
        id = "upload-data",
        children = html_div([
            "Drag and Drop or ",
            html_a("Select Files")
        ]),
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
        multiple = false
    ),
    html_div(id = "output-data-upload")
end

function dfp(df::DataFrame)
    ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
        @warn "No basename in metadata!"
        ti = raw""
    end
    if (any(x->occursin("year",x),names(df)))
        s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
        @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
    end
end

function parse_contents(contents, filename)
    content_type, content_string = split(contents, ",")
    decoded = base64decode(content_string)
    #df = CSV.File(IOBuffer(decoded)) |> DataFrame
    ms = ["-9999", "lin", "log", "--"]
    # df = CSV.File(IOBuffer(decoded); delim="\t", header=1, normalizenames=true, 
    #         missingstring=ms, types=Float64) |> DataFrame

    df = CSV.read(IOBuffer(decoded), DataFrame;
        delim="\t",header=1,
        silencewarnings=true,
        missingstring=ms,
        types=Float64)


    dropmissing!(df,1)
    dt2 = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
    select!(df, Not(1:4))
    df.date = dt2
    #plot(df[!,:date], df[!, 5], label = filename, xlabel = "Time", ylabel = "Value")
    
    ti = filename
    
    if (any(x->occursin("year",x),names(df)))
        s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
        @df df Plots.plot(:year,cols(s),legend = :topright, title=ti)
    else    
    s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
    @df df Plots.plot(:date,cols(s),legend = :topright, title=ti)
    end
end

callback!(app,
    Output("output-data-upload", "children"),
    Input("upload-data", "contents"),
    State("upload-data", "filename")
) do contents, filename
    if contents == nothing
        return nothing
    end

    return html_div() do
        html_h5(filename),
        dcc_graph(
            id = "plot-data-upload",
            figure = parse_contents(contents, filename)
        )
    end
end

run_server(app, debug=true)
