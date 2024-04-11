# 
function waread2(x::String)
    ms = ["-9999", "lin", "log", "--"]
    df = CSV.File(x; delim="\t", header=1, normalizenames=true, missingstring=ms, types=Float64) |> DataFrame
    dropmissing!(df,1)
    dt2 = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
    select!(df, Not(1:4))
    df.date = dt2
    metadata!(df, "filename", x, style=:note)
    return df
end

"""
basic tsv reader, takes arguments from CSV.File
df = CSV.File(x;kw...)|>DataFrame|>z->dropmissing(z,1)
"""
