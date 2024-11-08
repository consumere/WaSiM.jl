# 
function pyplot_df(df::DataFrame;log=false)
    x = df.date
    ln = (filter(x -> !occursin(r"date|month|year", x), names(df)))

       
    for col in ln
        y = df[!, Symbol(col)]
        PyPlot.plot(x, y, label=col)
    end

    if log
        PyPlot.yscale("log")
    end

    PyPlot.xlabel("Date")
    PyPlot.ylabel("")
    PyPlot.legend()
    ti = only(values(DataFrames.metadata(df)))
    PyPlot.title(ti)
    PyPlot.grid(true)
end