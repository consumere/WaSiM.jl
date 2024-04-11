# 
function old_waread(x::AbstractString)
        """
        skipping 6 lines - no dropmissing
            for Meteo Time Series
        """
        ms=["-9999","lin","log","-9999.0"]
        df = CSV.read(x,
        DataFrame,
        missingstring=ms,
        ntasks=8,
        skipto=6,
        limit=typemax(Int),
        delim="\t",
        silencewarnings=false,
        normalizenames=true,drop=(i, nm) -> i == 4) #|> dropmissing
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
        df=df[:,Not(1:3)]
        metadata!(df, "filename", x, style=:note);
    end

    