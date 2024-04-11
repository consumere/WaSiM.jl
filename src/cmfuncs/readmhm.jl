# 
function readmhm(x::AbstractString)
        "--- main reader ---"
        ms=["-9999","lin","log"]
        #x=lk
        df::DataFrame = CSV.read(x,DataFrame,
        header=false,
        missingstring = ms,
        types =  Dict(6=>Float64,
        1=>Int64,
        2=>Int64,
        3=>Int64,
        4=>Int64,
        5=>Int64), 
        skipto = 6,
        delim=" ",
        ignorerepeated=true,
        silencewarnings=false,
        normalizenames=true)
        DataFrames.metadata!(df, "filename", x, style=:note);
        #DataFrames.metadata!(df, "basename", basename(x), style=:note);
        #split(basename(x),".")[1]
        #    drop=(i, nm) -> i == 4) |> dropmissing
        df = df[!,Cols(1:3,end)]      #subset to dayres only
        newnames = ["YY","MM","DD",split(basename(x),".")[1]]
        rename!(df, newnames)
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,Not(1:3)]    
    end

    # 