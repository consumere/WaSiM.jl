# 
function fread(z::Union{Regex,String})
        if isa(z,Regex)
            v = filter(file -> occursin(z,file), readdir());
            z = v[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt|log"i),v)]|>first
        end 
        println("loading $z ...")   
        m = map(x->(x=>Int64),[:YY,:MM,:DD,:HH])
        ms = ["-9999","-9999.0","lin", "log", "--"]
        df = CSV.read(z,DataFrame;
            delim="\t", header=1, 
            normalizenames=true, 
            missingstring=ms,
            types=Dict(m),
            stripwhitespace=true)|>f->dropmissing(f,:YY)
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
        select!(df,Not(1:4))
        DataFrames.metadata!(df, "filename", z, style=:note)
    end


    readall = loadalldfs
    readmeteo = waread
    loaddf = waread

    