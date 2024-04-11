# 
function writewa(file::AbstractString, df::DataFrame)
        dout = copy(df)
        dout.YY = map(x ->year(x),dout.date)
        dout.MM = map(x ->month(x),dout.date)
        dout.DD = map(x ->day(x),dout.date)
        dout[!, "HH"] .= 0
        #df = select!(df,Symbol.(filter(x->!occursin("date",x), names(df))))
        #dout = select(df, Not(:date))
        #dout = dout[!,Cols([:YY,:MM,:HH,:DD],1:end-4)]
        dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
        #cls = propertynames(df)|>sort|>reverse
        #df = df[!,cls[2:end]] 
        CSV.write(file, dout, transform = (col, val) -> something(val, missing),delim="\t")  
        nothing
    end

    """
    newer version with copy df and switched func positions
    """
    