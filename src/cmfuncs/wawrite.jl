# 
function wawrite(df::DataFrame,file::AbstractString)
        """
        newer version with copy df and switched func positions
        """
        dout = copy(df)
        #dout[!,Cols(r"date")]
        #in("date",names(dout))
        if in("year",names(dout))
            @warn "yearcol found!"
            CSV.write(file, dout, 
            transform = (col, val) -> something(val, missing), delim="\t")  
            return
        end
        dout.YY = map(x ->year(x),dout.date)
        dout.MM = map(x ->month(x),dout.date)
        dout.DD = map(x ->day(x),dout.date)
        dout[!, "HH"] .= 0
        dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
        CSV.write(file, dout, 
        transform = (col, val) -> something(val, missing), delim="\t")  
        nothing
    end

    