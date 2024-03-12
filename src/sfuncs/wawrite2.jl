# 
function wawrite2(df::DataFrame,file::AbstractString;
        hourint::Int=24,
        misval::Union{Number,Missing}=-9999)
    dout = copy(df)
    if in("year",names(dout))
        @warn "yearcol found!"
        CSV.write(file, dout, 
        transform = (col, val) -> something(val, missing), delim="\t")  
        return
    end
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= hourint
    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    CSV.write(file, dout, 
    transform = (col, val) -> something(val, misval), delim="\t")  
    nothing
end



"""
writes df to file, no date conversion
"""
