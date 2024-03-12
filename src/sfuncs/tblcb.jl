# 
function tblcb(kd::DataFrame;kwargs...) 
        pretty_table(kd,header=uppercasefirst.(names(kd)); kwargs...)
        txt_str = sprint(io -> pretty_table(io, kd, header=uppercasefirst.(names(kd)); kwargs...))
        return clipboard(txt_str) 
    end

    """
    applies to quarterly data
    qrtr(pt::Union{String,DataFrame},
        fun=sum;agg=quarterofyear)
    example 
     qrtr(df;agg=year)|>dfp
     qrtr(df;fun=mean)|>dfp
     qrtr(skipyr(df);fun=mean)|>dfp
    """
    