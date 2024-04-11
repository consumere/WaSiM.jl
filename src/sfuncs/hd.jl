# 
function hd(x::DataFrame)
        if nrow(x)>20
            @info "headtail of df:"
            vcat(first(x,5),last(x,5))
        else
            first(x,5)
        end
    end

    """
    headtail of df using mapcols
    """
    