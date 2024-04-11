# 
function bargroup(x::Union{Regex,String,DataFrame};leg = :topright,fun=monsum, lcols=1)
        if isa(x,DataFrame)
            df = (x)
        else
            df = waread(x)
        end

        v = map(
            (x->occursin(r"date", x) & !occursin(r"year|month", x)),
            (names(df))
            )

        if fun==false
            df = hcat(df[:,Cols(r"date")],df[!,Not(Cols(r"date"))])
            @info "no 