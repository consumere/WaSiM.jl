# 
function wqlen()
        vw = glob("^wq")
        out::Vector{Float64} = []
        for w in vw
        m = Grep.grep("channel length",readlines(w))|>first
        m = split(m,":")|>last|>strip
        push!(out,parse.(Float64,String((m))))
        end
        out|>println
        A=round(sum(out)/1000,digits=2)
        println("total channel length is $A km")
    end

    