# 
function pfix()
    #files = filter(x -> occursin(Regex(m,"i"), x), readdir())
    files = readdir()
    fzs = [(file, filesize(file) / 2^20) for file in files]
    tot = round(sum(map(x->x[2],fzs));digits=3)
    #d=Dict(fzs) 
    #sort(collect(d), by = x -> x[2];rev=true)
    sv = sort(collect(fzs), by = x -> x[2];rev=true)
    sv = first(sv,8)
    for i in sv
        printstyled(rpad(i[1],40),
        lpad(round(i[2];digits=4),10)," MB\n",
        color=:yellow)
    end
    #sv = first(Dict(sv),4)
    #printstyled("biggest files :\n $sv \n",color=:yellow)
    lf = length(files)
    printstyled("$lf files ...\n Total Size: $tot MB\n",color=:green)
end
    
"""
varinfo(Main,r"ds")
"""
