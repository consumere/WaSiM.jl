# 
function mvwasim2(;ta=pwd(),pt="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05",kw...) 
        
        println("\nmoves all wq, xml and log files to from
        $pt  to current pwd")
        
        println("target dir is $ta");
        #@vv "af"
    
        af = filter(x -> occursin(r"wq", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i));kw...) #force=false
            println(basename(i)," --> ", ta)
        end
        #@rg "wq"
    
        af = filter(x -> occursin(r"xml", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i));kw...)
            println(basename(i)," --> ", ta)
        end
        af = filter(x -> occursin("modell", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i));kw...)
            println(basename(i)," --> ", ta)
        end
    
    end
    
    """
    takes first two cols of df and plots r2 QQ
    """
    