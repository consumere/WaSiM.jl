# 
function mvwasim2() 
        
        println("\nmoves all wq, xml and log files to from
        c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05
        to current pwd")
        ta=pwd()
        pt="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05";
        println("target dir is $ta");
        #@vv "af"
    
        af = filter(x -> occursin(r"wq", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i)))
            println(basename(i)," --> ", ta)
        end
        #@rg "wq"
    
        af = filter(x -> occursin(r"xml", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i)))
            println(basename(i)," --> ", ta)
        end
        af = filter(x -> occursin("modell", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i)))
            println(basename(i)," --> ", ta)
        end
    
    end

    