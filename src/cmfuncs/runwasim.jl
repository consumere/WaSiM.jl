# 
function runwasim(ctlfile)
        try
            exec = normpath("C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05/wasimvzo64.exe")
            run(`$exec $ctlfile`)
        catch e
            println("no valid input!")
        end
    end

    """
    usage: 
    
    nd = ctsum("thickness",infile)
    nd = ctsum("ksat",infile)
    nd = ctsum("Par_n",infile)
    nd = ctsum("theta_sat",infile)

    """
    