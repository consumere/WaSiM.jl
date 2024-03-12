# 
function wslp(winpt::AbstractString)
        """
        """
        #printstyled("needRAW like raw",color=:red)
        #return(raw$winpt)
        #quote($winpt);end
        #winpt=raw"D:\Wasim\regio\out\lowres\c5\loc5"
        winpt = replace(winpt, "\\"=> "/")
        #winpt = join(winpt,"'")
        #print("\"",winpt,"\"")
        #join(wsl_cmd,winpt)
        wsl_cmd = `wsl wslpath -m $winpt`
        run(wsl_cmd,winpt)
        #wsl_path = readchomp(pipeline(wsl_cmd))
        # Return the WSL path
        return wsl_path
    end

    