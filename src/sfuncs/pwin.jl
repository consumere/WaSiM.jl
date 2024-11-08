# 
function pwin()
    if platform == "windows"
        pt = wslpath()
        println("$pt in clipboard...")
        pt |> clipboard
    elseif platform == "osx"
        pt = pwd()
        pt |> clipboard
        println("$pt in clipboard...")
    else
        wsl_cmd = `wslpath -m .`
        wsl_path = readchomp(pipeline(wsl_cmd))
        # Return the WSL path
        clipboard(wsl_path)
        println("$wsl_path in clipboard...")
        #return wsl_path
    end
end