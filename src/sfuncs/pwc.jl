# 
function pwc()
    if platform == "windows"
        pt = wslpath()
        println("$pt in clipboard...")
        pt |> clipboard
    else
        pt = pwd()
        pt |> clipboard
        println("$pt in clipoard...")
    end
end

"""
return clipboard(x)
"""
