# 
function wslpath()
    # Run the `wslpath` command to convert the current directory to a WSL path
    #wsl_cmd = `wsl wslpath -a $(pwd)`
    wsl_cmd = `wsl wslpath -a .`
    wsl_path = readchomp(pipeline(wsl_cmd))
    # Return the WSL path
    return wsl_path
end

"""
import InteractiveUtils.clipboard
wslpath()|>clipboard
cb(wslpath())
"""
