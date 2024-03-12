# 
function npplat(;opener="c:/Program Files (x86)/Notepad++/notepad++.exe")
    fl = lat()
    try
        if platform == "windows"
            run(`$opener $fl`)
        elseif platform == "unix" || platform == "linux" && src_path == "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
            fl = fl
            run(`"/mnt/c/Program Files (x86)/Notepad++/notepad++.exe" $fl`)
        else
            fx = `wslpath -m $fl` #translate from a WSL path to a Windows path, with '/' 
            opcmd = `wslpath -a $opener` #force result to absolute path format
            wslpt = readchomp(pipeline(opcmd))
            run(`$wslpt $fx`)

        end
    catch        
          @error "could not open $fl via notepad++
              check if notepad++ is installed in $opener"
    end
end

"""
df to clipboard
"""
