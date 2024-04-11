# 
function findctl(snippet::Union{String,Regex};recurse=false, dir="D:/Wasim/regio/control",suffix=".ctl")
        owd = abspath(pwd())
        platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
        
        if platform != "windows"
            wsl_cmd = `wslpath $(dir)`
            dir = readchomp(pipeline(wsl_cmd))
        end
        
        cd(dir)
        if recurse
            #files = rglob(Regex(suffix,"i"))
            res = []
            for (looproot, dirs, filenames) in walkdir(dir)
                for filename in filenames
                    if (occursin(snippet,filename) && endswith(filename, suffix))
                        push!(res, joinpath(looproot, filename)) 
                    end
                end
            end
        else
            files = filter(file -> endswith(file, suffix), readdir(dir,join=true))
            res = filter(x->occursin(snippet,x),files)
        end    
        
        cd(owd)
        return res
    end

    """
    sums up Conda.ROOTENV
    win.julia_conda:            3.06 GB
    """
    