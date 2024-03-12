# 
function jlcnt(path=pwd(), level=0)
    n, s = count_files(path, 0)    
    printstyled("Directory: $path\n",color=:yellow)
    printstyled("Number of files: $n\n",color=:yellow)
    printstyled("Total size: $(round(s / (1024 * 1024), digits=2)) MB\n",color=:yellow)
end

"""
wslpath to windows
"""    
