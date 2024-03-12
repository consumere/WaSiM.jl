# 
function oplat()
    files = filter(isfile, readdir())
    sorted_files = sort(files, by = mtime, rev = true)
    if !isempty(sorted_files)
        lat = sorted_files[1]
        println("opening: $lat ...")
        if platform == "osx"
            run(`open $lat`)
        else
            try
                run(`cmd /c start $lat`)    
            catch
                @error "could not open $lat via cmd.exe"
            end
        end
    end
end   

"""
removes latest file
"""
