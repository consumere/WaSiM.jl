# 
function rmlat()
    files = filter(isfile, readdir(;sort=false))
    sorted_files = sort(files, by = mtime, rev = true)
    if !isempty(sorted_files)
        lat = sorted_files[1]
        println("This deletes the latest created file, i.e: ", lat)
        print("continue? (y/n): ")
        reply = readline(stdin)
        println(reply)
        if lowercase(string.(reply)) == "y"
            rm(lat, force=true)
            println("Deleted: ", lat)
        else
            @error "abort...."
        end
    end
end

