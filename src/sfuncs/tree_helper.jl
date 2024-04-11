# 
function tree_helper(root, indent)
    
    entries = filter(isdir, readdir(root))
    
    # #das geht auch, aber too much
    # entries = []
    # for (root, dirs, files) in walkdir(root)
    #     for dir in dirs
    #         #dstr = (joinpath(root, dir))
    #         dstr = dir
    #         push!(entries, dstr)
    #     end
    # end
    
    for (i, entry) in enumerate(entries)
        entry_path = joinpath(root, entry)  
        is_last = i == length(entries)
        entry_path = replace(entry_path, "/"=>"\\") #for win.
        
        entry_path = replace(entry_path, "\\mnt"=>"mnt") #
        # Calculate the depth of the current directory
        #depth = count(x -> x == '/', entry_path) - count(x -> x == '/', root)
        depth = count(x -> x == '\\', entry_path) - count(x -> x == '\\', root)
        
        # Print indentation
        if depth > 0
            print(indent)
            print(is_last ? "└── " : "├── ")
        end
        
        # Print the current entry name
        println(entry)
        
        # Recursively process subdirectories
        tree_helper(entry_path, indent * (is_last ? "    " : "│   "))
    end
end

