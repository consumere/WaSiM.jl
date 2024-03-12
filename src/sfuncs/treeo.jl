# 
function treeo(;root = pwd(), indent::String = "    ")
    println(root)
    for (root, dirs, files) in walkdir(root)
        # for file in files
        #     println(indent, "├── ", file)
        # end
        for d in dirs
            dp = joinpath(root, d)
            if Base.contains(dp, "\\")
                ln = splitpath(dp)|>length
                dp = replace(dp, root=>"","\\" => "──")
                println(indent,"├──",repeat("────",ln-1), dp)
            end
            #println(indent, "├── ", dp)
        end
    end
end

"""
5 biggest folders.
"""
