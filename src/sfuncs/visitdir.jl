# 
function visitdir(dir)
    push!(visited_dirs, dir)
    cd(dir)
    println("Visited directory: ", dir)
end

# 