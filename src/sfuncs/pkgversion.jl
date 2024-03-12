# 
function pkgversion(pkg::String)
    #pak = Pkg.dependencies()
    # filter(x-> x.second.name == pkg, Pkg.dependencies()) |> 
    #     x -> first(x)[2].version
    #info = filter(x-> x.second.name == pkg, Pkg.dependencies())
    #
    try
        pkg = strip(pkg)
        info = filter(x-> occursin(pkg,x.second.name),
            #Regex(x.second.name,"i")), 
            Pkg.dependencies())
        nam = first(info)[2].name
        ver = first(info)[2].version
        #nam = map(x->first(x)[2].name,info)
        #ver = map(x->first(x)[2].version,info)
        println("$nam: $ver")    
    catch
        @error ("not found: $pkg")
        return
    end
    
end


"""
newer version with copy df and switched func positions
"""
