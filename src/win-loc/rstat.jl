#rstat

# Retrieve command line arguments
args = ARGS
# Get the file pattern argument
xpattern = args[1]
#xpattern = "alb"
# List files matching the pattern
xlf = filter(
    x -> occursin(Regex(xpattern, "i"), x),    
    readdir())

# filter(
#     x -> 
#     (occursin(regand(xpattern,"nc"), x),
#     readdir()))

# xlf = filter(
#     x -> 
#     (occursin(Regex(xpattern, "i"), x) &
#     endswith(x,".nc")    , 
#     readdir()))
    
#(!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|png|svg",x)

xlf = filter(x->endswith(x,".nc"),xlf)

# Check for valid input and print information
if isempty(xpattern) || ismissing(xpattern) || ismissing(xlf)
    println("<<Julia stats from cmdline>>")
    println("usage rst <filepattern>!")
    println("match failed on \n*", xpattern, "\n... abort.")
    exit(1)
elseif isempty(xlf)
    println(xpattern, " is not a valid input!\n... abort.")
    exit(1)
end

println("eval on: ", xlf, "...\n")

using Rasters
using StatsBase
# Define the function to read raster files
rfun(x) = begin
    try
        y = Raster(x)
#        println("..working on it...")
        return y
    catch
        error("$(basename(x)) is not a valid input!")
        exit(1)
    end
end

# Describe the first file
#des = Rasters.summary(xlf[1])


#else
    for x in eachindex(xlf)
        #println("eval on: ", xlf[x], "...\n")
        r = rfun(xlf[x])
        println("grid: ",xlf[x])
        println("dims: ", summary(r))
        println("dims: ", describe(r))    
    end   
#end


#     println("eval on: ", xlf[1], "...\n")
#     r = rfun(xlf[1])
#     println("grid: ", xlf[1])
#     println("dims: ", summary(r))
#     println("dims: ", describe(r))
#     #setminmax!(r)
#     #println("crs: ", crs(r, unit = false, asproj4 = true))
#     #println("resolution: ", res(r))
#     #println("ncells: ", ncell(r))
#     #summary(r, size = ncell(r))
# end


#summarystats(r)

