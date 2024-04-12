#my Rasters implementation ...
# using PyCall # for calling Python functions
# #xpattern can also be ARGS ...
# xpattern=r"[A-z]"
# #cd("/mnt/d/Wasim/streu/out/v12")
# # get a list of files with .nc extension
# lf = pyimport("glob").glob("**/*.nc", recursive=true);
# # filter the files by xpattern
# lf = filter(x -> occursin(xpattern, x), lf)
using Rasters
using Statistics
using Glob
# Find all netCDF files in the current directory and its subdirectories
#lf = filter(x -> occursin(r"\.nc$", x), readdir(join=true))
lf = glob("*.nc")
#xpattern = "la"
# Filter the files using a pattern
#xpattern = ""  # replace with the desired pattern
xpattern = r"[A-z]"  # replace with the desired pattern
lf = filter(x -> occursin(xpattern, x) != nothing, lf)
#lf = lf[map(x -> occursin(xpattern, x) != nothing, lf)]

# Check if any files were found
if length(lf) <= 0
    error("No matches found. Aborting...")
else
    println("Sample output of your match:")
    println(sample(lf, 5, replace=true))
end

# Process each file
for i in 1:length(lf)
    a = Raster(lf[i])
    if sum(map(ismissing, a)) == length(a)
        rm(lf[i])
        println(lf[i], " empty...removed!")
    else
        a_min = minimum(skipmissing(a))
        a_max = maximum(skipmissing(a))
        if ismissing(a_min) || ismissing(a_max)
            println(lf[i], " has no minmax! values...skipped!")
        else
            # Process the file
            # ...
            #stats=Dict("min"=>a_min,"max"=>a_max);
            println(basename(lf[i]), " has non-missing values...\n",
            "min:",a_min," max:",a_max," var:",var(skipmissing(a)))
            #describe(a))

        end
    end
end

println("All done!")



# #In general, the code I provided uses built-in Julia functions and packages that are 
# #designed for efficient file and data processing. The Glob package, for example, uses optimized algorithms to search for files in a directory and its subdirectories. The NCDatasets package provides fast and memory-efficient access to netCDF data.
# #Therefore, it's possible that the code I provided is faster than your approach, 
# #but it also depends on the specific factors mentioned above. If you have a specific benchmark or performance requirements, you can try running both approaches on your specific dataset and measure the execution time.
# using Glob
# using NCDatasets
# using Statistics

# # Find all netCDF files in the current directory and its subdirectories
# lf = glob("*.nc")
# # Filter the files using a pattern
# #xpattern = ""  # replace with the desired pattern
# lf = filter(x -> occursin(xpattern, x), lf)

# # Check if any files were found
# if length(lf) <= 0
#     error("No matches found. Aborting...")
# else
#     println("Sample output of your match:")
#     println(sample(lf, 5, replace=true))
# end

# # Process each file
# for f in lf
#     ds = Dataset(lf[f])
#     vars = ds.group
#     if length(vars) <= 0
#         rm(lf[f])
#         println(lf[f], " empty...removed!")
#     else
#         a = ds[vars[1]][:]
#         a_min = minimum(skipmissing(a))
#         a_max = maximum(skipmissing(a))
#         if ismissing(a_min) || ismissing(a_max)
#             println(f, " has no non-missing values...skipped!")
#         else
#             # Process the file
#             # ...
#             println(f, " has non-missing values...")
#         end
#     end
#     close(ds)
# end

# println("All done!")










# ######
# # Find all netCDF files in the current directory and its subdirectories
# lf = filter(x -> occursin(r"\.nc$", x), readdir(".", join=true))
# # filter files that match a certain pattern
# #xpattern = "" # set your desired pattern here
# xpattern=r"[A-z]"
# lf = filter(x -> occursin(xpattern, x), lf)

# # Check if any files were found
# if length(lf) <= 0
#     error("No matches found. Aborting...")
# else
#     println("Sample output of your match:")
#     println(sample(lf, 5, replace=false))
# end

# using Rasters
# using Statistics

# # Process each file
# for i in eachindex(lf)
#     a = Raster(lf[i];missingval=-9999)
#     if sum(map(ismissing, a)) == length(a)
#         rm(lf[i])
#         println(lf[i], " empty...removed!")
#     else
#        println((minimum(a), maximum(a)))
#     end
# end

# println("All done!")



# using Rasters
# using Statistics

# # Find all netCDF files in the current directory and its subdirectories
# lf = filter(x -> match(r"\.nc$", x), readdir(join=true))

# # Filter the files using a pattern
# xpattern = ""  # replace with the desired pattern
# lf = filter(x -> match(xpattern, x) != nothing, lf)

# # Check if any files were found
# if length(lf) <= 0
#     error("No matches found. Aborting...")
# else
#     println("Sample output of your match:")
#     println(sample(lf, 5, replace=true))
# end

# # Process each file
# for i in 1:length(lf)
#     a = Raster(lf[i])
#     if sum(map(ismissing, a)) == length(a)
#         rm(lf[i])
#         println(lf[i], " empty...removed!")
#     else
#         a_min = minimum(skipmissing(a))
#         a_max = maximum(skipmissing(a))
#         if ismissing(a_min) || ismissing(a_max)
#             println(lf[i], " has no non-missing values...skipped!")
#         else
#             # Process the file
#             # ...
#             println(lf[i], " has non-missing values...")
#         end
#     end
# end

# println("All done!")






# using Rasters

# # find all files with extension ".nc"
# # for (root, dirs, files) in walkdir(dir)
# #     for file in files
# #         size += stat(joinpath(root, file)).size
# #     end
# # end
# # find all files with extension ".nc"
# lf = filter(x -> occursin(r"\.nc$", x), readdir(".", join=true))
# # filter files that match a certain pattern
# #xpattern = "" # set your desired pattern here
# xpattern=r"[A-z]"
# lf = filter(x -> occursin(xpattern, x), lf)
# # stop the program if no matching files are found
# if length(lf) <= 0
#     error("...no match found...abort")
# else
#     println("sample output of your match: ")
#     println(sample(lf, 5, replace=true))
# end

# # loop through each file and check if it's empty
# # for file in lf
# #     a = Raster(file)
# #     #setminmax!(a)
# #     if sum(minmax(a)) == 0
# #         rm(lf[i])
# #         println(lf[i], " empty...removed!")
# #     else
# #        println(file," checked...")
# #        println(minmax(a))
# #     end
# # end
# #for i in 1:length(lf)
# using Rasters
# using Statistics

# for i in eachindex(lf)
#     a = Raster(lf[i])
#     if sum(map(ismissing, a)) == length(a)
#     #if mapreduce(x -> ismissing(x) ? 0 : x, +, a) == 0 #geht auch
#         rm(lf[i])
#         println(lf[i], " empty...removed!")
#     else
#        println(minmax(a))
#     end
# end
# println("all done!")


# # # check if there are any matches
# # if isempty(lf)
# #     error("..no match found...abort")
# # else
# #     println("sample output of your match: ")
# #     println(rand(lf, 5)) # sample 5 files randomly with replacement
# # end

# # # load RASTERS
# # using Rasters

# # # loop over the files
# # for i in eachindex(lf)
# #     # read a file as a raster object using rast function from R
# #     #a = rcall(:rast, lf[i])
# #     a = read(Raster(lf[i]))   
# #     # set min and max values using setMinMax function from R
# #     #rcall(setMinMax, a)
    
# #     #Rasters.minimum(a)
# #     # get the min and max values as a Julia vector using rcopy function from RCall.jl
# #     describe(a)
# #     mm[3]


# #     # check if the sum of min and max values is zero
# #     if sum(mm) == 0
# #         # remove the file using rm function from Python's os module
# #         pyimport("os").rm(lf[i])
# #         println(lf[i], " empty...removed!")
# #     end 
# # end

# # println("all done!")


