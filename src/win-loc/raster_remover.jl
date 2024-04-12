println("removes empty grd files recursively...")

function recursive_glob_sufx(rootdir="."; suffix::AbstractString)
    results = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (endswith(filename, Regex(suffix))) && (!occursin(r"txt|yrly|png|svg",filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    return results
end

# "d:/Wasim/Testdata/Output_ascii/oldrun"|>cd
#"D:/Wasim/Docker/clusterbrend/output-gw/"|>cd
#"D:/temp/dgm100_sinn/input_fromDockerimage/sinn_noponds_out/"|>cd
ncf = recursive_glob_sufx(pwd(),suffix="[.]grd")

#ncf = filter(x->!occursin(Regex("stack","i"),x),ncf)

# Check if any files were found
if length(ncf) <= 0
    error("No matches grd files found. Abort now...")
else
    println("Sample output of your match:")
    println(rand((ncf), 4))
end

println("now precompiling Statistics and ArchGDAL ...")
using ArchGDAL
const AG = ArchGDAL
using Statistics
using DelimitedFiles
# using Rasters
# import NCDatasets

outname=basename(pwd())*"_nclogjl.txt"

header = ["filename" "variance" "standard deviation" "min" "max"]

open(outname, "w") do io
    writedlm(io, header)
end

# Define the name of the text file to store the list of files to remove
remove_list_file = "files_to_remove.txt"

open(outname, "a") do io # append
    # Initialize an empty array to store the names of files to remove
    files_to_remove = String[]

    for i in eachindex(ncf)
        xxv = AG.readraster(ncf[i])
        if occursin(Regex("stack","i"),ncf[i])
            a = AG.getband(xxv, 2)
        else
            a = AG.getband(xxv, 1)
        end

        if sum(map(ismissing, a)) == length(a) || unique(a) == [0, -9999]
            xxv = nothing
            a = nothing
            #close(xxv)  # Close the file if it's not `nothing`
            push!(files_to_remove, ncf[i])  # Add the file name to the list
            println(basename(ncf[i]), " is empty...marked for removal! \n")
        else
            a_min = minimum(skipmissing(a))
            a_max = maximum(skipmissing(a))
            a_var = var(skipmissing(a))
            a_std = std(skipmissing(a))
            fname = basename(ncf[i])
            ans = [string(basename(ncf[i])), a_var, a_std, a_min, a_max]
            ans = hcat(ans...)  # write as a 1x5 matrix using horizontal concatenation
            if a_min == a_max
                xxv = nothing
                if xxv !== nothing  # Check if xxv is not `nothing` before closing
                    close(xxv)  # Close the file if it's not `nothing`
                end
                push!(files_to_remove, ncf[i])  # Add the file name to the list
                println("--> min == max of ", fname, "...marked for removal!!!", "\n")
            else
                println(fname, "\tmin:", a_min, "\tmax:", a_max)
                writedlm(io, ans)
            end
        end
    end

    # Write the list of files to remove to the separate text file
    writedlm(remove_list_file, files_to_remove)
end  # end of io...




function sdf()
    
    function calculate_folder_size(directory)
        size = 0
        count = 0
        for (root, dirs, files) in walkdir(directory)
            for file in files
                size += stat(joinpath(root, file)).size
                count += 1
            end
        end
        return size, count
    end

    function print_folder_size(directory, size, count)
        size_gb = round(size / 1024^3, digits=3)
        printstyled(rpad("$directory: $size_gb GB", 40), color=:green)
        printstyled(lpad("($count files)", 20), "\n", color=:green)
    end

    printstyled("folder sizes on Julia...\n", color=:red)

    cwd = pwd()
    dirs = readdir(cwd)

    rs, rcnt = calculate_folder_size(cwd)

    print_folder_size(cwd,rs,rcnt)

    n = repeat(" - -", 10)
    printstyled(n * " subfolders of " * basename(cwd) * n, "\n", color=:yellow)

    for dir in dirs
        if isdir(dir)
            size, count = calculate_folder_size(joinpath(cwd, dir))
            print_folder_size(dir, size, count)
        end
    end
end

function remove_files_from_list(file_list_filename::AbstractString)
    # Read the list of files to remove from the specified file
    files_to_remove = readlines(file_list_filename)

    # Loop through the list and remove each file
    for file_to_remove in files_to_remove
        try
            rm(file_to_remove, force = true)
            println("Removed file: ", file_to_remove)
        catch err
            println("Error while removing file: ", file_to_remove)
            println(err)
        end
    end
end

remove_files_from_list("files_to_remove.txt")

sdf()

println("All done!")



##errors on stack
#fn = raw"D:\Wasim\Testdata\Output_ascii\hgeocols500.grd"
#AG.readraster(fn)
# using ArchGDAL
# const AG = ArchGDAL
# fn = raw"D:\Wasim\Testdata\Output_ascii\evarcol500.grd"
# fn = raw"D:\Wasim\Testdata\Output_ascii\tempcol500.grd"
# fn = raw"D:\Wasim\regio\rcm200\v6\rcm.dhm"
# fn = raw"D:\Wasim\Testdata\Output_ascii\AnnualTemperature"
# fn = raw"D:\Wasim\Testdata\Output_ascii\gwn_col500.grd"
# r = AG.readraster(fn)
# a = AG.getband(r,1)
# AG.setnodatavalue!(a,Float32(-9999.0f0))
# AG.setnodatavalue!(a,-9999)
# #println(minimum(a),maximum(a))
# println(extrema(a))
# using Plots
# #plot(a)
# contour(a)
# contourf(a)

# dx = a
# msk = 610.0
# msk = -10.0
# bitmat = dx .> msk
# dx_filtered = dx .* bitmat
# dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
# dx_output .= NaN
# dx_output[bitmat] .= dx_filtered[bitmat]
# dx_output|>contourf

# #dx = dx_output


# open(outname,  "a") do io #append
#     for i in eachindex(ncf)
#         xxv = AG.readraster(
#             ncf[i]
#             )
#         a = AG.getband(xxv,1)
       
#         if sum(map(ismissing, a)) == length(a)|| unique(a) == [0,-9999]
#             xxv = nothing
#             a = nothing
#             rm(ncf[i],force=true)
#             println(basename(ncf[i]), "is empty...removed! \n")
#         else
#             a_min = minimum(skipmissing(a))
#             a_max = maximum(skipmissing(a))
#             a_var = var(skipmissing(a))
#             a_std = std(skipmissing(a))
#             fname=basename(ncf[i])
#             ans=[string(basename(ncf[i])),a_var,a_std,a_min,a_max]
#             ans=hcat(ans...)   ## write as a 1x5 matrix using horizontal concatenation
#             if a_min == a_max
#                 xxv = nothing
#                 # a = nothing
#                 rm(ncf[i],force=true)
#                 println("--> min==max of ",fname,"...removed!!!","\n")
#             else
#                 println(fname, 
#                 "\tmin:",a_min,"\tmax:",a_max)
#                 writedlm(io, ans)
#             end
#         end
#     end  
#     end #end of io... 



# open(outname,  "a") do io #append
#     for i in eachindex(ncf)
#         #a = read(Raster(ncf[i]))
#         #i=4
#         xxv = AG.readraster(
#             ncf[i]
#             )
#         # a = read(Raster(xxv[:,:,1],(X,Y);    
#         #     missingval=-9999.0))
#         #mappedcrs=25832,
#         a = AG.getband(xxv,1)
#         #a = AG.setnodatavalue!(a,Float32(-9999))
#         #a = AG.setnodatavalue!(a,Int32(-9999))
#         #allequal(a)
#         #unique(a)
#         #[0,-9999]|>typeof

        
#         #println(minimum(a),maximum(a))
#         #nm = ncf[i]
#         # println("extrema of $i :",extrema(a))
#         # msk = -3000
#         # bitmat = a .> msk
#         # dx_filtered = a .* bitmat
#         # dx_output = Matrix{Float32}(undef, size(a, 1),size(a,2))
#         # dx_output .= NaN
#         # dx_output[bitmat] .= dx_filtered[bitmat]
#         # a = dx_output
        
#         # #heatmap(a, c=:matter)
#         # allequal(a)
#         # unique(a)
#         # if sum(map(ismissing, a)) == length(a)|| unique(a) == [0,-9999]
#         #     println("YES")
#         # end
        
#         if sum(map(ismissing, a)) == length(a)|| unique(a) == [0,-9999]
#             xxv = nothing
#             a = nothing
#             rm(ncf[i],force=true)
#             println(basename(ncf[i]), "is empty...removed! \n")
#         else
#             a_min = minimum(skipmissing(a))
#             a_max = maximum(skipmissing(a))
#             a_var = var(skipmissing(a))
#             a_std = std(skipmissing(a))
#             fname=basename(ncf[i])
#             ans=[string(basename(ncf[i])),a_var,a_std,a_min,a_max]
#             ans=hcat(ans...)   ## write as a 1x5 matrix using horizontal concatenation
#             if a_min == a_max
#                 xxv = nothing
#                 # a = nothing
#                 rm(ncf[i],force=true)
#                 println("--> min==max of ",fname,"...removed!!!","\n")
#             else
#                 println(fname, 
#                 "\tmin:",a_min,"\tmax:",a_max)
#                 writedlm(io, ans)
#             end
#         end
#     end  
#     end #end of io... 
 


# import ArchGDAL as AG
# AG.read(ncf[3])
# @vv "AG.read"

# xxv = AG.readraster(
#     ncf[1]
#     )
# my = Raster(xxv[:,:,1],(X,Y);mappedcrs=25832,missingval=-9999.0)
# using Plots
# plot(my)
#plot(Rasters.trim(my))

#flags = Dict(:tr => [1000,1000], :r => :near)
#plot(warp(my, flags)) #Attempt to create 0x0 dataset is illegal,sizes must be larger than zero.  

# flags = Dict(:tr => [50,50], :r => :near)
# plot(warp(my, flags))       #works!

# size(my)
# flags = Dict(:tr => [.5,.5], :r => :near) #increase res
# myz = warp(my, flags)
# size(myz)
# plot(myz)

#reproject(target::GeoFormat, x)
#reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)
#Run `using ArchGDAL` to use `reproject`
#using ArchGDAL
#import ArchGDAL
#rpj = Rasters.reproject(EPSG(25832),EPSG(4326),my)
# Process each file

#rm("d:\\Wasim\\Testdata\\Output_ascii\\oldrun\\evarcol500.grd")

# dataset = AG.read(ncf[i])
# R = AG.imread(dataset,1)
# R = Raster(dataset[:,:,1])


