println("removes empty NetCDFs ...")
# using Glob
# Find all netCDF files in the current directory and its subdirectories
#ncf = filter(x -> occursin(r"\.nc$", x), readdir(join=true))
#ncf = glob("*.nc")
#occursin(Regex(suffix),file)

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


ncf = recursive_glob_sufx(pwd(),suffix="nc")

# xpattern = r"[A-z]"  # replace with the desired pattern
# ncf = filter(x -> occursin(xpattern, x) != nothing, ncf)



# Check if any files were found
if length(ncf) <= 0
    error("No matches found. Aborting...")
else
    println("Sample output of your match:")
    #println(sample(ncf, 5, replace=true))
    println(rand((ncf), 8))
end


println("now precompiling Statistics and Rasters ...")
using Statistics
using Rasters
using DelimitedFiles
import NCDatasets


# s=dirname(pwd())
# outname=split(s,"/")[end]*"_nclogjl.txt"

outname=basename(pwd())*"_nclogjl.txt"

header = ["filename" "variance" "standard deviation" "min" "max"]
#header = zip("filename","variance","standard deviation","min","max")
#header = ("filename","variance","standard deviation","min","max")
open(outname, "w") do io
    writedlm(io, header)
end

#name(a)
#ans=zip(string(name(a)),a_var,a_std,a_min,a_max)

# v=(name(a),a_var,a_std,a_min,a_max)
# ans=zip(map(string, v))
# #ans=fname,zip(a_var,a_std,a_min,a_max)
# open(outname, "a") do io
#     writedlm(io, ans)
# end



#ans=[fname,a_var,a_std,a_min,a_max]

#writedlm("file.txt", ans)
#writedlm("file.txt",header)

# open("tst.txt", "w") do io
#     writedlm(io, ans)
# end



# Process each file
open(outname,  "a") do io #append
for i in eachindex(ncf)
    a = read(Raster(ncf[i]))
    if sum(map(ismissing, a)) == length(a)
        rm(ncf[i])
        println(basename(ncf[i]), "is empty...removed! \n")
    else
        a_min = minimum(skipmissing(a))
        a_max = maximum(skipmissing(a))
        a_var = var(skipmissing(a))
        a_std = std(skipmissing(a))
        #a_cor=cor(a_var,a_std)
        #quantile(0:20, [0.1, 0.5, 0.9])
        #a_q=quantile(skipmissing(a), [.01,.33,.66,.99])
        #ans=zip(a_var,a_std,a_q,a_min,a_max)
        fname=basename(ncf[i])
        #ans=zip(fname,a_var,a_std,a_min,a_max)
        #v=[name(a),a_var,a_std,a_min,a_max]
        #ans=zip(map(string, v))
        #ans=fname,zip(a_var,a_std,a_min,a_max)
        ans=[string(Rasters.name(a)),a_var,a_std,a_min,a_max]
        ans=hcat(ans...)   ## write as a 1x5 matrix using horizontal concatenation
        if a_min == a_max
            rm(ncf[i])
            println("--> min==max of ",fname,"...removed!!!","\n")
        # if ismissing(a_min) || ismissing(a_max)
        #     println(ncf[i], " has no minmax! values...skipped!")
        else
            # Process the file
            # ...
            #stats=Dict("min"=>a_min,"max"=>a_max);
            println(fname, 
            #" has non-missing values...\n",
            #" var:",a_var," min:",a_min," max:",a_max)
            "\tmin:",a_min,"\tmax:",a_max)
            #describe(a))
            #writedlm(io, [header; ans])
            writedlm(io, ans)
        end
    end
end
   
end #end of io... 

# v=[name(a),a_var,a_std,a_min,a_max]
# v=[string(name(a)),a_var,a_std,a_min,a_max]
# #ans=zip(map(string, v))
# transpose(v)
# open("output_file.txt", "w") do out_io
#       writedlm(out_io, zip(v), ',')
#     end


#numbers = rand(5)
#writedlm("geek.txt", numbers)
#For example, two vectors x and y of the same length can be written as two columns of 
#tab-delimited text to f by either writedlm(f, [x y]) or by writedlm(f, zip(x, y)).
#ans=zip(("fnam","\tvar:",var(skipmissing(a)),"\tmin:",a_min,"\tmax:",a_max))


#no printf!
function dd()
    cwd = pwd()
    dirs = readdir(".")
    osize = 0
    cnt = 0
    out = []
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
         cnt += 1
     end
    end 
    #@printf("%-40s (%d files) %15.4f GB\n","$(cwd):",cnt,osize/1024^3);
    os=round(osize/1024^2,digits=3)
    printstyled(rpad("$(cwd): $os MB\n",20),lpad("($cnt files)\n",10),color=:green);
end 

function jdd()
    cwd = pwd()
    dirs = readdir(".")
    for dir in dirs
        if isdir(dir)
            osize = 0
            cnt = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    osize += stat(joinpath(root, file)).size
                    cnt += 1
                end
            end
        os=round(osize/1024^2,digits=3)
        printstyled(rpad("$(dir): ",40),lpad("$os MB", 10),lpad("($cnt files)\n",30),color=:green);
        end
    end
    dd()
end

# printstyled("folder sizes on julia...\n",color=:red)
# println("folder sizes on julia...:\n",pwd(),":")
# dd()
# n=repeat(" - -",10)
# #println("\n",n*" subfolders of ",basename(pwd())*n,"\n")
# printstyled(n*" subfolders of ",basename(pwd())*n,"\n",color=:yellow)

jdd()

println("All done!")

#read logfile in again:
# fn = outname;
# df = CSV.read(fn,DataFrame,delim="\t")


# #In general, the code I provided uses built-in Julia functions and packages that are 
# #designed for efficient file and data processing. The Glob package, for example, uses optimized algorithms to search for files in a directory and its subdirectories. The NCDatasets package provides fast and memory-efficient access to netCDF data.
# #Therefore, it's possible that the code I provided is faster than your approach, 
# #but it also depends on the specific factors mentioned above. If you have a specific benchmark or performance requirements, you can try running both approaches on your specific dataset and measure the execution time.
# using Glob
# using NCDatasets
# using Statistics

# # Find all netCDF files in the current directory and its subdirectories
# ncf = glob("*.nc")
# # Filter the files using a pattern
# #xpattern = ""  # replace with the desired pattern
# ncf = filter(x -> occursin(xpattern, x), ncf)

# # Check if any files were found
# if length(ncf) <= 0
#     error("No matches found. Aborting...")
# else
#     println("Sample output of your match:")
#     println(sample(ncf, 5, replace=true))
# end

# # Process each file
# for f in ncf
#     ds = Dataset(ncf[f])
#     vars = ds.group
#     if length(vars) <= 0
#         rm(ncf[f])
#         println(ncf[f], " empty...removed!")
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


# using BenchmarkTools
# using Glob
# using NCDatasets
# using Statistics

# # Define the two approaches
# function approach1()
#     # Your approach here
# end

# function approach2()
#     # The code I provided here
# end

# # Define the directory and pattern to search for
# dir = "path/to/directory"
# pattern = ".*\.nc"

# # Benchmark the two approaches
# println("Benchmarking approach 1...")
# approach1_bench = @benchmark approach1($dir, $pattern)
# println("Approach 1 took $(minimum(approach1_bench.times)/1e6) ms")

# println("Benchmarking approach 2...")
# approach2_bench = @benchmark approach2($dir, $pattern)
# println("Approach 2 took $(minimum(approach2_bench.times)/1e6) ms")


# using PackageCompiler
#using Statistics, Rasters, DelimitedFiles
# pt,nm=(pwd(),"rasters-win.so")
# create_sysimage([:Statistics,:Rasters,:DelimitedFiles],sysimage_path=joinpath(pt,nm))
# create_sysimage([:Statistics,:Rasters],sysimage_path=joinpath(pt,nm))
