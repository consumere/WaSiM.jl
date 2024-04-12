# if length(ARGS) == 0
# 	println("need args! <file>...")
#     exit()
# end
#using Printf

#println("performing rmeq_rec from current dir...")
printstyled("performing rmeq_rec from $(pwd())...\n", color=:yellow, bold=true, underline=true)

#using DataFrames, CSV
import DataFrames: DataFrame, dropmissing!, nrow, ncol
import CSV: File


# function rmeq_rec(; rootdir = ".")
#     """
#     removes empty TS recursively; 
#     use with caution!
#     """
    
#     ext_regex = r".R|.py|.jl|.tex|.pl|.sh|.csv|.html|.xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt"i
#     ms = ["-9999", "lin", "log", "--"]
    
#     files::Vector{String} = []       #String[]
#     for (looproot, dirs, filenames) in walkdir(rootdir)
#         for filename in filenames
#             if !occursin(ext_regex, filename)
#                 push!(files, joinpath(looproot, filename))
#             end
#         end
#     end
    
#     for inF in files
#         if isfile(inF)
#             println("reading ", inF, "...")
#             df = CSV.File(inF; delim = "\t", header = 1,
#                           silencewarnings = true, 
#                           normalizenames = false, 
#                           missingstring = ms, 
#                           types = Float64) |> DataFrame
#             if (isempty(df) || nrow(dropmissing(df,ncol(df)))==0 || ncol(df)==0) 
#                 rm(inF)
#                 println(basename(inF), " removed!")
#             end
#         end
#     end
# end

function rmopt(; rootdir = ".")
    # Get the list of files
    #files = readdir()
    # files = filter(file -> 
    # !occursin(r"(wq_|pl$|sh$|fzt|ftz|log$|ini|otherdata|intern)", file)
    # , files)

    
    #ext_regex = r".R|.py|.jl|.tex|.pl|.wa|.sh|.csv|.html|.xml|fzt|ftz|png|svg|txt"i
    ext_regex = r"(wq_|pl$|sh$|fzt|ftz|log$|ini|otherdata|intern)"i
    ms = ["-9999", "lin", "log", "--"]
    
    files::Vector{String} = []       #String[]
    for (looproot, dirs, filenames) in walkdir(rootdir)
    for d in dirs
    println("reading $d ...")
    end
        for filename in filenames
            if !occursin(ext_regex, filename)
                push!(files, joinpath(looproot, filename))
            end
        end
    end

    #!occursin(r"pl|sh|csv|html|xml|fzt|ftz|log|ini|^wq|yrly|nc|tif|jpeg|png|svg|txt", file)
    # py|R|ftz_0|tex
    # Define the file extensions to keep
    keep_exts = [".ipynb", ".py", ".R", ".Rmd", ".log",
        ".tif", ".jpeg", ".png", ".svg",
        ".cpg", ".shx", ".dbf", ".prj", ".shp", ".tex", 
        ".csv", 
        ".html", ".ftz", ".ftz_0", ".txt", 
        ".list", ".nc", ".xml", ".sh", ".grd", ".yrly"]
        
    
        # Process each file
    for file in files
        # Check if the file extension is in the list of extensions to keep
        if !(splitext(file)[2] in keep_exts)
            try
            # Read the file as a DataFrame
            println("Processing $file ...")
            #ms = ["-9999", "lin", "log", "--"]
            ms = ["lin", "log", "--"] #!!!
            #df = CSV.read(file, DataFrame)
            #df = CSV.File(file; 
            df = File(file; 
                delim="\t", 
                header=1,
                silencewarnings=true, 
                normalizenames=false, 
                missingstring=ms, 
                types=Float64) |> DataFrame
            dropmissing!(df,ncol(df))
            # dropmissing!(df,4)
            # Calculate the sums of the 5th column and the last column
            y = sum(df[:, 5])
            x = sum(df[:, end])

            # Check the conditions
            if x <= (nrow(df) - 5) * -9999 && y <= (nrow(df) - 5) * -9999 || x == 0 && y == 0
                # Remove the file
                rm(file)
                printstyled("$file removed ...\n", color=:red)
            end
            catch e
                println("Error processing $file: $e")
            end
        end
    end
end


# function rmeq()
    # """
    # removes empty TS; 
    # use with caution!
    # """    
    # files = filter(file -> 
    # !occursin(r".tex|.pl$|.sh$|.csv|.html|.xml|.fzt|.ftz|.log|.ini|^wq|.yrly|.nc|.png|.svg|.txt", file)
    # , readdir())
    
    # ms = ["-9999", "lin", "log", "--"]
    # for inF in files
        # if isfile(inF)
            # df = CSV.File(inF; delim="\t", header=1,
            # silencewarnings=true, 
                # normalizenames=false, 
                # missingstring=ms, 
                # types=Float64) |> DataFrame
            # dropmissing!(df,ncol(df))
            # if nrow(df)==0
                # println(basename(inF)," removed!")
                # rm(inF)
            # end
        # end
    # end
# end

# function tff3(x::Vector{String})
#     for filename in x
#         if (
#             (!occursin(r"html|sh|txt|xml|qgk|fzt|ftz|log|ini|wq|yrly|nc|png|svg",
#                 filename))
#             )
        
#             b = Dict{String, Float64}()
#             m = Dict{String, Float64}()
#             h = Dict{String, Float64}()
#             cnte = Dict{String, Int64}()
#             ncols = 0

#             open(filename) do file
#                 first_line = true
                

#                 for (i, line) in enumerate(eachline(file))
#                     fields = split(line)
#                     if first_line
#                         #println("filename: $filename")
#                         printstyled("filename: $filename\n",color=:yellow)
#                         ncols = length(fields)
#                         datcols = ncols - 4
#                         println("no of fields: $datcols")
#                         println("year\t$(join(fields[5:end], "\t"))")
#                         first_line = false
#                         continue
#                     end

#                     if match(r"^\d{4}$", fields[1]) != nothing
#                         cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
#                         b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
#                         m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
#                         h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
#                     end
#                 end
#             end

#             if ncols <= 4
#                 for key in sort(collect(keys(b)))
#                     println("$key\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
#                 end
#             elseif ncols == 5
#                 for key in sort(collect(keys(b)))
#                     println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
#                 end
#             elseif ncols >= 6
#                 for key in sort(collect(keys(b)))
#                     println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", m[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", m[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
#                 end
#             end
#         end
#     end
# end

#rmeq_rec()
rmopt()


# fullpath = pwd()
# nms::Vector{String} = []
# #ARGS="qg"
# xm = ARGS[1]
# println("looking for match of *$xm* ...")

# for (subdir, _, filenames) in walkdir(fullpath)
#     for filename in filenames
#         if (occursin(Regex(xm,"i"),filename) && (!occursin(r"sh|yrly|nc|png|svg",filename)) )
#         pt = joinpath(subdir, filename)
#         #sizes[fullpath] = stat(fullpath).size
#         push!(nms,pt)
#         end
#     end
# end

# ti= try 
#     #tff3(nms)
#     rmeq(nms)
# catch 
#     wd = pwd();
#     @error "
#     no recursive match for $ARGS on $wd !
#     ...exiting now"; 
#     exit()
# end

#in wsl:
#julia --startup-file=no $(wslpath "C:\Users\Public\Documents\Python_Scripts\julia\tffrec.jl") rg*
#julia --startup-file=no $(wslpath "C:\Users\Public\Documents\Python_Scripts\julia\rmeq.jl")

#julia --startup-file=no --optimize=0 --threads 8  "C:\Users\Public\Documents\Python_Scripts\julia\rmeq.jl"

