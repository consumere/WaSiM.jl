#qglin

if length(ARGS) == 0
	println("need args! <upper bound!>...")
    exit()
end

# if length(ARGS) == 1
# 	lyr=1
# 	println("skipping to layer 1...")
# else
# 	lyr=parse(Float16,ARGS[1]);
# end

bound = parse(Float16,ARGS[1]);


function recursive_glob_prfx(rootdir=".", prefix="")
    results::Vector{String} = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg|jpg",filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    results
end

files = recursive_glob_prfx(pwd(),"qgk")

# function newf(dir::String, prefix="")
#     dirs = readdir(dir)
#     routes::Vector{String} = []
#     for directory in dirs
#         if isfile("$dir/" * directory) & (startswith(directory, prefix)) & (!occursin(r"txt|yrly|nc|png|svg|jpg",directory))
#             push!(routes, "$dir/$directory")
#         else
#             if ~(directory in routes)
#                 newread = dir * "/$directory"
#                 newrs = newf(newread)
#                 [push!(routes, r) for r in newrs]
#             end
#         end
#     end
#     routes
# end

# function route_from_dir(dir::String)
#     dirs = readdir(dir)
#     routes::Vector{String} = []
#     for directory in dirs
#         if isfile("$dir/" * directory)
#             push!(routes, "$dir/$directory")
#         else
#             if ~(directory in routes)
#                 newread = dir * "/$directory"
#                 newrs = route_from_dir(newread)
#                 [push!(routes, r) for r in newrs]
#             end
#         end
#     end
#     routes
# end

# files = filter(x->occursin("qgk",x),route_from_dir(pwd()))
    
function qglin(upper::AbstractFloat,files::Vector{String})
        #files = glob("qgk*", ".") # find all files starting with qgk in current directory
        for file in files
            open(file) do f # open each file
                for line in eachline(f) # read each line
                    if startswith(line, "LIN") || startswith(line, "lin")   # check if line starts with LIN or lin
                        fields = split(line, "\t")                          # split line by tabs
                        k = parse(Float64, fields[end])                     # parse last field as float
                        if k > upper && k < 1                               # check if k is in range
                            #println(join([file, line], "\t"))               # print file name and line separated by tab
                            #printstyled(join([basename(file), line], "\t"),"\n",color=:green)
                            printstyled("$file\n",color=:light_white)
                            #println("$file:")
                            printstyled("$line\n",color=:green)
                        end
                    end
                end
            end
        end
end

qglin(bound,files)

# This function takes an upper bound as an argument and prints the lines that match the criteria from the files that match the pattern. I hope this helps.😊  
    # Quelle: Unterhaltung mit Bing, 08/04/2023(1) . https://bing.com/search?q=translate+bash+to+julia Zugegriffen 08/04/2023.
    # (2) shell - Julia from bash script - Stack Overflow. https://stackoverflow.com/questions/48476159/julia-from-bash-script Zugegriffen 08/04/2023.
    # (3) How to translate bash to .bat - General Usage - Julia ... - JuliaLang. https://discourse.julialang.org/t/how-to-translate-bash-to-bat/81652 Zugegriffen 08/04/2023.
    # (4) How to translate bash to .bat - General Usage - Julia ... - JuliaLang. https://discourse.julialang.org/t/how-to-translate-bash-to-bat/81652?page=2 Zugegriffen 08/04/2023.

