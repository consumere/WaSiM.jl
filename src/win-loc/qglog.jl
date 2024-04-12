#qglog

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
   
function qglog(upper::AbstractFloat,files::Vector{String})
        for file in files
            open(file) do f # open each file
                for line in eachline(f) # read each line
                    if startswith(line, "LOG") || startswith(line, "log")   # check if line starts with LIN or lin
                        fields = split(line, "\t")
                        if length(fields)>2
                            k = parse(Float64, fields[end])
                        else
                            k = 1
                        end
                        if k > upper && k < 1                               # check if k is in range
                            #printstyled("$file\n",color=:light_white)
                            println("$file:")
                            #line = strip(line,' ') 
                            line = replace(line,"\t"=>" ",". "=>".") 
                            printstyled("$line\n",color=:green)
                        end
                    end
                end
            end
        end
end

qglog(bound,files)
