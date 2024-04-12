#using DataFrames,CSV, FileIO
#using FileIO

function recursive_glob_prfx(rootdir=".", prefix="")
    results::Vector{String} = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            #if (startswith(filename, prefix)) && (!endswith(filename, "png")) && (!endswith(filename, "yrly"))      #NUR so gings.. :(
            if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    return results
end

# paths = []
# push!(paths, dirs) 
# return paths
# println(paths)
#if startswith(filename,"q") and !endswith("png")
#    show(filename) 
# results=[]
# if !endswith(filename, "png")
#     push!(results, filename)
# end

#and !endswith(filename,'png')
files = recursive_glob_prfx(pwd(),"qgk")
#for f in files;println(f);end 
#println(files < startswith("q") ? "less than" : "not less than")

#or
#recursive_glob_prfx(pwd(),"qgk")
#for z in files;println(z);end



for z in files
        println("generating\t",basename(z),"...")
        m=(filter(line -> occursin(r"^[LIN. R]|^[LOG. R]|^CO",line),readlines(open(z))))
        #for l in m;show(replace(l,"\t"=>" "));end
        for l in m
#            x=l
            x=replace(l,r"\s+"=>"\t")
            x=replace(x,".\t"=>" ")
#            x=replace(l,"\t"=>" ")
#            x=x*"\n" #concat strings
#            x=replace(x,'/"/'=>'\n')
#            show(x)
            println(x)
#            println(x,"\n")
        end
        #println(replace(m, "\t" => " "))
end

#6 ∉ a



#[replace(i,"\t"," ") for i in m]
#println(replace("GFG is a CS portal.", "CS" => "Computer Science"))



# function qall(files)
#     for file in files
#         try
# #            df = readtable(file, header=false)
#             df = DataFrame(CSV.File(file, header=false,delim="\t",missingstring="-9999"))
#             pattern = r"^[LIN. R]|^[LOG. R]|^CO"
#             mask = [occursin(regex(pattern), row[1]) for row in eachrow(df)]
#             dd = df[mask, :]
#             new = convert(Array{String, 1}, names(dd)[4:end])
#             insert!(new, 1, "timestep")
#             insert!(new, 1, "score")
#             dd[!, :basin] = new
#             dd = rename(dd, :basin => :index, :1 => :score, :2 => :timestep)
#             dd = dd[:, [:score, :timestep, 3:end]]
#             dd = dropna!(dd)
#             println(dd)
#         catch e
#             continue
#         end
#     end
# end

# #qall(files)
