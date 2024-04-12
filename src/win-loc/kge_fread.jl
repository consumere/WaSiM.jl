
if length(ARGS)<2 #|| directory_path == nothing
	println("need directory_path at arg1 and file ending at arg2 for valdat!")
    exit()
end


# directory_path = ARGS[1];
# ext = ARGS[2];

directory_path = pwd();
ext = ARGS[1];

#if directory_path == nothing
#if length(directory_path)<1 #|| directory_path == nothing

#if ext == nothing
if length(ext)<1
	#println("need directory_path at arg1 and file ending at arg2 for valdat!")
	#printstyled("need directory_path at arg1 and file snippet at arg2 for valdat!\n",color=:red,underline=true)
	printstyled("need file snippet at arg1 for valdat!\n",color=:red,underline=true)
	printstyled("try jlkge qout|sort -nr|tee kge.log \n",color=:yellow,bold=true)
    exit()
end

#println("eval on ",pwd(),"/$directory_path...")

printstyled("eval on $directory_path...\n",color=:green)


using DataFrames, CSV, Statistics

# function nse(predictions::Vector{Float64}, targets::Vector{Float64})
#     return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
# end

function kge2(observed, simulated)
    r = cor(observed, simulated)
    α = std(simulated) / std(observed)
    β = mean(simulated) / mean(observed)
    return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
end

function kge_read(path::AbstractString, ext::AbstractString)
    files = readdir(path)
    for file in files
        file_path = joinpath(path, file)
        if isfile(file_path) && 
            contains(file, ext) &&
            (!occursin(r"xml|qgk|mon|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|[0-9]|jpg|zip|tar",file))
            dd = CSV.read(file_path,
            DataFrame,
            missingstring="-9999",
            normalizenames=true,
            silencewarnings=true,
            limit=typemax(Int),
            types = Float64,
            delim="\t")
            observed  = dd[:,5]
            simulated = dd[:,6]
            kge_value = kge2(observed, simulated)
#            nse_value = nse(observed, simulated)
            nm = basename(file)
            pt = abspath(file_path)
            println(replace("KGE value is $kge_value on $nm in $pt", "\\"  => "/"))
#            printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
        elseif isdir(file_path)
            dfs_in_subdir = kge_read(file_path, ext)
        end
    end
end
# directory_path = "D:/Wasim/Goldbach/revision/v5"
# ext = "qout_gb"

kge_read(directory_path,ext)

# println(replace("GFG is a CS portal.", "CS" => "Computer Science"))
# println(replace("GeeksforGeeks is a CS portal.", "GeeksforGeeks" => "GFG"))
printstyled("\ndone!\n",color=:green,bold=true)