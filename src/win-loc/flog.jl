
@time using DataFrames
@time using CSV: File, write
#@time using DataFrames: DataFrame, ByRow, Cols, Not, permutedims, dropmissing, rename!, nrow, hcat

#jlno $jlpt/flog.jl
#julia --threads auto -q --startup-file=no $jlpt/flog.jl
#julia.exe --threads auto -q --startup-file=no $(wslpath -m $jlpt/flog.jl)


# function jlog(){ 
#     vers=$(whoami)
#     [[ "$vers" == ubu ]] && { 
#       printf "julia on $vers\n"
#       julia --startup-file=no --color=yes --threads auto -q --optimize=0 $jlpt/flog.jl $*;
#       return 1
#     }
#       julia.exe --startup-file=no --color=yes --threads auto -q $(wslpath -m $jlpt/flog.jl) $*;
#       return 1
#     }

"""
function rglob(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())
with filter for files that contain pattern and 
    not contain yrly, nc, jpg, png, svg    
"""
function rglob(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())
    #check it rootdir is a dir
    rootdir = isdir(rootdir) ? rootdir : throw(ArgumentError("rootdir is not a directory!"))    
    results = []
    if (any(x->isdir(x), readdir(rootdir)))
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for dir in dirs
                printstyled("check dir: $dir\n",color=:light_red)
            end
            for filename in filenames
                if pattern isa AbstractString
                    if (occursin(Regex("$(pattern)", "i"), filename) && (!occursin(r"yrly|nc|jpg|png|svg",filename)))
                        push!(results, joinpath(looproot, filename)) 
                    end
                else
                    if (occursin(pattern, filename) && (!occursin(r"yrly|nc|jpg|png|svg",filename)))
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
        end
    else
        printstyled("no dirs in $rootdir !\n", color=:light_red)
        for filename in (filter(x->isfile(x), readdir(rootdir; join=true)))
            if pattern isa AbstractString
                if (occursin(Regex("$(pattern)", "i"), filename))
                    push!(results, filename) 
                end
            else
                if (occursin(pattern, filename))
                    push!(results, filename) 
                end
            end
        end
    end
    return results
end

function qall(;recursive=false)
    if recursive
        files = rglob("qgk")
    else
        files = filter(x->occursin(r"qgk",x) & isfile(x), readdir())
    end
    outdf::Vector{DataFrame} = []
    for file in files
        x = file
        try
            df = DataFrame(File(x, header=1, 
                                delim="\t",
                                skipto=366,
                                ntasks=1,
                                types = String,
                                silencewarnings=true,
                                ignorerepeated=true,
                                ignoreemptyrows=true,
                                stripwhitespace=true))
            println(x)
            pattern = r"^[LIN. R]|^[LOG. R]|^CO"
            mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
            dx = df[mask, :]
            dx = permutedims(dx) |>dropmissing
            basins = []
            for i in names(df)[5:end]
                push!(basins,i)
            end
            insert!(basins, 1, "score")
            insert!(basins, 2, "timestep")
            dx[!, "basin"] = basins

            cn = (dx[1,:])
            rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
            dout = dx[3:end,:]
            for col in names(dout)[1:end-1]
                dout[!, col] = parse.(Float64, replace.(dout[!, col], "," => ""))
            end          
            dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
            dout.basin = parse.(Int64,dout[!, :basin])
            dout.nm .= replace(file,".\\"=>"")
            push!(outdf,dout)
        catch
            @warn("$x can not be loaded as a DataFrame! ")
            # Skip files that can't be loaded as a DataFrame
            continue
        end
    end
    return(outdf)
end

"""
find LOG. R-SQUARE > .4 recursivley
uses v = qall(;recursive=true)
"""
function findlog(;lb=.4)
    v = qall(;recursive=true)
    k = try 
        #map(x->subset(x,3 => ByRow(>(lb) && <(1.0))),v)
        map(x -> subset(x, :3 => ByRow(y -> y > lb && y < 1.0)), v)
        catch
            @error "no df on lowerbound! "
        return
    end
    
    if sum(map(nrow,k))==0
        @error "no values at lowerbound found! "
        return
    end

    df = try 
        reduce(vcat,k) 
        catch
            @error "vcat failed! "
        return
    end
   
    if nrow(df) <= 0
        @error "no values at lowerbound found!"
        return
    end

    # df = try 
    #     subset(df,3 => ByRow(<(1.0)))
    #     catch
    #         println(first(k))
    #         @error "no df on upperbound! "
    #     return
    # end
    
    
    return df
end


if length(ARGS) == 0
    bound = .4
    @info "no lower bound given, using default .4"
else
    bound = parse(Float64,ARGS[1])
end

#@time dfs = findlog(;lb=bound)
dfs = try
    findlog(;lb=bound)
    catch e
#    @error "$e \n findlog failed!"
    @error "findlog returns nothing..."
    nothing
end

#return 
#if nrow(dfs) <= 0 || isempty(dfs) || dfs isa Nothing
if dfs isa Nothing
#    println(dfs)
    printstyled("no values at lowerbound found! \n", color=:light_red)
    println("exiting...")
    exit()
else
	sort!(dfs,3;rev=false)
    println(dfs)
	dm = pwd()|>splitpath|>last
	outfile="logjl-$dm.txt"

	write(outfile, dfs, 
		transform = (col, val) -> something(val, missing), 
		delim="\t")  

    printstyled("written to $outfile ! \n", color=:light_green)    
end


