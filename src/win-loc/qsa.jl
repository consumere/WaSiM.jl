#prints routingtable to console
#alias jlq

using DataFrames
#using CSV: read
import CSV
import PrettyTables: pretty_table
#import DelimitedFiles
# @time import Glob
# function recursive_glob_prfx(rootdir=".", prefix="")
#     return Glob.glob(prefix*"*", rootdir)
# end

"""
writes df to file, no date conversion
"""
function writedf(df::Union{DataFrame,String},file::Union{DataFrame,String})
    if df isa String
        @info "write DataFrame to $df !"
        CSV.write(df, file,
        transform = (col, val) -> something(val, missing),
            delim="\t")
        return
    end
    @info "write DataFrame to $file !"
    CSV.write(file, df, 
    transform = (col, val) -> something(val, missing),
        delim="\t")
    nothing
end

function qgkgrep(prefix::AbstractString)
    rootdir=pwd();
    results = []
    if (any(x->isdir(x),readdir()))
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                if (occursin(Regex(prefix,"i"),filename) && !occursin(r"txt|yrly|nc|png|svg",filename))
                    push!(results, joinpath(looproot, filename)) 
                end
            end
        end
    else
        #printstyled("no subdirs in $rootdir !\n",color=:light_red)
        printstyled("eval on: $rootdir !\n",color=:light_red)
        for filename in (filter(x->isfile(x),readdir(;join=false)))
            if (occursin(Regex(prefix,"i"),filename) && !occursin(r"txt|yrly|nc|png|svg",filename))
                printstyled("collecting $filename...\n",color=:light_yellow)
                push!(results, filename) 
            end
        end
    end
    return results
end

#files = recursive_glob_prfx(pwd(), "qgk")

# files = rglob("qgk")
# #if isfile(file_path) && (!occursin(r"xml|fzt|ftz|log|ini|wq|txt|yrly|nc|tif|jpg|png|svg",file))
# filter!(x -> !occursin(r"xml|fzt|ftz|log|ini|wq|txt|yrly|nc|tif|jpg|png|svg", x), files)

files = qgkgrep("qgk")
if isempty(files)
    printstyled("no qgk files found!\n",color=:light_red)
    exit()
end

function qall(files)
    for file in files
        x = file
        try
            df = CSV.read(x, DataFrame,delim="\t";
                          header=true,
                          types=String,
                          silencewarnings=true,
                          skipto=364)
            # df = DelimitedFiles.readdlm(x, '\t', String, '\n';
            #                         header=false,skipstart=364)
            # df = DataFrame(df,:auto)
            # printstyled("$x\n", bold=true, color=:light_blue)
            pattern = r"^[LIN. R]|^[LOG. R]|^CO"
            mask = occursin.(pattern, df[!, 1])
            ddd = (df[mask, :])
            
            new = names(ddd)[5:end]
            insert!(new, 1, "basin")
            insert!(new, 2, "timestep")
            #ddd = permutedims(select(ddd,Not(2)))
            ddd = permutedims(ddd)
            dropmissing!(ddd)
            ddd.basin = new
            #rename!(ddd, Symbol("x1") => :basin)
            #setindex!(ddd, :basin)
            select!(ddd, :basin, :)
            #ddd = ddd[2:end, :]
            #columns = vec(view(ddd, 1, :))
            df = permutedims(ddd)
            columns = uppercasefirst.(getproperty(df,propertynames(df)[1]))
            rename!(ddd, Symbol.(columns))
            kd = ddd[3:end, :]
            #println(ddd)
            #pretty_table(kd,header=uppercasefirst.(names(kd));)
            return kd
        catch e
            println(e)
            @warn "skipping $x"
        end
    end
end

dfs  = qall(files)

# function wread(x::String;skip=3)
#     df = DelimitedFiles.readdlm(x, '\t', Float64, '\n';
#         header=false,skipstart=skip)
#     df = DataFrame(df,:auto)
#     for i in 5:size(df,2)
#         df[!,i]=replace(df[!,i],-9999.0 => missing)
#     end 
#     for i in 5:size(df,2)
#         replace!(df[!,i],-9999.0 => missing)
#     end
#     for i in 1:3
#         df[!,i]=map(x ->Int(x),df[!,i])
#     end
#     #and parse dates...
#     df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
#     df=df[:,Not(1:4)]
#     metadata!(df, "filename", x, style=:note);
# end

"""
windows path to wsl
"""
function towsl(file_path::Union{String, Nothing})
    if isnothing(file_path)
        #return error("File path cannot be nothing")
        drive = lowercase(pwd()[1])
        rst = replace(pwd()[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    else
        drive = lowercase(file_path[1])
        rst = replace(file_path[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    end
    return wsl_path
end


"""
non-recursively search for control file in current directory
"""
function ctlx()
    matches::Vector{Any} = []
    for file in readdir(".")
        if endswith(file, ".xml")
            path = joinpath(pwd(), file)
            open(path) do f
                for line in eachline(f)
                    if occursin("compiling symbols in control file ", line)
                        fields = split(line)[8:end]
                        println(join(fields, " "))
                        out = join(fields, " ")
                        push!(matches, out)
                    end
                    if occursin("looking for starting date in ", line)
                        fields = split(line)[9:end]
                        str = replace(join(fields, " "),r"\"/>.*" => " -> ")
                        printstyled("discharge file is: $str\n",color=:light_green)
                    end
                end
            end
        end
    end
    if !isempty(matches)
        fl = first(matches)
        fl = split(fl) |> last
        fl = split(fl, "\"") |> first
        if (!occursin("regio",fl) && occursin("regio",pwd()))
            fl = replace(fl,"control"=>"D:/Wasim/regio/control")
        elseif (!occursin("brend",fl) && occursin("brend",pwd()))
            fl = replace(fl,"control"=>"D:/Wasim/Tanalys/DEM/brend_fab/control")            
        elseif (!occursin("temp",fl) && occursin("saale",pwd()))
            fl = replace(fl,"control"=>"D:/temp/saale/control")
        end
        return string(fl)
    else
        @warn "no control file or xml found!..."
        return ""
    end
end

# z = ctlx()

function routeg(input_file::String, output_file::String)
    open(output_file, "w") do output
        line_num = 0
        in_range = false

        for line in eachline(input_file)
            line_num += 1

            if line_num > 50 && contains(line, "routing_model")
                in_range = true
                line = replace(line, r"ß" => "ss")
                # line = replace(line, r"[\/]" => "_")
                # line = replace(line, r"_" => "-")
                line = replace(line, r"[,,]" => "")
                line = replace(line, r"\xc4" => "Ae")
                line = replace(line, r"\xd6" => "Oe")
                line = replace(line, r"\xdc" => "Ue")
                line = replace(line, r"\xe4" => "ae")
                line = replace(line, r"\xf6" => "oe")
                line = replace(line, r"\xfc" => "ue")
                line = replace(line, r"\xdf" => "ss")
                println(output, line)
            elseif in_range && contains(line, "timeoffset")
                in_range = false
                line = replace(line, r"ß" => "ss")
                # line = replace(line, r"[\/]" => "_")
                # line = replace(line, r"_" => "-")
                line = replace(line, r"[,,]" => "")
                line = replace(line, r"\xc4" => "Ae")
                line = replace(line, r"\xd6" => "Oe")
                line = replace(line, r"\xdc" => "Ue")
                line = replace(line, r"\xe4" => "ae")
                line = replace(line, r"\xf6" => "oe")
                line = replace(line, r"\xfc" => "ue")
                line = replace(line, r"\xdf" => "ss")
                println(output, line)
            elseif in_range
                line = replace(line, r"ß" => "ss")
                # line = replace(line, r"[\/]" => "_")
                # line = replace(line, r"_" => "-")
                line = replace(line, r"[,,]" => "")
                line = replace(line, r"\xc4" => "Ae")
                line = replace(line, r"\xd6" => "Oe")
                line = replace(line, r"\xdc" => "Ue")
                line = replace(line, r"\xe4" => "ae")
                line = replace(line, r"\xf6" => "oe")
                line = replace(line, r"\xfc" => "ue")
                line = replace(line, r"\xdf" => "ss")
                println(output, line)
            end
        end
    end
end

"""
dfroute(;ofl="route.txt")
reads from routeg(infile, ofl) and returns a DataFrame with the following columns:
"""
function dfroute(;ofl="routejl.txt")
    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"Basin",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"_>.*" => "",r"-" => "_"),df[:,3])
    sort!(df, :Basin)
    return df    
end

# function newctl(input_file::String)
#     df = DataFrame(sim=String[], obs=String[], name=String[])
#     line_num = 0
#     in_range = false

#     for line in eachline(input_file)
#         line_num += 1

#         if line_num > 50 && contains(line, "routing_model")
#             in_range = true
#             line = replace(line, r"ß" => "ss")
#             line = replace(line, r"[,,]" => "")
#             line = replace(line, r"\xc4" => "Ae")
#             line = replace(line, r"\xd6" => "Oe")
#             line = replace(line, r"\xdc" => "Ue")
#             line = replace(line, r"\xe4" => "ae")
#             line = replace(line, r"\xf6" => "oe")
#             line = replace(line, r"\xfc" => "ue")
#             line = replace(line, r"\xdf" => "ss")
#             push!(df, [line, "", ""])
#         elseif in_range && contains(line, "timeoffset")
#             in_range = false
#             line = replace(line, r"ß" => "ss")
#             line = replace(line, r"[,,]" => "")
#             line = replace(line, r"\xc4" => "Ae")
#             line = replace(line, r"\xd6" => "Oe")
#             line = replace(line, r"\xdc" => "Ue")
#             line = replace(line, r"\xe4" => "ae")
#             line = replace(line, r"\xf6" => "oe")
#             line = replace(line, r"\xfc" => "ue")
#             line = replace(line, r"\xdf" => "ss")
#             push!(df, [line, "", ""])
#         elseif in_range
#             line = replace(line, r"ß" => "ss")
#             line = replace(line, r"[,,]" => "")
#             line = replace(line, r"\xc4" => "Ae")
#             line = replace(line, r"\xd6" => "Oe")
#             line = replace(line, r"\xdc" => "Ue")
#             line = replace(line, r"\xe4" => "ae")
#             line = replace(line, r"\xf6" => "oe")
#             line = replace(line, r"\xfc" => "ue")
#             line = replace(line, r"\xdf" => "ss")
#             push!(df, [line, "", ""])
#         end
#     end

#     return df
# end

# newctl(raw"D:\Wasim\regio\control\rcm-c9_win_spin3.ctl")

z = ctlx()
if isempty(z)
    printstyled("no control file found!\n",color=:light_red)
    pretty_table(dfs,header=uppercasefirst.(names(dfs));)
    exit()
end

ofl = "routejl.txt"

if Sys.iswindows()
    routeg(z, ofl)
else
    zwsl = towsl(z)
    routeg(zwsl, ofl)
end
dx = dfroute(ofl=ofl)
dfs.Basin = parse.(Int64,dfs.Basin)
kd  = innerjoin(dfs, dx, on=:Basin)
# pretty_table(kd,header=uppercasefirst.(names(kd)); alignment=:l)#see below
# rm(ofl) #optional

dm = pwd()|>splitpath|>last 

#see: https://github.com/ronisbr/PrettyTables.jl/issues/127
open("routing-table-$dm.txt", "w") do f
    pretty_table(
        f,
        kd;header=uppercasefirst.(names(kd)), backend = Val(:text), alignment=:l
    )
end
@info "routing-table-$dm.txt written!"
#Base.write("routing-table-$dm.txt",sprint(io -> pretty_table(kd;header=uppercasefirst.(names(kd)), backend = Val(:text))))
writedf("routing-df-$dm.txt",kd)

Sys.iswindows() ? println("routingtable written to: ",pwd()) : println("routingtable written to: ",towsl(pwd()))
#Sys.iswindows() ? println("ja") : println("nein")

# tbl="routing-table-$dm.txt"
# Sys.iswindows() ? run(`type $tbl`) : run(`cat $tbl`)

# #output to screen

# lines = readlines("$tbl");
# for line in lines
#     println(line)
# end
pretty_table(kd,header=uppercasefirst.(names(kd)); alignment=:l)