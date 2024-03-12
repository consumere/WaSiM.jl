# 
function lf()
    directory = pwd()
    files = readdir(directory)  
    file_times = Dict{String, Dates.DateTime}()
    for file in files
        file_path = joinpath(directory, file)
        stat_info = stat(file_path)
        file_times[file] = Dates.unix2datetime(stat_info.mtime)
    end   
    file_times = collect(file_times)
    xf = sort(file_times, by = x -> x[2], rev=true)
    sz = []                 #::Vector{Float64}
    for i in 1:length(xf)
       push!(sz,round(stat(xf[i][1]).size/1024^2,digits=4))
    end
    merged_vector = []
    for (item1, item2) in zip(xf,sz)
        merged_item = (item1.first, item1.second,item2)
        push!(merged_vector, merged_item)
    end
    df = DataFrame(merged_vector);
    rename!(df,["latest_file","modified","size_in_mb"]);
    df.modified=map(x -> Dates.format(x, "yyyy-mm-dd HH:MM"),df.modified)
    return(df)
end

##very nice meta programming... 
macro vv(s) vgjl(s);end
#@vv "unter"
macro vpy(s) vgpy(s);end
#@vpy "climate"
macro vr(s) vgr(s);end
#@vr "climate"
macro vct(s) vgctl(s);end
#@vct "das ist"
macro rg(s) rglob(s);end
macro glb(s) glob(s);end
macro gl(s) glob(s)|>first;end
macro hd(df) df[1:4,:];end
#fastplot
macro fp(s) dfp(Regex(s));end
macro flog(s) dfl(Regex(s));end
macro ncrm() ncrem=src_path*"/ncremover.jl";include(ncrem);end
macro rasrm() remer=src_path*"/raster_remover.jl";include(remer);end
macro nco(s) first(nconly(s));end
macro dfo(s) first(dfonly(s));end
macro wajs() pt=src_path*"/wajs.jl";include(pt);end
macro jljs() pt=src_path*"/wajs.jl";include(pt);end
macro cmk() pt=src_path*"/cairomakie.jl";include(pt);end
macro rcall() pt=src_path*"/rcall.jl";include(pt);end
macro rgof() pt=src_path*"/RCall_gof.jl";include(pt);end
macro pj() pt=src_path*"/pyjl.jl";include(pt);end
macro pyjl() pt=src_path*"/pyjl.jl";include(pt);end

# using TOML
# #TOML is short for Tom’s Obvious Minimal Language and is a configuration file format 
# #that should be “easy to parse into data structures in a wide variety of languages”
# macro toml_str(s::String)
#     TOML.parse(s)::Dict{String, <:Any}
# end

# toml"""
# ki11 = 6.6
# dr11 = 1
# kb11 = 0.4
# q011 = 0.02
# """

# dx = toml"""ki11 = 6.6 
# dr11 = 1
# kb11 = 0.4
# q011 = 0.02"""|>DataFrame
# #DataFrame(dx)
# findall
# findfirst
# findlast
# findnext
# findprev

# fx=raw"D:\Wasim\regio\control\rcm200_r6.ctl"
# for (i, line) in enumerate(eachline(fx))
#     if findfirst("dr11", line) !== nothing
#         println("Line $i: $line")
#     end
# end

