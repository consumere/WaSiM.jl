#mackie time Series
if length(ARGS) <= 0
	println("need args! <timeseries file>...")
    exit()
end

#julia --threads auto -q C:\Users\Public\Documents\Python_Scripts\julia\tsmackie.jl $x

file=ARGS[1];

if endswith(file,".nc")
	println("need <timeseries file>...")
    exit()
end
pt = split(file)[1]
m = match(r".*[.]",basename(file))

#outfile = contains(basename(file),".") ? string(m.match,"png") : basename(file)*".png"
outfile = contains(basename(file),".") ? string(m.match,"svg") : basename(file)*".svg"

println("loading...\n",pt,"\nand save it to ",outfile,"...\n")
using CairoMakie #, Makie    
println("using CairoMakie ...")
# import Makie:wong_colors      #<-only for tsp, tsp2
# using DataFrames, CSV, Statistics, Dates, Distributions
#using DataFrames, Dates #, Grep , Printf
using Dates
import DelimitedFiles:readdlm
import DataFrames:DataFrame
import DataFrames:metadata
import DataFrames:metadata!
import DataFrames:Not
import DataFrames:select

function wread(x::String;skip=3)
    """
    Read wasim ts with DelimitedFiles.readdlm, skipto line 3 
    no header column
    """
    #df = DelimitedFiles.readdlm(x, '\t', Float64, '\n';
    df =readdlm(x, '\t', Float64, '\n';
        header=false,skipstart=skip)
    df = DataFrame(df,:auto)
    for i in 5:size(df,2)
        df[!,i]=replace(df[!,i],-9999.0 => missing)
    end 
    for i in 5:size(df,2)
        replace!(df[!,i],-9999.0 => missing)
    end
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df=df[:,Not(1:4)]
    metadata!(df, "filename", x, style=:note);
end

# function loaddf(path::AbstractString)
#     ms="-9999"
#     df = CSV.read(path,DataFrame,
#     missingstring=ms,
#     delim="\t",comment="-",
#     silencewarnings=false,                                         
#     normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
#     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#     df=df[:,Not(1:3)]
# end

#df = loaddf(file)
df = wread(file)

function tsp(df::DataFrame)
    """
    plots timeseries
    """
    dt = df.date
    fig = Figure(resolution=(600, 400), 
            fonts=(;regular = "consolas")            
            )

    
    tempo = string.(dt)
    lentime = size(df, 1)
    slice_dates = range(1, lentime, step=lentime ÷ 8)
    tit = DataFrames.metadata(df) |> only |> last |> basename
    
    ax3 = Axis(fig[1, 1], 
               title = replace(tit, r"_|so" => " "),
               xlabel = "Date",
               ylabel = first(names(df)))
    
    cols = names(df)
    filter!(x -> !occursin(r"date|year", x), cols)
    colors = wong_colors(length(cols))
    
    for (i, col) in enumerate(cols)
        vals = select(df, col) |> Matrix |> vec
        lines!(ax3, 1:lentime, vals; color = colors[i], linewidth = 0.85)
    end
    
    ax3.xticks = (slice_dates, tempo[slice_dates])
    ax3.xticklabelrotation = π / 4
    ax3.xticklabelalign = (:right, :center)
    
    return fig
end

function tsp2(df::DataFrame)
    """
    plots timeseries for each column separately in one figure
    """
    dt = df.date
    tempo = string.(dt)
    lentime = size(df, 1)
    slice_dates = range(1, lentime, step=lentime ÷ 8)
    tit = DataFrames.metadata(df) |> only |> last |> basename
    cols = filter(x -> !occursin(r"date|year", x), names(df))

    fig = Figure(resolution=(800, 400), 
    fonts=(;regular = "consolas")            
    )
                 
    
    for (i, col) in enumerate(cols)
        ax = Axis(fig[1, i], 
                   title = replace(tit, r"_|so" => " "),
                   xlabel = "Date",
                   ylabel = col)
        
        vals = select(df, col) |> Matrix |> vec
        lines!(ax, 1:lentime, vals; linewidth = 0.85)
        
        ax.xticks = (slice_dates, tempo[slice_dates])
        ax.xticklabelrotation = π / 4
        ax.xticklabelalign = (:right, :center)
    end
    
    return fig
end

function tsp3(df::DataFrame)
    """
    plots timeseries for each column separately in one figure
    """
    dt = df.date
    tempo = string.(dt)
    lentime = size(df, 1)
    slice_dates = range(1, lentime, step=lentime ÷ 8)
    tit = metadata(df) |> only |> last |> basename
    cols = filter(x -> !occursin(r"date|year", x), names(df))

    fig = Figure(resolution=(1024, 800).*1.1, 
    fonts=(;regular = "consolas")            
    )   
                 
    for (i, col) in enumerate(cols)
        println("add facet $col ...")
        row, cl = fldmod1(i, 2)
        ax = Axis( fig[row, cl],
                   title = replace(tit, r"_|so" => " "),
                   xlabel = "Date",
                   ylabel = col)
        
        vals = select(df, col) |> Matrix |> vec
        lines!(ax, 1:lentime, vals; linewidth = 0.85)
        
        ax.xticks = (slice_dates, tempo[slice_dates])
        ax.xticklabelrotation = π / 4
        ax.xticklabelalign = (:right, :center)
    end
    
    return fig
end

plt = tsp3(df)

println("saving plot to ",outfile,"...")
#savefig(plt, outfile)
#save(filename, data...) saves the contents of a formatted file, trying to infer the format from filename.
CairoMakie.save(outfile,plt)
println("done! ...")

# #######make image ##############
# cd("D:/Wasim/Docker/JuliaImagesWin")
# #activate .  
# add PackageCompiler
# add CairoMakie
# add Dates
# add DataFrames,DelimitedFiles

# using PackageCompiler
# using CairoMakie
# using Dates
# using DataFrames,DelimitedFiles

# pt,nm=(pwd(),"winsys_cmakie.so")
# pts = joinpath(pt,nm)
# create_sysimage([:CairoMakie, 
# :DataFrames,:DelimitedFiles,
# :Dates],sysimage_path=pts)

# :readdlm,
# :DataFrame,
# :metadata,
# :metadata!,
# :Not,
# :select,

# pt,nm=(pwd(),"winsys_rasters.so")
# pts = joinpath(pt,nm)
# create_sysimage(sysimage_path=pts)