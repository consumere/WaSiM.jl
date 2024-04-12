#plotly time Series
if length(ARGS) <= 0
	#println("need args! <timeseries file regex!>...")
	println("need args! <timeseries file!>...")
    exit()
end

file=ARGS[1];
#file="/mnt/d/Wasim/Goldbach/revision/fab30/v1/sb05fab.v0.2021"


if endswith(file,".nc")
	println("need <timeseries file regex>...")
    exit()
end

pt = split(file)[1]
m = match(r".*[.]",basename(file))
#outfile = contains(basename(file),".") ? string(m.match,"png") : basename(file)*".png"
outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"
#println("loading...\n",pt,"\nand save it to ",outfile,"...\n")

printstyled(
    "loading...\n",pt,"\nand save it to ",outfile,"...\n",
    color=:green,
    bold=true)


#using Query, DataFrames, CSV, Dates, PlotlyJS
using PlotlyJS, DataFrames, CSV, Dates
#set_theme("ggplot2")
#PlotlyJS.set_theme("seaborn")
#theme="plotly_dark")



function readdf(x::AbstractString)
    ms=["-9999","lin","log"];
    df = CSV.read(x,
    DataFrame,
    missingstring=ms,
    ntasks=4,
    limit=typemax(Int),
    types = Float64,
    delim="\t",
    silencewarnings=true,
    normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
    df.YY=map(x ->Int(x),df.YY);
    df.MM=map(x ->Int(x),df.MM);
    df.DD=map(x ->Int(x),df.DD);
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
    metadata!(df, "filename", x, style=:note);
end

# function readdf(x::Regex)
#     """
#     readdf(x::Regex)
#     reads first match of regex wasim timeseries
#     """
#     rootdir=pwd()
#     results = []
#     for (looproot, dirs, filenames) in walkdir(rootdir)
#         for filename in filenames
#             if (occursin(x,filename)) && (!occursin(r"txt|yrly|nc|png|svg|grd",filename))
#                 push!(results, joinpath(looproot, filename)) 
#             end
#         end
#     end
#     x = first(results)
#     ms=["-9999","lin","log"]
#     df = CSV.read(x,DataFrame,missingstring=ms,
#     ntasks=4,
#     limit=typemax(Int),
#     types = Float64,
#     delim="\t",
#     silencewarnings=true,
#     normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
#     df.YY=map(x ->Int(x),df.YY);
#     df.MM=map(x ->Int(x),df.MM);
#     df.DD=map(x ->Int(x),df.DD);
#     df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
#     df=df[:,Not(1:3)]
#     metadata!(df, "filename", x, style=:note);
# end

function dfpjs(df::DataFrame;)
    nrows=size(df)[2]-1 
    ti = DataFrames.metadata(df)|>only|>last|>basename
    fig = PlotlyJS.make_subplots(
        shared_xaxes=true, 
        shared_yaxes=true    
        );
    for i in 1:nrows;
        PlotlyJS.add_trace!(fig, 
        PlotlyJS.scatter(
            #theme="seaborn",
            theme="ggplot2",
            x=df.date, 
            y=df[:,i],
            name=names(df)[i]));
    end
    fact,logy = 1.2,1
    if logy == true
        PlotlyJS.relayout!(fig,yaxis_type="log",
        height=600*fact,width=900*fact,
        title_text="Series of "*ti)
    else
        PlotlyJS.relayout!(fig,
        height=600*fact,width=900*fact,
        title_text="Series of "*ti)
    end
    fig
    #display(fig)
end

#out = dfpjs(readdf(Regex(file)))

inpath=abspath(file)
df=readdf(inpath)
println(describe(df))

#println("saving plotly plot to",outfile,"...")

printstyled("saving plotly plot to $outfile...",color=:green,bold=true)

savefig(dfpjs(df),outfile)

printstyled("\n ...done!\n",color=:green,underline=true)
