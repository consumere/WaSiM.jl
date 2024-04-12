#plotly time Series
if length(ARGS) <= 0
	println("need args! <timeseries file>...")
    exit()
end

#julia --threads auto -q C:\Users\Public\Documents\Python_Scripts\julia\tsplot_win.jl $x


file=ARGS[1];

if endswith(file,".nc")
	println("need <timeseries file>...")
    exit()
end

#path="D:/Wasim/regio/out/lowres/c4_rev/tso"
pt = split(file)[1]
m = match(r".*[.]",basename(file))
#outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"
#outfile = contains(basename(file),".") ? string(m.match,"png") : basename(file)*".png"
outfile = contains(basename(file),".") ? string(m.match,"svg") : basename(file)*".svg"
println("loading...\n",pt,"\nand save it to ",outfile,"...\n")

using StatsPlots, DataFrames, CSV, Dates

using Printf

function loaddf(path::AbstractString)
    ms="-9999"
    df = CSV.read(path,DataFrame,
    missingstring=ms,
    delim="\t",comment="-",
    silencewarnings=false,                                         
    normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    df=df[:,Not(1:3)]
end

df = loaddf(file)
println(describe(df))
# s = propertynames(df)[1:end-1]
# # p1 = @df df plot(:date, cols(s),size=plotsize)
# # p2 = @df df density(cols(s), legend = :topleft,size=plotsize)
# @df df plot(:date, cols(s),size=plotsize)
# @df df plot(:date, cols(s), yaxis=:log,size=plotsize)
# @df df corrplot(cols(s),size=plotsize)  
# ##@df df bar!(cols(s),size=plotsize,legend=false)
# x,y = df[!,1],df[!,2]
# density(x,legend = :topleft,size=plotsize)
# density!(y,legend = :topleft,size=plotsize)

# plot(qqplot(df,cols(s)), corrplot(df, cols(1:2)))

# using Distributions 
# # compare with a Cauchy distribution fitted to y; pass an instance (e.g. Normal(0,1)) to compare with a specific distribution
# # the :R default line passes through the 1st and 3rd quartiles of the distribution
# plot( 
# qqplot(df[!,1],df[!,2], qqline = :fit), 
# qqplot(Cauchy,df[!,2]), 
# qqnorm(df[!,2], qqline = :R) 
# )

# @df df qqplot(cols(s), legend = :topleft,size=plotsize)

#for qouts
#df=loaddf(fn)
#@df df corrplot(cols(1:2), grid = false)      
#out = pline(file)

#fact=0.88
fact=2
plotsize = (1600*fact,800*fact)

function vio(df::DataFrame)
	str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
	ln = Symbol.(filter(x->!occursin("date",x),names(df)))
	@df df StatsPlots.violin(str,cols(ln),linewidth=0.1,plotsize=plotsize)
end



println("saving violin plot to",outfile,"...")
savefig(vio(df),outfile)
println("done! ...")


# using PlotlyJS, DataFrames, CSV, Dates
	
	# function pline(path::AbstractString)
	    # ms="-9999"
	    # df = CSV.read(path,DataFrame,
	    # missingstring=ms,
	    # delim="\t",comment="-",
	    # silencewarnings=false,                                         
	    # normalizenames=true,drop=(i, nm) -> i == 4) |> dropmissing
	    # df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
	        # df=df[:,Not(1:3)]
	    # nrows=size(df)[2]-1
	    # st=[]
	    # for i in 1:size(df)[2]-1; push!(st,string(propertynames(df)[i]));end
	    # p = make_subplots(rows=nrows, cols=1, 
	    # shared_xaxes=true, 
	    # shared_yaxes=false,
	    # vertical_spacing=0.05,
	    # )
	    # for i in 1:nrows;
	            # add_trace!(p, 
	                        # scatter(x=df.date, y=df[:,i],
	            # name=st[i]),   row=i,     col=1);
	    # end
	# end
# out = pline(file)
# println("saving plotly plot to",outfile,"...")
# savefig(out,outfile)
# println("done! ...")

#######make image ##############
# using StatsPlots, DataFrames, CSV, Dates

# pt,nm=(pwd(),"winsys_ts_plots.so")
# pts = joinpath(pt,nm)
# create_sysimage([:StatsPlots, :DataFrames, :CSV, :Dates],sysimage_path=pts)

