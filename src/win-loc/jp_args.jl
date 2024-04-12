
path=ARGS[1];

function asnullable(x)
   isa(x, Nullable) ? x : Nullable(x)
end
# asnullable(path)

#if path == nothing
#println("error!")
#end

println("file: ",path);
# println("file start stop skip sep");
println("file start stop skip");


#if isnan(skip)
start=parse(Int, ARGS[2])
stop=parse(Int, ARGS[3])

#skip=nothing
skip=parse(Int, ARGS[4])
if skip == nothing
	skip=5
	println("skipping 5")
end

#sep=parse(String, ARGS[5])
#sep=nothing
# sep=ARGS[5]
# if sep == nothing
	# sep="\t"
	# println("sep is tab!")
# end

println("..loading packages")
using DataFrames, CSV, Dates  
println("..loading dataframe")


# df = DataFrame(CSV.File(path, header=skip, delim=sep,missingstring="-9999"));
# df = DataFrame(CSV.File(path, header=skip,missingstring="-9999"));

df = CSV.read(path, DataFrame, header=skip,missingstring="-9999")

df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD,"-",df.HH),"yyyy-mm-dd-HH");

println(path,"\neval on:\n",(propertynames(df)[5:end-1]),"\n on years",start," to ",stop)
df=df[in(start:stop).(year.(df.date)),Not(1:4)]

println("..loading StatsPlots")
using StatsPlots
# gr()
plotsize = (1600,800)
@df df plot(:date, propertynames[1:end-1],size=plotsize)

#plot(df[in(2007:2008).(year.(df.date)), :date],df[in(2007:2008).(year.(df.date)), :4], label = names(df)[4], color = "blue")

fn=string(path,"_jplot",".png")
savefig(fn)
println("saved as: ",fn)