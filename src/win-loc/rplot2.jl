
using DaemonMode

if length(ARGS) == 0
	println("need args! <file><opt:layer>...")
    exit()
end

if length(ARGS) == 1
	lyr=1
	println("skipping to layer 1...")
else
	lyr=parse(Int,ARGS[2]);
end

#runexpr("using CSV, DataFrames")
#fname = "tsp_50.csv";
# runexpr("""begin
      # df = CSV.File("$fname") |> DataFrame
      # println(last(df, 3))
  # end""")

file=ARGS[1];

runexpr("using Rasters,Plots")
runexpr("""begin
println("load",pwd()," ",$file," ","an subset to Layer t=",$lyr,"...")
s=file
m=match(r".*[.]",s)
outfile = string(m.match,"png")
ts=read(Raster($file,missingval=-9999))
x = ts[t=$lyr]
print("saving countour plot to",outfile,"...")
plotsize = (1600,800)
p = contourf(x; dpi=300,size=plotsize,c=cgrad(:thermal)) 
savefig(p, outfile)
print(outfile,"... saved! \n"); nothing
  end""")
