
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

#outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"

file=ARGS[1];
#file="/mnt/d/Wasim/regio/out/v6/gcloud/vaporcm.mit.nc"
#file="/mnt/d/temp/saale/out_smf200/v2/windsmf200.2001.nc"
#lyr = isempty(ARGS) ? Int(1) : ARGS[2]
using Rasters,Plots
println("load",pwd()," ",file," ","an subset to Layer t=",lyr,"...")
#ts[t=lyr]|>plot
s=file
m=match(r".*[.]",s)
outfile = string(m.match,"png")

#x = read(Raster(dem_path))  
#ts=read(Raster(file,missingval=-9999))
ts=read(Raster(file,missingval=0))
x = ts[t=lyr]
#print("saving countour plot to ",outfile,"...")

printstyled("saving countour plot to $outfile...",color=:light_magenta,underline = true,blink = false,bold=true) 



#x = ts[t=1]
#x = x .> 0
#x=trim(x)

#mean_height = map(x -> mean(skipmissing(replace_missing(mask(dem_model; to=x)))), 	vect.geom)
#y = clamp(x.data,0)
# x
# x[Where(x)]
#w=classify(x)
#mean_val = mean(filter(!isnan, x))

# rmask = ts .> 0
# rmask = x .> 0 
# #g = mask(x;rmask)
# rmask = x .<= 0
# g = x .* rmask
# #g = x .- rmask
# #g = x .> 1
# g = x .> rmask
# g |> Plots.plot

# A = x
# #A[A.name(Where(x -> x > 15)), Y(Where(x -> x in (19, 21)))]
# A[A.name(Where(x -> x > 0))] |> Plots.plot
# A = ts

# A[A.data(Where(x -> x > 0))] |> Plots.plot

# w= replace_missing(x)
# w|>Plots.plot
#using Statistics
#d = zonal(mean, ts;of=:t)
#crs(x)
#mask(x,with=x>0)

# println(
# describe(x) #NOT in Rasters, but in DataFrames
# )

plotsize = (1600,800)
#p = contourf(x; dpi=300, size=(800, 400))
plt = contourf(x; 
	dpi=300, 
	size=plotsize,
	xlabel="",
	ylabel="",
	c=cgrad(:thermal))
# default(show = true)
# #p = plot(x,size=plotsize)
# display(p)
#savefig(p, "t.png")
savefig(plt, outfile)
#println(outfile," ... saved! \n"); nothing

printstyled("\n$outfile\t  ... saved! \n",color=:green,bold=true) ; nothing


# Julia also provides * for string concatenation:
# julia> greet * ", " * whom * ".\n"
# "Hello, world.\n"


# using Rasters
# # Create a mask to keep only values greater than -400
# mask = z .> -400
# # Apply the mask to the Raster
# g = z .* mask
# # Drop the spatial reference variable
# delete!(g, :spatial_ref)
# # Plot the raster at a specific time slice
# plot(g[20,:,:])




###########sysimage################################
# julia --startup-file=no

# using Pkg
# # Pkg.add("PackageCompiler")
# using PackageCompiler
# using Rasters,Plots
# pt,nm=("/home/ubu/.julia/JuliaImages","sys_raster_plots.so")
# pts = joinpath(pt,nm)
# create_sysimage([:Rasters,:Plots],sysimage_path=pts)

# ✔ [04m:20s] PackageCompiler: compiling incremental system image
# #sys_raster_plots.so             Tue Apr  4 15:17:29 2023        342.81MB  .      
#julia --startup-file=no --sysimage "/home/ubu/.julia/JuliaImages/sys_raster_plots.so" /mnt/c/Users/Public/Documents/Python_Scripts/julia/rplot.jl

###does not work for windows.##########################
#(base) PS D:\temp\saale\output\jan23\coarse-pest\coarse-pest\pestout> $f="C:/Users/Public/Documents/Python_Scripts/julia/tsplot.jl"
#(base) PS D:\temp\saale\output\jan23\coarse-pest\coarse-pest\pestout> julia --startup-file=no --sysimage $s $f .\windsmf180.2009.nc
# ERROR: could not load library "\\wsl.localhost\Ubuntu-20.04\home\ubu\.julia\JuliaImages\sys_js_ts_plots.so"
# %1 is not a valid Win32 application.

#	  

# using PackageCompiler
# using Rasters,Plots
# pt,nm=(pwd(),"winsys_raster_plots.so")
# pts = joinpath(pt,nm)
# create_sysimage([:Rasters,:Plots],sysimage_path=pts)


#The cairoplugin.dll file is not in the PATH environment variable.
# You can check if the GR bin directory is in your PATH by running ENV["PATH"] in the Julia REPL. 
# 	It should contain a path like C:\Users\YourName\.julia\packages\GR\SomeID\deps\gr\bin. 
# 	If not, you may need to add it manually by following these steps:
# 	 https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/

