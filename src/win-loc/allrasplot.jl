if length(ARGS) == 0
	println("need args! <fileregex>...")
    exit()
end

# if length(ARGS) == 1
# 	lyr=1
# 	println("skipping to layer 1...")
# else
# 	lyr=parse(Int,ARGS[2]);
# end

function nconly(x1::AbstractString)
    v = readdir();
    v = v[broadcast(x->endswith(x,"nc"),v)];
    z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
    return(z)
    end
    
fn = ARGS[1];
#"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m2"|>cd
#fn="te"
files=nconly(fn)

if length(files) == 0
	println("no nc files in current dir present...")
    exit()
end

println("This plots all $files \n loading packages","...")


using Rasters,Plots


#ts[t=lyr]|>plot
lyr=1;
plotsize = (1600,800);
for file in files
    s=file
    m=match(r".*[.]",s)
    outfile = string(m.match,"png")
    println("load ",file," with missingval=0 and subset to Layer t=",lyr,"...")
    #ts=read(Raster(file,missingval=-9999))
    ts=read(Raster(file,missingval=0))
    x = ts[t=lyr]
    println("saving countour plot to ",outfile,"...")
    p = contourf(x; 
    	dpi=300, 
    	size=plotsize,
    	xlabel="",
    	ylabel="",
        title=replace(basename(file),".nc"=>""),
    	c=cgrad(:thermal;scale=:log,rev=true))
    #Plots.plot(xr;c=cgrad(:thermal;scale=:log,rev=true),
    savefig(p, outfile)
    println(outfile,"... saved! \n"); nothing
end

print("all done!")

#"D:/Wasim/Docker/JuliaImagesWin"
# using PackageCompiler
# using Rasters,Plots
# pt,nm=(pwd(),"sys_raster_plots.so")
# pts = joinpath(pt,nm)
# create_sysimage([:Rasters,:Plots],sysimage_path=pts)
##✔ [04m:24s] PackageCompiler: compiling incremental system image