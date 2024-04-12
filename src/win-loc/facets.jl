
if length(ARGS) == 0
	println("need args! <file>...")
    exit()
end

# if length(ARGS) == 1
# 	lyr=1
# 	println("skipping to layer 1...")
# else
# 	lyr=parse(Int,ARGS[2]);
# end

#outfile = contains(basename(file),".") ? string(m.match,"html") : basename(file)*".html"

file=ARGS[1];
#file=parse(String,ARGS[1])
println("searching for regex ",file," in ",pwd(),"...")

grids = filter(x -> occursin(Regex(file,"i"),x), readdir())
grids = filter(x->endswith(x,"nc"),grids)

if (length(grids)>0)
    file = first(grids)
    @warn("taking first match of $grids\n -> $file")
else
    dir=pwd()
    @warn "No netcdf match found in $dir !"
    exit()
end

outfile = string(file*".png")
#outfile = string(file*".svg")     #works!

#m=match(r".*[.]",s)
#outfile = string(m.match,"png")
printstyled("saving facet plot to $outfile...",
color=:light_magenta,underline = true,
blink = false,bold=false) 

using Rasters, Plots, Statistics

import NCDatasets

function facets(ext::AbstractString, outname::String)
    """
    like stackplot, but with grep
    """
           if isfile(file)
            r=read(Raster(file,missingval=0,mappedcrs=EPSG(25832)));
            fact=1.5        #.75
            if (r.dims[end]|>length == 1)
                @warn("only one layer available...")
#                println(DataFrames.describe(r))
                # a_min = minimum(skipmissing(r))
                # a_max = maximum(skipmissing(r))
                # a_var = var(skipmissing(r))
                # a_std = std(skipmissing(r))
                # ans=[string(name(r)),a_var,a_std,a_min,a_max]
                # #ans=hcat(ans...)  
                # ans=hcat(["name","var","std","min","max"],ans)
                # println(ans)
                try
                    a_min = minimum(skipmissing(r))
                    a_max = maximum(skipmissing(r))
                    a_var = var(skipmissing(r))
                    a_std = std(skipmissing(r))
                    ans=[string(name(r)),a_var,a_std,a_min,a_max]
                    ans=hcat(["name","var","std","min","max"],ans)
                    println(ans)
                catch
                @warn "all data of selection is NaN - exiting now!"
                return
                end
                p=Plots.plot(r;
                xlabel="",
                ylabel="",
                title=replace(basename(file),".nc"=>""),
                #c=cgrad(:thermal),
                c=cgrad(:matter),
                size=(1200*fact, 800*fact));
            else
                @warn("excluding first layer...")
                ee = Int(r.dims[3][end])
                rn = r[t=2:ee];    #subset till end
                p=Plots.plot(rn;
                xlabel="",
                ylabel="",
                #title=replace(basename(i),".nc"=>""), # title cause problems
                c=cgrad(:thermal),
                size=(1200*fact, 800*fact));
                #display(p)
            end
           savefig(p,outname)
           #println(basename(outname)," saved!");
           printstyled("$outname saved! \n",
           color=:green,underline = true,
           blink = false,bold=true) 
   end
end

facets(file,outfile)