cd("E:/qq/jlcor/")
using Printf
using Statistics

function process_files(;suffix=".winlist")
    #files = filter(x -> contains(x, "cor") && endswith(x, ".nc"), readdir())
    #files = filter(x -> occursin(r"p+.*nc",x) && endswith(x, ".nc"), readdir())
    files = filter(x -> occursin(r"daily+.*nc",x) && endswith(x, ".nc"), readdir())

    if any(x -> occursin("cor", x) && endswith(x, suffix), readdir())
        println("Corresponding files already exist. Exiting now!")
        return
    end
    

    for var in files
        gridfile = string(replace(var,".nc" => ".gridfile"))
        global onam = string(replace(var,".nc" => suffix))

        #cmds = `cdo griddes`
        
        #run(`cdo griddes $var > $gridfile`)
        #run(pipeline(cmds))
        #pipeline(`cdo griddes `,$var, `>`, $gridfile)
        open(gridfile, "w") do out
            #write(out, String(read(pipeline(`cdo griddes $var`, `cat`))))
            write(out, String(read(pipeline(`cdo griddes $var`))))
        end
        

        open(onam, "w") do file
            println(file, "GRID LIST of $var")
            println(file, "YY\tMM\tDD\tHH\t9999")
            print(file, "YY\tMM\tDD\tHH\t")
            println(file, read_grid_value(gridfile, "xfirst"))
            #print("\n")
            print(file, "YY\tMM\tDD\tHH\t")
            println(file, read_grid_value(gridfile, "yfirst"))
           # print("\n")
            println(file, "YY\tMM\tDD\tHH\tListe")

            ts2 = "ts2"  # Replace this with the actual variable or data source
            lines = readlines(ts2)
            for (i, line) in enumerate(lines)
                global param = split(var, '-')[1]

                #awk '{OFS=FS}{$NF=FS"D:/remo/qm/corgrids/jlcor/outw.nc<param>"i-1+1;i++;print}' ts2 >> $onam
                #line = replace(line, "outw.nc" => var)
                #line = replace(line, "param" => param)
                #/jlcor/pre-cor.nc<pre>0
                line = replace(line, "00" => string("\tD:/remo/qm/corgrids/jlcor/",var,"<",param,">",i-1))

                println(file, replace(line, r"\s+" => "\t"))
            end
        end

        println("$param written to $onam !")
        println("$onam done!")
    end

    println("Lists for WaSiM Input are written to the current directory. All done!")
end

function read_grid_value(gridfile, pattern)
    lines = readlines(gridfile)
    for line in lines
        #if contains(pattern, line)
        if contains(line,pattern)
            #match_obj = occursin(r"$pattern\s*:\s*([^\s]+)", line)
            #return isnothing(match_obj) ? "N/A" : match_obj.captures[1]
            match_obj = split(line,"=")|>last
            return isnothing(match_obj) ? "N/A" : strip(match_obj)
        end
    end
    return "N/A"
end

######################xr method

#gridlist for WaSiM
# Reihenfolge ist sehr wichtig
# 1. xfirst
# 2. yfirst
# 3. gridres

import xarray as xr

function xrlist(;cwd::AbstractString=pwd(),gridres::Float64=8317.01,xmatch="qq.",suffix=".winlist")
    files = filter(x -> contains(x, xmatch) && endswith(x, ".nc"), readdir())

    if any(x -> occursin(xmatch, x) && endswith(x, suffix), readdir())
        println("Corresponding files already exist. Exiting now!")
        return
    end

    for var in files
        #gridfile = string(replace(var,".nc" => ".gridfile"))
        global onam = string(replace(var,".nc" => suffix))
        # Load the NetCDF file into an xarray Dataset
        ds = xr.open_dataset(var)
        println("Processing $var ...")

        # # Access the coordinate information
        # lon_values = ds["lon"].values
        # lat_values = ds["lat"].values

        xf = ds["lon"].values|>first
        #xl = ds["lon"].values|>last
        yf = ds["lat"].values|>first
        #yl = ds["lat"].values|>last

        # gridfile_content = ["xfirst:  $xf",
        #                     "xlast:   $xl",
        #                     "yfirst:  $yf",
        #                     "ylast:   $yl"]

        gridfile_content = ["GRID\tLIST\tof\t$var", #var is var in files
                            "YY\tMM\tDD\tHH\t$xf",
                            "YY\tMM\tDD\tHH\t$yf",
                            "YY\tMM\tDD\tHH\t$gridres",
                            "YY\tMM\tDD\tHH\tListe\n"]

        # open("x_gridfile.txt", "w") do out
        #     write(out, join(gridfile_content, "\n"))
        # end    
        
        open(onam, "w") do out
            write(out, join(gridfile_content, "\n"))
        end

        timevec = ds["time"].values
        #typeof(timevec|>first)
        dts = nothing
        try
            dts = [Dates.DateTime.(t.year, t.month, t.day, t.hour, t.minute) for t in timevec]
        catch
            @warn "timevec is not ctftime"
            dts = timevec.astype("M8[m]").tolist()
        end      
        # #string(s.year, s.month, s.day, s.hour, s.minute)
        # #stateini/jlcor/pre-cor.nc<pre>0
        #dk = DataFrame(time=dts,link="stateini/jlcor/pre-cor.nc",param="<pre>",cnt=(1:length(dts)).-1)
        param = ds.keys()|>collect|>last
        dk = DataFrame(YY=year.(dts),MM=month.(dts),DD=day.(dts),HH=hour.(dts),
            #link="stateini/jlcor/$var",
            #link="D:/remo/qm/corgrids/jlcor/$var",
            link = joinpath(cwd,var),
            param = "<$param>",
            cnt = (1:length(dts)).-1)
        
        
        dout = hcat(dk[:, 1:4], [join(row, "") for row in eachrow(dk[:, 5:end])])
        #onam
        #append to onam file
        CSV.write(onam,dout, 
            transform = (col, val) -> something(val, missing),
            delim="\t",append=true)  

            println("$param written to $onam !")
            println("$onam done!")
    end

    println("Lists for WaSiM Input are written to the current directory. All done!")
end

function xrlist2(;cwd::AbstractString=pwd(),xcrd="x",ycrd="y",tvec="t",gridres::Float64=8317.01, xmatch="qq2",suffix=".winlist")
    files = filter(x -> contains(x, xmatch) && endswith(x, ".nc"), readdir())
    if any(x -> occursin(xmatch, x) && endswith(x, suffix), readdir())
        println("Corresponding files already exist. Exiting now!")
        return
    end

    for var in files
        global onam = string(replace(var,".nc" => suffix))
        ds = xr.open_dataset(var)
        println("Processing $var ...")
        # # Access the coordinate information
        #xf = ds["lon"].values|>first
        xf = ds[xcrd].values|>first
        #yf = ds["lat"].values|>first
        #xf = round(xf, digits=3)
        yf = ds[ycrd].values|>first
        #yf = round(yf, digits=3)
        gridfile_content = ["GRID\tLIST\tof\t$var", #var is var in files
        "YY\tMM\tDD\tHH\t$xf",
        "YY\tMM\tDD\tHH\t$yf",
        "YY\tMM\tDD\tHH\t$gridres",
                            "YY\tMM\tDD\tHH\tListe\n"]
        open(onam, "w") do out
            write(out, join(gridfile_content, "\n"))
        end
        #timevec = ds["time"].values
        timevec = ds[tvec].values
        dts = nothing
        try
            dts = [Dates.DateTime.(t.year, t.month, t.day, t.hour, t.minute) for t in timevec]
        catch
            @warn "timevec is not ctftime"
            dts = timevec.astype("M8[m]").tolist()
        end      
        param = ds.keys()|>collect|>last
        dk = DataFrame(YY=year.(dts),MM=month.(dts),DD=day.(dts),HH=hour.(dts),
            #link="stateini/jlcor/$var",
            #link="D:/remo/qm/corgrids/jlcor/$var",
            link = joinpath(cwd,var),
            param = "<$param>",
            cnt = (1:length(dts)).-1)
        dout = hcat(dk[:, 1:4], [join(row, "") for row in eachrow(dk[:, 5:end])])
        #append to onam file
        CSV.write(onam,dout, 
            transform = (col, val) -> something(val, missing),
            delim="\t",append=true)  
            println("$param written to $onam !")
            println("$onam done!")
    end
    println("Lists for WaSiM Input are written to the current directory. All done!")
end

function grid_resolution(ds)
    x_resolution = ds["x"][2] - ds["x"][1]
    y_resolution = ds["y"][2] - ds["y"][1]
    return (x_resolution, y_resolution)
end

ds = xr.open_dataset("tas-qq.nc")
#check grid resolution
xres = grid_resolution(ds)
first(xres)
last(dts)
xrlist()
npplat()
pwd()|>cb

#now test
#D:\Wasim\regio\control\rcm200_y4-loc.ctl

#to unixlist
pwd()
v = glob("2.winlist")
for x in v
    y = replace(x,".winlist" => ".unixlist")
    cp(x,"../unixlists/$y")
end
cd("../unixlists")
# replace path E:\qq\jlcor\ in each unixlist to stateini/
# like this but in julia
# perl -i.bak -pe 's#E:\qq\jlcor\#stateini/#g' *list

function replace_path_in_file(filename::String, old_path::String, new_path::String)
    # Read the file
    content = read(filename, String)

    # Replace the old path with the new path
    new_content = replace(content, old_path => new_path)

    # Write the new content back to the file
    open(filename, "w") do f
        write(f, new_content)
    end
    println("Replaced paths in $filename")
end

# Use the function
map(fl->replace_path_in_file(fl, "E:\\qq\\jlcor\\", "stateini/"), readdir())
npplat()

#model output is all NaN #maby 10-07-02 nochmal
cd("D:/Wasim/regio/out/rc200/y4/loc/")
dfr(r"sb")
dfr(r"wu")

xrlist(xmatch="qq2",suffix=".winlist")

##
pwd()
cd("E:\\qq\\rainfarm")
nconly("p")
xr = pyimport("xarray")
xrlist(xmatch="pre_biljl",suffix=".winlist")
l=lat()
npp(l)
cp(l,replace(l,".winlist" => ".unixlist"))
N=replace(l,".winlist" => ".unixlist")
replace_path_in_file(N, "E:\\qq\\rainfarm\\", "stateini/")
npplat()

da = xr.open_dataset("pre_biljl.nc")
da.keys()|>collect
#da = da.rename(Dict("lon_2"=>"x","lat_2"=>"y"))
da = da.rename(Dict("lon"=>"x","lat"=>"y"))
grid_resolution(da)
