cd("D:/remo/genRE/")
using Printf
using Statistics
@pyjl
using PyCall
xr = pyimport("xarray")
rio = pyimport("rioxarray")
r = xr.open_dataset(raw"D:/remo/genRE/pre-gen.nc")
r = xr.open_dataset(raw"D:\remo\genRE\genRE_main_25832.nc")
r["pr"].isel(time=3).plot()
function grid_resolution(ds)
    x_resolution = ds["x"][2] - ds["x"][1]
    y_resolution = ds["y"][2] - ds["y"][1]
    return (x_resolution, y_resolution)
end

grid_resolution(r)

#now make a list.

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
        gridfile_content = [
            "GRID\tLIST\tof\t$var", #var is var in files
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


collect(r.keys())
collect(r.variables)
xrlist2(xmatch="25832",suffix=".winlist",tvec="time",gridres=1200.0)
v=lat()
npp(v)
v
cp(v,"genRE_main_24h.winlist")
pwc()
#wsl
#p024 genRE_main_24h.winlist



# #to unixlist
# pwd()
# v = glob("2.winlist")
# for x in v
#     y = replace(x,".winlist" => ".unixlist")
#     cp(x,"../unixlists/$y")
# end
# cd("../unixlists")



# perl -i.bak -pe 's#E:\qq\jlcor\#stateini/#g' *list


