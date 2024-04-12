using Printf
using Statistics
#cd("/mnt/d/remo/genRE")
#a="genRE_daily.nc"
a="genRE_precipitation_hour_1.1.nc"
# run(`cdo griddes $var > $gridfile`)
# run(`cdo showtimestamp $a|tr -s ' ' '\n'|head -1`)
a=readdir()[2]
#run(pipeline(`cdo showtimestamp $a`,stdout="tst"))
#run(pipeline("tst", `tr -s ' ' '\n'`))

#######das geht gut
# result = read(pipeline("tst", `tr -s ' ' '\n'`), String)
# output_array = split(result[2:end-1], '\n')
# using Dates
# oa = parse.(DateTime, output_array)
# oa = Dates.format.(oa, "yyyy mm dd HH MM")
# oa = split.(oa, ' ')

#example for pipeline(pipeline(a,b),c)
#cmd= pipeline(`cdo showtimestamp $a`, stdout="xst"),pipeline(`tr -s ' ' '\n'`, stdin="xst", stdout="tst")

###make ts2
# # there is a cdo command already inside this loop 
open("ts2", "w") do out
    #write(out, String(read(pipeline(`cdo griddes $var`, `cat`))))
    ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n")
    #ts_content = replace(ts_content, r"-|T|:" => "\t") #this is for daily... 
    ts_content = replace(ts_content, r"-|:" => "\t")
    ts_content = replace(ts_content, "T" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "\n" => "", count=1)
    write(out, ts_content)
end

# function make_ts(a::String)
#     ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
#     ts_content = replace.( ts_content, r"\s+" => "\n")
#     ts_content = replace.( ts_content, r"-|:" => "\t")
#     ts_content = replace.( ts_content, r"T" => "\t")
#     ts_content = replace.( ts_content, r"^\s+\S+$" => "")
#     ts_content = replace.( ts_content, r"^\s+.*$" => "")
#     x = [join(split(line, "\t")[1:5], "\t") for line in ts_content]
#     open("ts2", "w") do out
#         aout = join(x, "\n")
#         aout = replace.(aout, "\n" => "", count=1)
#         write(out, aout)
#     end
# end
##make_ts(a)
run(`head -5 ts2`)
run(pipeline(`cut -f1-5 ts2`, stdout="ts.txt"))
run(`head -5 ts.txt`)

"""
what::Regex=r"daily";suffix=".winlist",tsfile="ts.txt"
"""
function process_files(what::Regex=r"daily";suffix=".winlist",tsfile="ts.txt")
    #files = filter(x -> contains(x, "cor") && endswith(x, ".nc"), readdir())
    #files = filter(x -> occursin(r"p+.*nc",x) && endswith(x, ".nc"), readdir())
    files = filter(x -> occursin(what,x) && endswith(x, ".nc"), readdir())

    if any(x -> occursin(what, x) && endswith(x, suffix), readdir())
        println("Corresponding files already exist. Exiting now!")
        return
    end

    for var in files
        gridfile = string(replace(var,".nc" => ".gridfile"))
        global onam = string(replace(var,".nc" => suffix))
        open(gridfile, "w") do out
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
            lines = readlines(tsfile)                  ##
            for (i, line) in enumerate(lines)
                global param = split(var, '-')[1]
                line = replace(line, "\t00" => string("\tD:/remo/genRE/",var,"<",param,">",i-1))
                println(file, replace(line, r"\s+" => "\t"))
            end
        end

        println("$param written to $onam !")
        println("$onam done!")
    end

    println("Lists for WaSiM Input are written to the current directory. All done!")
end


k  = glob("gri")|>first
run(`cat $k`)
read_grid_value(k, "yfirst")
read_grid_value(k, "xfirst")

vardes pre_daily.nc 
griddes pre_daily.nc 
variables from ./pre_daily.nc:
['time', 'longitude', 'latitude', 'pr']



function read_grid_value(gridfile, pattern)
    lines = readlines(gridfile)
    for line in lines
        if contains(line,pattern)
            match_obj = split(line,"=")|>last
            return isnothing(match_obj) ? "N/A" : strip(match_obj)
        end
    end
    return "N/A"
end

process_files(r"hour")


#now a tst for the daily nc
a=nconly("32632")|>only
open("dts", "w") do out
    ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n")
    ts_content = replace(ts_content, r"-|T|:" => "\t") #this is for daily... 
    #ts_content = replace(ts_content, r"-|:" => "\t")
    ts_content = replace(ts_content, "T" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "\n" => "", count=1)
    write(out, ts_content)
end
@doc process_files(r"32632";tsfile="dts")
process_files(r"32632";tsfile="dts")
x=lat()
readlines(x)[1:15]

@edit pwc()
#run(pipeline(`wslpath -m $x`, stdout=clipboard))

Main.pyjl.npp("dts")
k=lat()
Main.pyjl.npp(k)


genRE_daily.nc


#now a tst for the pre_daily.nc
a=nconly("pre-daily")|>only
open("ts2", "w") do out
    ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n")
    ts_content = replace(ts_content, r"-|T|:" => "\t") #this is for daily... 
    #ts_content = replace(ts_content, r"-|:" => "\t")
    ts_content = replace(ts_content, "T" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "00\t00" => "00")       #like cut -f1-5 ts2
    ts_content = replace(ts_content, "\n" => "", count=1)
    write(out, ts_content)
end

process_files(r"pre-daily";tsfile="ts2")
x=lat()
pyjl.npp(x)
rm(x)

var=r"pre_daily"|>glob|>first
split(var, '-')[1]

vgjl("rdr")

a=pre-daily.nc
griddes $a
ncdump -h $a
#ncpdq -x -v Regular_longitude_latitude -v latitude -v longitude --rdr=t,x,y $a pr-obs.nc
ncpdq --rdr=t,x,y $a pr-obs.nc
process_files(r"pr-obs";tsfile="ts2")
x=lat()
pyjl.npp(x)
da = xr.open_dataset("pr-obs.nc")
#da.sel(time=da.t.dt.year < 2010).mean("x").mean("y").groupby("t.year").sum("t").plot()
#["pr"].mean("longitude").mean("latitude").plot()

display(da)
#lower left of x
da["pr"].isel(t=1).plot()
#xll_corner coordinate of the gif-data
da["pr"].isel(t=101).plot(kind="surface")

o=da["pr"].isel(t=101).squeeze() 
o.plot.surface() 
## get grid information from the netcdf file with R 
# r = rast("D:/remo/genRE/daily_4326.nc")
# ext(r)

#-x, --xcl, --exclude    Extract all variables EXCEPT those specified with -v
#ncpdq -x -v time,x,y --rdr=t,x,y genRE_daily_32632.nc gen.nc
ncpdq -C -x -v lon,lat --rdr=time,x,y genRE_daily_32632.nc gen.nc

cdo -v selyear,2000/2015 genRE_daily_32632.nc ot.nc
vardes ot.nc
griddes ot.nc
descpy ot
ncrename -O -v time,t ot.nc
descpy ot
mv -v ot.nc pr-dly.nc
griddes pr-dly.nc


#ncpdq -C -x -v lon,lat,crs --rdr=time,x,y ot.nc gen.nc #too much memory
#ncrename -O -v time,t gen.nc
xrtsc gen.nc


a=nconly("pr-dly")|>only
open("ts2", "w") do out
    ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n")
    ts_content = replace(ts_content, r"-|T|:" => "\t") #this is for daily... 
    #ts_content = replace(ts_content, r"-|:" => "\t")
    ts_content = replace(ts_content, "T" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "00\t00\t00" => "12\t00")       #like cut -f1-5 ts2
    ts_content = replace(ts_content, "\n" => "", count=1)
    write(out, ts_content)
end
process_files(r"pr-dly";tsfile="ts2")
x=lat()
pyjl.npp(x)

perl -pe 's#D:/remo/genRE#stateini#g' pr-dly.winlist > pr-dly.list 

xrcrds genRE_precipitation_hour_1.1.nc 543972.9 5564373 -1


function calculate_bounding_box(lon::Float64, lat::Float64, cell_size_meters::Float64)
    # Approximation for meters per degree at the equator
    meters_per_degree = 111320.0

    # Calculate degrees for the given cell size
    degrees = cell_size_meters / meters_per_degree

    # Calculate bounding box
    lon_min = lon - degrees / 2
    lon_max = lon + degrees / 2
    lat_min = lat - degrees / 2
    lat_max = lat + degrees / 2

    return lon_min, lon_max, lat_min, lat_max
end

# Example usage:
dlon = 10.0  # Replace with your desired longitude
dlat = 49.0  # Replace with your desired latitude
cell_size_meters = 1200.0

lon_min, lon_max, lat_min, lat_max = 
calculate_bounding_box(dlon, dlat, cell_size_meters)

println("Bounding Box:")
println("Longitude: $lon_min to $lon_max")
println("Latitude: $lat_min to $lat_max")



meters_per_degree = 111320.0
degrees = 1200.0 / meters_per_degree
cb(degrees)
run(`cdo griddes pr-obs.nc`)
pwc()

descpy pr-obs
ncrename -O -v time,t -v longitude,x -v latitude,y pr-obs.nc #err

#ncpdq -C -x -v lon,lat,crs --rdr=time,x,y ot.nc gen.nc #too much memory

griddes daily_4326.nc
descpy daily_4326

cdo chname,daily_4326,pr daily_4326.nc pre-gen.nc
ncrename -O -v time,t -v longitude,x -v latitude,y pre-gen.nc
descpy pre-gen
#ncpdq -C -x -v lon,lat,crs --rdr=time,x,y
a=nconly("pre-gen")|>only
open("ts2", "w") do out
    ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n")
    ts_content = replace(ts_content, r"-|T|:" => "\t")
    ts_content = replace(ts_content, "T" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "11\t00" => "12\t00")       #like cut -f1-5 ts2
    ts_content = replace(ts_content, "\n" => "", count=1)
    write(out, ts_content)
end
process_files(r"pre-gen";tsfile="ts2")
x=lat()
pyjl.npp(x)


cut -f1-5 pre-gen.winlist > tmp
mv -v tmp pre-gen.winlist

function loadcdo()
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/ClimateTools_functions.jl"
    include(pt)
end
loadcdo()
using ClimatePlots
using ClimateTools
using Dates
const tt = ClimateTools
const cm = ClimatePlots
using PyCall
pygui(true)
cd("/mnt/d/remo/genRE")
ncrename -O -v daily_4326,pr daily_4326.nc
#pr = tt.load("pre-daily.nc","pr")
#pr = tt.load("daily_4326.nc","pr")
variables from ./genRE_precipitation_hour_1.1.nc:
['y', 'x', 'lat', 'lon', 'crs', 'pr', 'time']
variables from ./pre-gen.nc:
['t', 'x', 'y', 'crs', 'pr']
pr = tt.load("pr-obs.nc","pr")