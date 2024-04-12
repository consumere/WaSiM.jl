using Printf
using Statistics

a="genRE_daily.nc"
# Run cdo showtimestamp $a and redirect the output to "ts"
open("ts2", "w") do out
    #write(out, String(read(pipeline(`cdo griddes $var`, `cat`))))
    ts_content = String(read(pipeline(`cdo showtimestamp $a`)))
    ts_content = replace(ts_content, r"\s+" => "\n") #s/ /\n/g
    #ts_content = replace(ts_content, r"\s+" => "\t")
    ts_content = replace(ts_content, r"-|T|:" => "\t")
    ts_content = replace(ts_content, r"^\s+\S+$" => "")
    ts_content = replace(ts_content, r"^\s+.*$" => "")
    ts_content = replace(ts_content, "\n" => "", count=1)
    #ts_content = replace(ts_content, r"\s+" => "\t")
    write(out, ts_content)
end

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

process_files()

gf="pre-cor.gridfile"
read_grid_value(gf, "xfirst")
rm("pre-cor.winlist2")

a = "f_rsds.nc"
pwc()

@pyimport xarray as xr
# Run cdo xarray redirect the output to "ts"
# Load the NetCDF file into an xarray Dataset
ds = xr.open_dataset(a)
# Access the time attribute
time_values = ds["time"].values
# Print or process the time values as needed
#print(time_values)
dts = [Dates.DateTime.(t.year, t.month, t.day, t.hour, t.minute) for t in time_values]

time_values[1]
t = time_values[1]
#ot = [string.(t.year, t.month, t.day, t.hour, t.minute) for t in time_values]
ot = [join.(hcat(t.year, t.month, t.day, t.hour, t.minute), "\t") for t in time_values]
xo = reduce(vcat, ot)
#ot = [hcat(t.year, t.month, t.day, t.hour, t.minute) for t in time_values]
#ot = join.(string.(ot), "\t")
#convert(Dates.DateTime, ot[1])

#Dates.DateTime(string.(ot[1]|>vec))
#Dates.DateTime(string.(ot[1]), "yyyy mm dd HH MM")
#datetime_matrix = map(s -> Dates.DateTime(s, "yyyy mm dd HH MM"), join(string.(ot[1]), "\t"))




# Convert time values to a string
time_strings = string.(time_values[1])



format_cftime(time_strings)

# Custom function to format cftime values
function format_cftime(cftime_value)
    return Dates.format(Dates.DateTime(cftime_value), "yyyy mm dd HH MM SS")
end

# Extract the date and time components
time_strings = [join(split(format_cftime(value), ' ')) for value in time_values]




######################xr method

import xarray as xr

function xrlist(;xmatch="cor",suffix=".winlist")
    files = filter(x -> contains(x, xmatch) && endswith(x, ".nc"), readdir())
    #files = filter(x -> occursin(r"p+.*nc",x) && endswith(x, ".nc"), readdir())

    if any(x -> occursin(xmatch, x) && endswith(x, suffix), readdir())
        println("Corresponding files already exist. Exiting now!")
        return
    end

    for var in files
        gridfile = string(replace(var,".nc" => ".gridfile"))
        global onam = string(replace(var,".nc" => suffix))

        # open(gridfile, "w") do out
        #     #write(out, String(read(pipeline(`cdo griddes $var`, `cat`))))
        #     write(out, String(read(pipeline(`cdo griddes $var`))))
        # end
        # Load the NetCDF file into an xarray Dataset
        ds = xr.open_dataset(var)

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
                            "YY\tMM\tDD\tHH\t9999",
                            "YY\tMM\tDD\tHH\t$xf",
                            "YY\tMM\tDD\tHH\t$yf",
                            "YY\tMM\tDD\tHH\tListe\n"]

        # open("x_gridfile.txt", "w") do out
        #     write(out, join(gridfile_content, "\n"))
        # end    
        
        open(onam, "w") do out
            write(out, join(gridfile_content, "\n"))
        end

        timevec = ds["time"].values
        #typeof(timevec|>first)
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
            
            link="D:/remo/qm/corgrids/jlcor/$var",
            param="<$param>",
            cnt=(1:length(dts)).-1)
        
        
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
        

cd("/mnt/d/remo/qm/corgrids/jlcor")
xrlist(;xmatch="cor",suffix=".winlist")




#xrlist(;xmatch="pre",suffix=".winlist")


filter(x -> contains(x, "c") && endswith(x, ".nc"), readdir())

ti = ds["time"].values
dts = [Dates.DateTime.(t.year, t.month, t.day, t.hour, t.minute) for t in ti]

da = xr.open_dataset("pre-cor.nc")
pt = da["time"].values

[DateTime(d, dateformat"yyyy-mm-dd") for d ∈ a] # preferred
#a = convert(Array,ti.tolist())
[DateTime(d, dateformat"yyyy-mm-dd") for d ∈ a] # preferred
#DateTime.(ti.tolist(),dateformat"yyyy mm dd HH MM SS")
DateTime.(ti.year)

np_datetime = first(ti)
        
        # @pyimport numpy as np
        # @pyimport datetime
        # datetime.datetime.utcfromtimestamp.((np_datetime - np.datetime64("1970-01-01T00:00:00")) ./ np.timedelta64(1, "s"))

#ncrename -O -d latitude,y -d longitude,x -d time,t pre-cor.nc

ncf = "pre-cor.nc"
ds = xr.open_dataset(ncf)
if "string1" in keys(ds.dims)
    ds = ds.drop_dims("string1")
end
if "Regular_longitude_latitude" in ds.variables
    ds = ds.drop_vars("Regular_longitude_latitude")
end
# Rename dimensions
ds = ds.rename_dims(Dict("latitude" => "y", 
    "longitude" => "x", "time" => "t"))
# Reorder dimensions
ds = ds.transpose("t", "x", "y")
# Update attributes
ds.x.attrs["units"] = "m East of reference point"
ds.y.attrs["units"] = "m North of reference point"
ds.pre.attrs["missing_value"] = -9999.
#outfile = replace(ncf, "_raw.nc" => ".nc")
outpt = ncf
ds.dims
ds.keys()
ds.to_netcdf(outpt)


cp("rh-cor.nc","rh-cor_raw.nc")
ncf = "rh-cor_raw.nc"
ds = xr.open_dataset(ncf)
if "string1" in keys(ds.dims)
    ds = ds.drop_dims("string1")
end
if "Regular_longitude_latitude" in ds.variables
    ds = ds.drop_vars("Regular_longitude_latitude")
end

# Rename dimensions
ds = ds.rename_dims(Dict("lat" => "y", 
    "lon" => "x", "time" => "t"))
# Reorder dimensions
ds = ds.transpose("t", "x", "y")
# Update attributes
ds.x.attrs["units"] = "m East of reference point"
ds.y.attrs["units"] = "m North of reference point"
ds.rh.attrs["missing_value"] = -9999.
#outfile = replace(ncf, "_raw.nc" => ".nc")
outpt = "rh-cor.nc"
ds.to_netcdf(outpt)

#on ubuntu... tas-cor etc. readin geht hier nicht.
cd("/mnt/d/remo/qm/corgrids/wind")
ls()
ob = xr.open_dataset("sfcWind-cor.nc")
k = ob.keys()|>collect|>last
ob[k].mean("t").plot() #this is xr contourf
ob.close()
cd("/mnt/d/remo/qm/corgrids/tas")
ls()
ob = xr.open_dataset("tas-cor.nc")
k = ob.keys()|>collect|>last
ob[k].mean("t").plot() #this is xr contourf
ob[k].mean("x").mean("y").plot() #this is xr contourf
ob.close()
op()

# #bei x32
# exception occured when reading from netCDF file stateini/jlcor/tas-cor.nc. Internal message is: NetCDF: Start+count exceeds dimension bound
# file: netcdfcpp/ncVar.cc  line:1614
# Exception in Modellapplication::run() caught!
# Exception 26 caught at top level!
# Program aborted with exit code 26



pwc()
cd("/mnt/d/remo/genRE")
xrlist(;xmatch="daily",suffix=".winlist")

#pt=/mnt/d/Wasim/main_basin.geojson
(5.238845967940431e6, 5.514078458420032e6)
 (1.206493525543713e6, 1.7250161769856259e6)
gdal_translate -of NetCDF -projwin ulx uly lrx lry genRE_daily.nc output.nc

ulx = 5.238845967940431e6
uly = 1.7250161769856259e6
lrx = 5.514078458420032e6
lry = 1.206493525543713e6
input_file = "genRE_daily.nc"
output_file = "genRE_crop.nc"
run(`gdal_translate -of NetCDF -projwin $ulx $uly $lrx $lry $input_file $output_file`)
run(`cdo griddes $output_file`)

pre = tt.load(output_file,"pr") #mist
#pre = tt.load("genRE_daily.nc","pr") #mist
cm.plot(annualsum(pre))
cm.plot(annualmax(pre))
contourf(pre)

#ubu@wphg066:/mnt/d/remo/genRE$ xrcrds genRE_daily.nc 543972.9 5564373 -1
cd("d:/remo/genRE")
xrcrds genRE_precipitation_hour_1.1.nc 543972.9 5564373 -1
@pj
fn="genRE_daily.wa_ctlp"
df = dfr(fn)
df = pyjl.pyread(fn)
pyjl.pyplot_df(df)
#pyjl.xrp("genRE_daily.nc")
@pyjl
lat()

#hist.interp_like(ref, method="nearest")