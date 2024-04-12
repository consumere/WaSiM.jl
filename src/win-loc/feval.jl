###################wasim output evaluation############################
if length(ARGS) == 0
	println("need args! <file>...")
    exit()
end
#wrkdir = ARGS[1]

wrkdir = "D:/Wasim/regio/coarse-pest/v5/jun/winout/"
#setup()
wrkdir|>cd

#mvwasim2 () 
println("\nmoves all wq, xml and log files to from
c:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05
to current pwd")
ta=pwd()
pt="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05";
println("target dir is $ta");
#@vv "af"

af = filter(x -> occursin(r"wq", x), readdir(pt,join=true))
for i in af
    mv(i,joinpath(ta,basename(i)))
    println(basename(i)," --> ", ta)
end
#@rg "wq"

af = filter(x -> occursin(r"xml", x), readdir(pt,join=true))
for i in af
    mv(i,joinpath(ta,basename(i)))
    println(basename(i)," --> ", ta)
end
af = filter(x -> occursin("modell", x), readdir(pt,join=true))
for i in af
    mv(i,joinpath(ta,basename(i)))
    println(basename(i)," --> ", ta)
end

@gl "mo"

ctl()

ctlfile=raw"D:\Wasim\regio\coarse-pest\controls-win\pst_v5_v2.ctl"
infile = ctlfile
outfile = "route.txt"
routeg(infile, outfile)
#run(`wsl -e perl -i -pe 's/ß/ss/g;s/[\/]/_/g;s/_/-/g;s/[,,]//g;s/\xc4/Ae/g;s/\xd6/Oe/g;s/\xdc/Ue/g;s/\xe4/ae/g;s/\xf6/oe/g;s/\xfc/ue/g;s/\xdf/ss/g' $outfile`)
run(`cat route.txt`)
#readdlm("route.txt")
#readlines("route.txt",keep=false)
x = filter(file -> occursin(r"qges",file), readdir())
cp(first(x),"qfile")
qdf = qall()
##grep the specdis_kmu
pt=pwd()*"/route.txt"
CSV.read(pt,DataFrame,header=false,skipto=6,delim="\t",footerskip=6)
#readdlm(pt)
#specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/saale_spec_ordered.09"
#specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/s200obs"
sim = r"qges"|>glob|>first|>readdf
obs = readf(specfile)
simcols = names(sim)
obscols = names(obs)
df = CSV.read(pt,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
rename!(df,1=>"sim",2=>"obs",3=>"name")

df.name=map(x->replace(x*"-qoutjl",r"#" => ""),df[:,3])  
cnter = 0
for i in eachrow(df)
    cnter += 1
    obsc = i
    try
        dm =innerjoin(
            sim[!,Cols(cnter,end)],
            obs[!,Cols(obsc[2],end)],    
            on=:date)
        onam = obsc[3]
        wawrite(dm,onam)
        println("$onam saved!")
    catch
        onam = obsc[3]
        @warn "merge is empty on $onam ! ..."
    end
end

ftp("Pop")

using PrettyTables
PrettyTables.pretty_table(hcat(qdf,df))
ftp("Unterwe")
#ggofjl()
gofbatch_nosvg()
kgegrep()
rsqgrep()
"tjl"|>globdf
dfp(r"glob")
ftp("Schwein")
dpr(r"Schwein")
dpr(r"qoutjl")

ds = kge_df3()
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)|>PrettyTables.pretty_table

###plots bars of scores...
@time kgeval()
@time nseval()
@time nsevalraw()

r"wind"|>dfp
r"WLF"|>dpr

baryrmean(sim)

#tdiff
td=tdiff()
wawrite(td,"tdiff-jl.txt")
dfp(td)
# tw = copy(td)
# tw = tw[!,Cols(1,2,end)]
# dropmissing!(tw)
# qqp(tw)

tdn = tdifnc()
write("tdiff-jl.nc",tdn;force=true)
rpm(r"tdiff-j",msk=150.0)

waba()
using Images
img=load("waba-jl.png")
r"qbas"|>dfp
r"gwst"|>dfp
r"unterw"i|>baryrsum
"unterw"|>ftp
r"prec"|>baryrsum
regand("rad","mit")|>rplot
z=regand("rad","nc")|>readras
describe(z)

#file remover.
@ncrm

"Schw"|>ftp

r"wind"|>glob|>first|>abspath

wabadf = globdf("wa")|>last|>readf
describe(wabadf)
##annonymous functions
#descr -> subset to minval without datecol -> subset all values >0 
#describe(wabadf)|> r->r[1:nrow(r)-1,Cols(1,3)]|> x-> DataFrames.subset(x,:min => ByRow(>(0.)))
describe(wabadf)|> r->r[1:nrow(r)-1,Cols(1,3)]|> x-> filter(2 => x -> (x) .> 0, x)
#lastbefore(x) = x[end-1]
describe(wabadf)|> r->r[1:nrow(r)-1,Cols(1,4)]|> x-> filter(
    last(names(x)) => x -> (x) .> 0, x)

## subsets rowwise
filter(:KGE => x -> min(x) .> 0, ds)
DataFrames.subset(ds,:NSE => ByRow(>(0.)))
subset(ds,:NSE => ByRow(>(0.)))

###subset and sort by last column
filter(:KGE => x -> min(x) .> 0, ds)|> v -> sort(v,ncol(v),rev=true)
###subset and sort by 2nd column and not rev.
filter(:KGE => x -> min(x) .> 0, ds)|> v -> sort(v,2)
###subset by value and sort by 2nd column and not rev.
filter(2 => x -> min(x) .> .25, ds)|> v -> sort(v,2)

###subset by value and sort by 2nd column and PLOT
filter(2 => x -> min(x) .> .1, ds)|> v -> sort(v,2)|> 
    k -> @df k plot(k.name,cols(propertynames(k)[2:end]), 
    xrotation = 40, xlabel = "", ylabel = "[score]", 
    title = "model efficiecy")


# k = filter(2 => x -> min(x) .> .25, ds)|> v -> sort(v,2)
# @df k plot(k.name,cols(propertynames(k)[2:end]))




filter(2 => x -> min(x) .> .1, ds)|> v -> sort(v,2)|> 
    k -> @df k scatter(k.name, cols(propertynames(k)[2:end]), 
    xrotation = -25, xlabel = "", ylabel = "[score]", 
    title = "model efficiency",
    xticks =  (1:length(k.name),[length(name) > 8 ? SubString(name,1,8) : name for name in k.name]),
    xtickfont = font(6))

short_names = [length(name) > 12 ? SubString(name, 1, 12) : name for name in k.name]
@df k plot(k.name, cols(propertynames(k)[2:end]), 
    xrotation = -25, xlabel = "", ylabel = "[score]", 
    title = "model efficiency",
    xticks = (1:length(k.name), short_names),
    xtickfont = font(6))    


v=raw"D:\Wasim\regio\out\lowres\linlog"
dx=CSV.read(v,DataFrame,header=false,delim="\t")
@vv "corrpl"
@vv "s = Symbol"
#s = Symbol.(filter(x->!occursin(r"year|date",x),names(df)))
dropmissing!(dx)
s = Symbol.(names(dx)[5:10])
@df dx corrplot(cols(s), grid = false)

s = Symbol.(names(dx)[3:end])
@df dx heatmap(cols(s), grid = false)



fl=@gl "nclog"
allstats = CSV.read(fl,DataFrame,header=true)
boxplot(allstats[!,3],yaxis=:log,legend=false,title=names(allstats)[3])

@df allstats boxplot(allstats.filename,allstats.max)
@df allstats plot(allstats.max,xlabel=allstats.filename)

allstats.filename|>PrettyTables.pretty_table


###############rcm200 r3
"D:/Wasim/regio/out/rc200/r3/"|>cd
ds = kge_df3()
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)|>PrettyTables.pretty_table

###plots bars of scores...
@time kgeval()
@time nseval()
@time nsevalraw()
ggofjl()
rx = tdifnc()

@nco "rad"
rplot(r"rad_rcm_1200")
pwd()
dfp(r"glo")
@gl "wab"
@gl "png"

vgctl("rcm200/v2")


using ArchGDAL, Plots

