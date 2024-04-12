pt="/mnt/d/Wasim/regio/out/rc200/v4/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])

d="qfile"
e = map(x->replace(x*"_qout",r"# " => ""),dout[:,3])
cmd = DataFrame(a="pyx ",b=dout[!,1],c=dout[!,2],d=d,e=e)
cd(dirname(pt))
# writedf("routecmd",cmd)

#oder inside!
#using PyCall
#python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py

a="python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py"

#println(strip(cmd[1,1], '\"'))

cmd = DataFrame(a=strip(a, '\"'),b=dout[!,1],c=dout[!,2],d=d,e=e)
onam="routecmd.sh"
hdr = "#/bin/bash\n"
write(onam, hdr)
#CSV.write(onam,hdr,delim=" ",quotestrings=false)
CSV.write(onam,cmd,
header=false,
#quotestrings=false,
#openquotechar=Char(' '),
#closequotechar=Char(' '),
transform = (col, val) -> something(val, missing),
append=true,
delim=" ")  

#cmd
#Base.unlink(onam)

# my_string = "\"Hello, world!\""
# # trim the quotes
# trimmed_string = strip(my_string, '\"')
# # print the trimmed string
# println(trimmed_string)



# leider geht das nicht mit qoutechar.
#deshalb in bash:
perl -i -nlwe 'tr/"//d; print if length' routecmd.sh 
sfvar qges
cp -v $ff qfile
. routecmd.sh 
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")
@df ds plot(:NSE,:KGE)
#@df ds plot(cols(s),yaxis=:log)
@df ds plot(cols(s))

#s=ds[!,Cols(2:3)]
s = Symbol.(filter(x->!occursin(r"name",x),names(ds)))
#@df ds qqplot(cols(s),qqline = :fit)
@df ds corrplot(cols(s))

vgjl("qq")
vgjl("Symbol.(")
vgjl("Cols")


###nice plots....
using Unitful
plot(Plots.fakedata(100, 10) * u"km", layout=4, 
    palette=[:grays :blues :heat :lightrainbow], 
    bg_inside=[:orange :pink :darkblue :black])

plotlyjs()
l = @layout([a{0.1h};b [c;d e]])
plot(randn(100, 5) * u"km", layout=l, 
    t=[:line :histogram :scatter :steppre :bar], 
        leg=false, ticks=nothing, border=:none)



################v2 test##########
#
pt="/mnt/d/Wasim/regio/out/rc200/v2/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
d="qfile"
e = map(x->replace(x*"_qout",r"# " => ""),dout[:,3])
cmd = DataFrame(a="pyx ",b=dout[!,1],c=dout[!,2],d=d,e=e)
cd(dirname(pt))
# == "/mnt/d/Wasim/regio/out/rc200/v2"|>cd
#oder inside!
#using PyCall
#python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py
a="python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py"
cmd = DataFrame(a=strip(a, '\"'),b=dout[!,1],c=dout[!,2],d=d,e=e)
onam="routecmd.sh"
hdr = "#/bin/bash\n"
write(onam, hdr)
CSV.write(onam,cmd,
header=false,
transform = (col, val) -> something(val, missing),
append=true,
delim=" ")  
# leider geht das nicht mit qoutechar.
#deshalb in bash:
#cmd2 = `head route.txt`
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#cmd2 = `find -maxdepth 1 -name "qges*" -exec cp -v {} qfile \;`
x = rglob("qges")
cp(only(x),"qfile")
run(cmd2)
#cp -v $ff qfile
#. routecmd.sh 
# Execute a shell script file
run(`./routecmd.sh`)
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")
@df ds plot(:NSE,:KGE)

@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)
@df ds plot!(:KGE,
seriestype=:scatter)


Plots.plot(ds.KGE)
Plots.plot!(ds.NSE)

M=Matrix(ds[:,Not(1)])

# @df ds[:,Not(1)] Plots.plot(
#     ds.NSE,ds.KGE,
#     layout=2, 
#     palette=[:blues :heat ], 
#     bg_inside=[:orange :grey ])

using Unitful
Plots.plot(M, layout=2, 
labels = ["NSE","KGE"],
    palette=[:blues :heat ], 
    bg_inside=[:orange :darkblue ])



@df ds plot(1:length(ds.name) ,
cols=[:NSE,:KGE],
#seriestype=:bar,
seriestype=:scatter,
xlabel=:name)


################v0 test##########
#
pt="/mnt/d/Wasim/regio/out/rc200/v0/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
cd(dirname(pt))
a="python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py"
cmd = DataFrame(a=strip(a, '\"'),b=dout[!,1],c=dout[!,2],d=d,e=e)
onam="routecmd.sh"
hdr = "#/bin/bash\n"
write(onam, hdr)
CSV.write(onam,cmd,
header=false,
transform = (col, val) -> something(val, missing),
append=true,
delim=" ")  
# leider geht das nicht mit qoutechar.
#deshalb in bash:
#cmd2 = `head route.txt`
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#cmd2 = `find -maxdepth 1 -name "qges*" -exec cp -v {} qfile \;`
x = rglob("qges")
cp(only(x),"qfile")
#. routecmd.sh 
# Execute a shell script file
run(`./routecmd.sh`)
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")
@df ds plot(:NSE,:KGE)
@df ds qqplot(:NSE,:KGE)

@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)
@df ds plot!(:KGE,
seriestype=:scatter)

Plots.plot(ds.KGE)
Plots.plot!(ds.NSE)

#dy = filter(:year => x -> x > 1950, dy)

filter(:KGE => x -> x > .2, ds)
sort(ds, :KGE,rev=true)
sort(ds, 2,rev=true)


pt="/mnt/d/remo/cordex/eobs/tx_ens_mean_crop.nc"
r = Raster(pt;name=:tx,lazy=true)
r
# xr = r.Where
#z = r[Ti(At(10))]
z = r[tx=2,Ti=2]
Plots.plot(z)
nm = r.dims|>last|>length

r[tx=nm,Ti=nm]|>contourf
r[tx=nm,Ti=nm-365]|>contourf
r[tx=nm,Ti=nm-365]|>contourf

r[tx=1,Ti=1]|>contourf
r[tx=1,Ti=2]|>contourf

xr = r[tx=1,Ti=2]
#xr[X(At(2))]


r.dims|>first|>length
r.dims|>summary


pt="/mnt/d/remo/cordex/eobs/qq_ens_spread_crop.nc"
r = Raster(pt;name=:qq,lazy=true)
nm = r.dims|>last|>length
r[tx=nm,Ti=nm]|>contourf
r[tx=nm,Ti=nm]|>contourf

vgjl("Where")

r[X(Rasters.Where(x -> x >= 9 && x < 10)),Ti=nm]|>Plots.plot
r[X(Rasters.Where(x -> x >= 9 && x < 10)),Y(Rasters.Where(x -> x >= 50 && x < 50.5)),Ti=nm-700]|>Plots.plot
using PlotlyJS
#unload("PlotlyJS")
plotly()
#gr() #geht in pipe... plotly aber nicht!
se = r[X(Rasters.Near(9.2))]
Plots.plot(se,legend=false)

se = r[X(Rasters.Near(9.2)),Y(Rasters.Near(50.00))]
Plots.plot(se,legend=false)

pt="/mnt/d/remo/cordex/eobs/fg_ens_spread_0.1deg_reg_v26.0e.nc"
r = Raster(pt;name=:fg,lazy=true)
se = r[X(Rasters.Near(9.2)),Y(Rasters.Near(50.00))]
Plots.plot(se,legend=false)

using NCDatasets
NCDatasets.cfvariable
ds = NCDataset(pt);
@show cfvariable(ds,"fg", fillvalue = -9999)[:]


##########v5################
pt = "/mnt/d/Wasim/regio/out/rc200/v5/route.txt"

#
pt="/mnt/d/Wasim/regio/out/rc200/v2/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
d="qfile"
e = map(x->replace(x*"_qout",r"# " => ""),dout[:,3])
cmd = DataFrame(a="pyx ",b=dout[!,1],c=dout[!,2],d=d,e=e)
cd(dirname(pt))
# == "/mnt/d/Wasim/regio/out/rc200/v2"|>cd
#oder inside!
#using PyCall
#python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py
a="python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py"
cmd = DataFrame(a=strip(a, '\"'),b=dout[!,1],c=dout[!,2],d=d,e=e)
onam="routecmd.sh"
hdr = "#/bin/bash\n"
write(onam, hdr)
CSV.write(onam,cmd,
header=false,
transform = (col, val) -> something(val, missing),
append=true,
delim=" ")  
# leider geht das nicht mit qoutechar.
#deshalb in bash:
#cmd2 = `head route.txt`
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#cmd2 = `find -maxdepth 1 -name "qges*" -exec cp -v {} qfile \;`
x = rglob("qges")
cp(only(x),"qfile")
run(cmd2)
#cp -v $ff qfile
#. routecmd.sh 
# Execute a shell script file
run(`./routecmd.sh`)
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")
@df ds plot(:NSE,:KGE)

@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)
@df ds plot!(:KGE,
seriestype=:scatter)


Plots.plot(ds.KGE)
Plots.plot!(ds.NSE)

M=Matrix(ds[:,Not(1)])

# @df ds[:,Not(1)] Plots.plot(
#     ds.NSE,ds.KGE,
#     layout=2, 
#     palette=[:blues :heat ], 
#     bg_inside=[:orange :grey ])

using Unitful
Plots.plot(M, layout=2, 
labels = ["NSE","KGE"],
    palette=[:blues :heat ], 
    bg_inside=[:orange :darkblue ])



@df ds plot(1:length(ds.name) ,
cols=[:NSE,:KGE],
#seriestype=:bar,
seriestype=:scatter,
xlabel=:name)


################v0 test##########
#
using CSV,DataFrames,Plots,Dates


pt="/mnt/d/Wasim/regio/out/rc200/v5/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
cd(dirname(pt))
a="python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py"
d="qfile"
e = map(x->replace(x*"_qout",r"# " => ""),dout[:,3])
cmd = DataFrame(a=strip(a, '\"'),b=dout[!,1],c=dout[!,2],d=d,e=e)
onam="routecmd.sh"
hdr = "#/bin/bash\n"
write(onam, hdr)
CSV.write(onam,cmd,
header=false,
transform = (col, val) -> something(val, missing),
append=true,
delim=" ")  
# leider geht das nicht mit qoutechar.
#deshalb in bash:
#cmd2 = `head route.txt`
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#cmd2 = `find -maxdepth 1 -name "qges*" -exec cp -v {} qfile \;`
#x = rglob("qges")
x = filter(file -> occursin(r"qges",file), readdir())
cp(only(x),"qfile")
#. routecmd.sh 
# Execute a shell script file
run(`./routecmd.sh`)
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")
@df ds plot(:NSE,:KGE)
@df ds qqplot(:NSE,:KGE)

@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)
@df ds plot!(:KGE,
seriestype=:scatter)

Plots.plot(ds.KGE)
Plots.plot!(ds.NSE)



lk="/mnt/d/Wasim/regio/out/rc200/v4/route.txt"
lk2="/mnt/d/Wasim/regio/out/rc200/v4/tab1"
# qbb > tab1
# pblank tab1

qbb| perl -nlwe 'tr/"//d;s/\h+/,/g; print if length' > tab1

# fcat route.txt 
# fcat tab1

df = CSV.read(lk,DataFrame,header=false,skipto=8,footerskip=1)
#df2 = CSV.read(lk2,DataFrame,header=false,skipto=4,footerskip=0)
#df2 = CSV.read(lk2,DataFrame,        header=3,        skipto=4,delim=",")
rename!(df,1=>"wasimbasin",2=>"specfile")
df2 = CSV.read(lk2,DataFrame;header=2)
df2 = df2[2:nrow(df2),:]                #drop a line from DF
e = map(x->replace(x*"_qout",r"# " => ""),df[:,3])

A=hcat(e,df[!,1:2])
od = hcat(A,df2)
rename!(od,1=>"name")
vgjl("parse")
#hcat(e,df[!,1:2],df2,makeunique=true)
#map(vcat,[e,df[!,1:2],df2])
vgjl("warn")
od.score = parse.([Int],od.score)
if (od.wasimbasin!=od.score)
    @warn "oh"
    printstyled("NOOOOOO!\n\n",color=:red)
end


dec qsb
qbb > tab2
lk3="/mnt/d/Wasim/regio/out/rc200/v4/tab2"
#df2 = CSV.read(lk3,DataFrame;header=2,delim=" ",downcast=true)

df2 = df2[2:nrow(df2),:]                #drop a line from DF


dfp("/mnt/d/Wasim/regio/out/rc200/v5/rgexrcm.v5.2010")
dfp("/mnt/d/Wasim/regio/out/rc200/v5/temprcm.v5.2010")



##############saale################
"/mnt/d/temp/saale/output/v6"|>cd
using CSV,DataFrames,Plots,Dates

#bash
routeg /mnt/d/temp/saale/control/smf180.v6.ctl > route.txt
pt=pwd()*"/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
cd(dirname(pt))
a="python /mnt/c/Users/Public/Documents/Python_Scripts/mergeQ.py"
d="qfile"
e = map(x->replace(x*"_qout",r"# " => ""),dout[:,3])
cmd = DataFrame(a=strip(a, '\"'),b=dout[!,1],c=dout[!,2],d=d,e=e)
onam="routecmd.sh"
hdr = "#/bin/bash\n"
write(onam, hdr)
CSV.write(onam,cmd,
header=false,
transform = (col, val) -> something(val, missing),
append=true,
delim=" ")  
# leider geht das nicht mit qoutechar.
#deshalb in bash:
#cmd2 = `head route.txt`
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#cmd2 = `find -maxdepth 1 -name "qges*" -exec cp -v {} qfile \;`
#x = rglob("qges")
x = filter(file -> occursin(r"qges",file), readdir())
cp(only(x),"qfile")
#. routecmd.sh 
# Execute a shell script file
run(`./routecmd.sh`)
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")
@df ds plot(:NSE,:KGE)
@df ds qqplot(:NSE,:KGE)
@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)
@df ds plot!(:KGE,
seriestype=:scatter)

Plots.plot(ds.KGE)
Plots.plot!(ds.NSE)

sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)