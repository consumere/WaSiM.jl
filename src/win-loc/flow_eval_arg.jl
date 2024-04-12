##############saale################

wrkdir="/mnt/d/Wasim/regio/out/c4/c4_rev"
cd(wrkdir)
#routeg ../../../control/rcm-c4-rev.ctl > route.txt
using CSV,DataFrames,Plots,Dates
#bash
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
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#x = glob("qges")
x = filter(file -> occursin(r"qges",file), readdir())

##headers were broken, so bashfunc:
# reghdrs
cp(only(x),"qfile")

#dy = waread("qfile")
#dfp("qfile")
#. routecmd.sh 
# Execute a shell script file
run(`./routecmd.sh`)
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")

sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

##some plotting...
@df ds plot(:NSE,:KGE)
@df ds qqplot(:NSE,:KGE)
describe(ds)

@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)

@df ds plot!(:KGE,
seriestype=:scatter)

Plots.plot(ds.KGE)
Plots.plot!(ds.NSE)

##doublecheck
fgof2="cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R"
gofyr='cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof_yr.R"'

#run(fgof2*" Wolfsmuenster_qout")
run(`cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R Arnstein_qout`)
run(`grep -rIHn -A1 -B3 "KGE"`)
ds
##year
run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof_yr.R" Wolfsmuenster_qout 2`)
#run(`grep -rIH -A2 -B2 --include=*_output.txt VE`)




"/mnt/d/Wasim/regio/out/cl"|>cd
m=rglob(r"qout")
m=filter(x->!occursin(r"xml|txt|mon|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",x),m)
ds = map(x->readdf(x),m)
a = map(x->kge2(x),ds)



function kge_df(ext::Regex)
    """
    should be recursive
    """
    m=rglob(ext)
    files=filter(x->!occursin(r"xml|txt|mon|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",x),m)
    v = []
    for file in files
        if isfile(file)
            # dd = CSV.read(file,DataFrame,missingstring="-9999",delim="\t")
            # observed  = dd[:,5]
            # simulated = dd[:,6]
            dd  =   readdf(file)
            observed  = dd[:,end-2]
            simulated = dd[:,end-1]
            kge_value = kge2(observed, simulated)
            nse_value = nse(observed, simulated)
            nm = basename(file)
            println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
            printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
            push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
            v = DataFrame(v)
        end
    end
    return(v)
end

ds=kge_df(r"qout")
x = readdf(m)





##############rc200################

wrkdir="/mnt/d/Wasim/regio/out/rc200/r2"
cd(wrkdir)
#routeg ../../../control/rcm200_r2.ctl > route.txt
using CSV,DataFrames,Plots,Dates
#bash
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
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#x = glob("qges")
x = filter(file -> occursin(r"qges",file), readdir())
x
cp(first(x),"qfile")

r"rge"|>dfp
"SCNTEMPrcm.2017.nc"  |>facets
"SCNLIQOUT_rcm.r2.2017" |>dfp

using PlotlyJS
"SCNLIQOUT_rcm.r2.2017" |>dfpjs



f=glob("qgk")|>first
#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qsa.py $f`
#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py`
cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py > teval`
# run(cmd2 > "temp") #geht net
#output = readlines()
run(cmd2)


using DataFrames

function qall()
    files = rglob("qgko")
    for file in files
        # Load the file into a DataFrame
        #x = "qgkofab.m6.2010"
        x = file
        try
            df = DataFrame(CSV.File(x, header=false, 
                                delim="\t",
                                ignorerepeated=true,
                                types = String
                                )) 
                                #types=[String, String, String])                                )
            println(x)
            pattern = r"^[LIN. R]|^[LOG. R]|^CO"
            #m=match(r".*[.]",s)
            #outfile = string(m.match,"png")
            #string.(df[i,:])
            # first(eachrow(df[!,1]))
            # m=[]
            # for i in eachrow(df[!,1])
            #     k=i
            #     n=(filter(line -> occursin(pattern,line),k))
            #     push!(m,m)
            # end
            # m
            mask = [occursin(pattern, df[i, 1]) for i in 1:nrow(df)]
            dx = df[mask, :]
            dx = permutedims(dx) |>dropmissing
            
            #mapcols!(x -> parse(Float64, x), dx)
            #mapcols(x -> x.^2, dx)

            #basins = copy(df[1,5:end])
            #AsTable(basins)
            basins = []
            for i in copy(df[1,5:end])
                push!(basins,i)
            end
            #size(basins)
            insert!(basins, 1, "score")
            insert!(basins, 2, "timestep")
            dx[!, "basin"] = basins


            cn = (dx[1,:])
            rename!(dx, 1 => cn[1], 2 => cn[2], 3 => cn[3], 4 => cn[4])
            dout = dx[3:end,:]
            for col in names(dout)[1:end-1]
                dout[!, col] = parse.(Float64, replace.(dout[!, col], "," => ""))
            end          
            
            #mapcols!(x -> parse(Float64, x), dout)
            #dout = hcat(dx[!,Cols(r"bas")],dx[:,Not(Cols(r"bas"))])
            dout.basin = map(x -> parse(Int, x), dout.basin)
            dout = hcat(dout[:,Cols("basin")],dout[:,Not(Cols(r"bas"))])
            return(dout)
        catch
            @warn("error! ")
            # Skip files that can't be loaded as a DataFrame
            continue
        end
    end
end

#f=glob("qgk")|>first
#xd = qall(f)

ds = qall()
typeof(ds)

s = Symbol.(filter(x->!occursin(r"bas",x),names(ds)))
@df ds plot(cols(s),:basin)
@df ds plot(:basin,cols(s))
@df ds plot(:basin,cols(s),yaxis=:log)
@df ds plot(:basin,cols(s[2]),label=string.(s[2]))

#dy = waread("qfile")
#dfp("qfile")
#. routecmd.sh 
# Execute a shell script file
run(`./routecmd.sh`)
#jlkge qout|sort -nr|tee kge.log
#kge_read(pwd(),"qout")
ds = kge_df("qout")

sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

##some plotting...
@df ds plot(:NSE,:KGE)

@df ds qqplot(:NSE,:KGE)

describe(ds)

@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)

@df ds plot!(:KGE,
seriestype=:scatter)

Plots.plot(ds.KGE)
Plots.plot!(ds.NSE)

##doublecheck
fgof2="cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R"
gofyr='cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof_yr.R"'

#run(fgof2*" Wolfsmuenster_qout")
run(`cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R Arnstein_qout`)
run(`grep -rIHn -A1 -B3 "KGE"`)
ds
##year
run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof_yr.R" Wolfsmuenster_qout 2`)
#run(`grep -rIH -A2 -B2 --include=*_output.txt VE`)




"/mnt/d/Wasim/regio/out/cl"|>cd
m=rglob(r"qout")
m=filter(x->!occursin(r"xml|txt|mon|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",x),m)
ds = map(x->readdf(x),m)
a = map(x->kge2(x),ds)



function kge_df(ext::Regex)
    """
    should be recursive
    """
    m=rglob(ext)
    files=filter(x->!occursin(r"xml|txt|mon|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|zip|tar",x),m)
    v = []
    for file in files
        if isfile(file)
            # dd = CSV.read(file,DataFrame,missingstring="-9999",delim="\t")
            # observed  = dd[:,5]
            # simulated = dd[:,6]
            dd  =   readdf(file)
            observed  = dd[:,end-2]
            simulated = dd[:,end-1]
            kge_value = kge2(observed, simulated)
            nse_value = nse(observed, simulated)
            nm = basename(file)
            println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
            printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
            push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
            v = DataFrame(v)
        end
    end
    return(v)
end

ds=kge_df(r"qout")
x = readdf(m)


"/mnt/d/Wasim/regio/out/rc200/v3"|>cd
"/mnt/d/Wasim/Goldbach/revision/fab150/v6"|>cd
"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m7"|>cd
"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m6"|>cd
"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m7/old"|>cd
ds = qall()
ds
#typeof(ds)
s = Symbol.(filter(x->!occursin(r"bas",x),names(ds)))
@df ds plot(cols(s),:basin)
@df ds plot(:basin,cols(s))
@df ds plot(:basin,cols(s),yaxis=:log)
@df ds plot(:basin,cols(s[2]),label=string.(s[2]))


@df ds corrplot(cols(s))


#label=["Simulated" "Observed"], xlabel="Date", ylabel="Value")

ds = qall_num()

v = (names(ds[!,Not(1,4,5)]))
a = reshape(v, 1, 2)

Plots.plot(ds.basin,[ds[!,2], ds[!,3]]) , 
label=a, xlabel="Date", ylabel="[mm/day]",legend = :topleft)

r2 = round(cor(df[!,2], df[!,1])^2, digits=2)
annotate!(last(df.date), 0.85*maximum(df[!,1]),
text("R² = $r2", 10, :black, :right))




############v5 ################

wrkdir="/mnt/d/Wasim/regio/out/rc200/v5/evap"
cd(wrkdir)
#routeg ../../../control/rcm200_r2.ctl > route.txt
#/mnt/d/Wasim/regio/control/rcm200_v5_win.ctl

using CSV,DataFrames,Plots,Dates
#bash
pt=pwd()*"/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
# cd(dirname(pt))
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
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
cmd2 = `tail routecmd.sh`
run(cmd2)
#x = glob("qges")
x = filter(file -> occursin(r"qges",file), readdir())
x
cp(first(x),"qfile")

r"rge"|>dfp
r"temp"  |>glob

#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py`
#run(cmd2)

run(`./routecmd.sh`)

ds = kge_df("qout")
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

##some plotting...
@df ds plot(:NSE,:KGE)
@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)

##doublecheck
fgof2=`cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R  Arnstein_qout`
run(fgof2)

gofyr='cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof_yr.R"'

#run(fgof2*" Wolfsmuenster_qout")
run(`cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R Arnstein_qout`)
run(`grep -rIHn -A1 -B3 "KGE"`)
ds
##year
run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof_yr.R" Wolfsmuenster_qout 2`)


##
#afnse qout hydt > evallog
#gofbatch #!


# fdx -t f "glo|regex" 

# fdx -d 1 -t f "glo|rgex" -E "*Asp*" -E "*nc"
# fdx -d 1 -t f "glo|rgex" -E ".nc"
# #find . -maxdepth 1 -type f -name "global_*" -o -name "rgex*" -not -regex '.*.\(png\|svg\|html\|ftz\|ftz_0\|list\|nc\|xml\|sh\|grd\|yrly\)$' -print

# radc () 
# { 
#     fdx -d 1 -t f "global*|rgex*" -E "*Asp*" -E "*nc" | while read f; do
#         if [ ! -f "$f" ]; then
#             echo "File "$1" does not exist";
#         else
#             printf '%s \n%s\n' "$(tff "${f}")";
#         fi;
#     done
# }

wrkdir="/mnt/d/Wasim/regio/out/rc200/v5/station"
cd(wrkdir)
#bash
pt=pwd()*"/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
# cd(dirname(pt))
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
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
x = filter(file -> occursin(r"qges",file), readdir())
x
cp(first(x),"qfile")
r"rge"|>dfp
r"global_radiation"|>dfp
#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py`
#run(cmd2)
run(`./routecmd.sh`)
ds = kge_df("qout")
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

##some plotting...
@df ds plot(:NSE,:KGE)
@df ds plot(:NSE ,
seriestype=:bar,
xlabel=:name)

##doublecheck
#gofbatch
# fgof2=`cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R  Arnstein_qout`
# run(fgof2)

#gofyr='cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof_yr.R"'
#run(fgof2*" Wolfsmuenster_qout")
#run(`cmd.exe /c Rscript D:/Fernerkundungsdaten/Klassifikation/R-Sessions/fgof2.R Arnstein_qout`)
run(`grep -rIHn -A1 -B3 "KGE"`)

#rename -v 's/´┐¢/ss/' *

r"Bad_Brueckenau"|>dpr

pt="/mnt/d/Wasim/regio/out/rc200/v5/station/qd__rcm.2012.nc"
cd(dirname(pt))
facets(pt|>basename)
r"tsoil"|>glob
r"tsoil"|>facets
r"avg"|>dfp
r"loc"|>dfp
df = r"loc"|>readdf


using Plots, CSV, DataFrames

x="/mnt/d/Wasim/regio/out/rc200/v5/station/Wechterswinkel_qout"


function theplot(x::AbstractString)
    df = DataFrame(CSV.File(x))
    df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
    ndf = df[!,Not(1:4)]
    rename!(ndf, [:Observed, :Simulated, :Date])
    overall_pearson_r = cor(ndf[!, :Observed], ndf[!, :Simulated])
    r2 = overall_pearson_r^2
    nse_score = nse(ndf[!, :Observed], ndf[!, :Simulated])
    kge_score = kge2(ndf[!, :Observed], ndf[!, :Simulated])
    #ti = "Time Series of $(uppercase(first(split(basename(x), '-'))))"
    ti = first(split(basename(x),"_"))
    #subs = "Pearson r: $(round(overall_pearson_r, digits=2))\nPearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    subs = "Pearson R²: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))"
    #p = plot(title=[ti, subs], ylabel="[unit/day]", xlabel="modeled time", yscale=:log, legend=:topleft)
    p = plot(title=ti, ylabel="[mm/day]", xlabel="modeled time", 
        yscale=:log, legend=:topleft)
    plot!(p, ndf[!, :Date], ndf[!, :Simulated], line=:dash, color=:blue, label="Modeled")
    plot!(p, ndf[!, :Date], ndf[!, :Observed], color=:red, label="Observed")
    annotate!(
    #nrow(ndf), 0.95*maximum(ndf.Observed),
    :bottomright,
    text("$subs", 10, :black, :right))
    return p
end


"/mnt/d/Wasim/regio/out/rc200/v5/station/"|>cd
x=glob("qout")
x=first(x)
theplot(x)

ds = qall()


function evp(ext::AbstractString)
    files = readdir(".")
    v=[]
    for file in files
        file_path = joinpath(pwd(), file)
        if isfile(file_path) && endswith(file, ext)
           println("eval on $file")
           pl=theplot(file);
           push!(v,pl)	   
        end
    end
    return(v)
end

pts = evp("qout")
display(pts[5])
display(pts[9])
for file in filter(x->endswith(x,"qout"),readdir())
    pl=nothing   
    pl=theplot(file);
    Plots.savefig(pl,file*"_eval.png")
end


v=filter(x->endswith(x,"qout"),readdir())
theplot(v[end-1])
theplot(v[end])
theplot(v[3])
theplot(first(filter(x->occursin(r"Unter",x),v)))

theplot(first(filter(x->occursin(Regex("lach","i"),x),v)))

print(v)

theplot(first(
    filter(
    x->occursin(Regex("hafen","i"),x),filter(x->endswith(x,"qout"),readdir()))))


ftp(z::AbstractString) = theplot(first(filter(x->occursin(Regex(z,"i"),x),filter(x->endswith(x,"qout"),readdir()))))
ftp("arn")
ftp("steinbach")

"Roemershofe"|>ftp

dfp(r"td")
dfp(r"sb1")
rp(r"te")
rp(r"rad")

"/mnt/d/Wasim/regio/out/rc200/v5/station2"|>cd
pt=pwd()*"/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])
# cd(dirname(pt))
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
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
x = filter(file -> occursin(r"qges",file), readdir())
x
cp(first(x),"qfile")
r"rge"|>dfp
r"global_radiation"|>dfp
#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py`
#run(cmd2)
run(`./routecmd.sh`)
ds = kge_df("qout")
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

"Roemershofe"|>ftp
"Brueckenau"|>ftp
"wechters"|>ftp

fd()
r"gwst"|>dfp
r"wab"|>glob


"waba-input.wa"|>vio

using PyPlot
pyplot()
"waba-input.wa"|>dfp

plotlyjs()
"Brueckenau"|>ftp
r"Brueckenau"|>dfpjs

ddx = "waba-input.wa"|>dfread


##some plotting...
@df ds plot(:NSE,:KGE)
@df ds plot(:NSE ,seriestype=:bar,xlabel=:name)


using DelimitedFiles
f="wurzrcm.v5.2012"
data = String(read(f, String)) 
# assuming the file contains text


####lowres########
"/mnt/d/Wasim/regio/out/lowres/c9/spin3"|>cd

ddx = "waba-input.wa"|>dfread
barp(ddx)


# umlauts route.txt
pt=pwd()*"/route.txt"
df = CSV.read(pt,DataFrame,header=false,skipto=8,footerskip=1)
dsub = df[!,1:2]
dout = hcat(dsub .+ 3,df[!,3])

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
cmd2 = `perl -i -nlwe 'tr/"//d; print if length;warn"done!"' routecmd.sh`
run(cmd2)
x = filter(file -> occursin(r"qges",file), readdir())
x
cp(first(x),"qfile")
r"rge"|>dfp
r"global_radiation"|>dfp
#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py`
#run(cmd2)
run(`./routecmd.sh`)
ds = kge_df("qout")
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

##--> gofbatch 

plot(ds[!,1],ds[!,end],seriestype=:bar,
xlabel = "", ylabel = "KGE", legend = false,
xrotation = 45)

bar(ds.name, ds.NSE, xlabel = "Name", ylabel = "NSE", legend = false, 
    title = "NSE vs Name", xrotation = 45, fmt = :png, size = (800, 600), 
    bar_width = 0.6, color = :blue)

scatter(ds.NSE, ds.KGE,
    xlabel = "NSE",
    ylabel = "KGE",
    title = "NSE vs KGE",
    legend = false)




ds[!,Not(:name)]|>Matrix|>corrplot

# M = Matrix(ds[!,Not(:name)])
#corrplot(M)
# M = Matrix(ds)
#corrplot(M[:,Not(1)])

s = Symbol.(filter(x->!occursin(r"name",x),names(ds)))
@df ds corrplot(cols(s), grid = false)

"sinn"|>ftp
"wolf"|>ftp
"arnst"|>ftp

vgctl("höheren")

facets(r"tsoi")




function gofbatch()
    println("batch R Script for GOF")
    arr = filter(x -> isfile(x) && endswith(x, "_qout") && !occursin(r"\.(png|svg|txt|html|ftz|ftz_0|list|nc|xml|sh|grd|yrly)$", x), readdir())
    for i in arr
        println(i)
    end
    if length(arr) < 1
        println("No match!\nneed qout files")
        return 1
    end
    for i in arr
        x = basename(i)
        run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof3.R" $x`)
    end
end

gofbatch()


using Grep

output = DelimitedFiles.readdlm(file,'\t', String)
Grep.grep(r"KGE.*[0-9].[3-9]",output)


function grep_KGE(path::AbstractString)
    #@printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), readdir(path))
        #output = read(file, String)
        output = DelimitedFiles.readdlm(file,'\t', String)
        #match = occursin(r"KGE.*[0-9].[3-9]", output)
        match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
        if !isempty(match)
            #@printf("%s: %s\n", file,match)
            #@printf("%s:", first(split(file,"_qout")))
            fn = first(split(file,"_qout"))
            for line in match
                #@printf("\t%s\n", line)
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
            end
        end
    end
end

grep_KGE(".")

line="KGE       -3.51"
strip(line)
join(split(line), " ") 

function search_for_KGE_output_files(root_path::String)
    files = filter(isfile, readdir(root_path))
    output_files = filter(file -> endswith(file, "_output.txt"), files)
    result = ""
    for file in output_files
        contents = read(file, String)
        if occursin(r"KGE.*[0-9].[3-9]", contents)
            result *= string(file, "\n")
        end
    end
    println(result)
end
search_for_KGE_output_files(".")


"smf180"|>vgr


cd("/mnt/d/Wasim/regio/out/lowres/c7")