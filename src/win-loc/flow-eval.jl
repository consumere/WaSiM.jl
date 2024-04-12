####flow-eval########
#wrkdir="/mnt/d/Wasim/regio/out/lowres/c7"
wrkdir="/mnt/d/Wasim/Tanalys/DEM/brendpest/out_v7/spin2"

#"/mnt/d/Wasim/regio/out/lowres/c9/spin3"|>cd
# ddx = "waba-input.wa"|>dfread
# barp(ddx)
cd(wrkdir)
#ctlpath=abspath("../../../control")
ctlpath=abspath("../../control")
cd(ctlpath)
pt=splitdir(wrkdir)|>last
filter(file -> occursin(Regex(pt,"i"),file), readdir())
v = filter(file -> occursin(Regex(pt,"i"),file), readdir())|>first
abspath(v)
cd(wrkdir)
#cp(abspath(v),"route.txt")

v=abspath(v)
cpcmd = `cp -v $v .`
run(cpcmd)

#readdlm("/mnt/d/Wasim/regio/out/lowres/c5/c5-2023/route.txt")
#10:23:11 cp -v "/mnt/d/Wasim/regio/out/lowres/c5/c5-2023/route.txt" .
#rpt="/mnt/d/Wasim/regio/out/lowres/c5/c5-2023/route.txt"
#cp(rpt,"route.txt")
# umlauts route.txt
cd(wrkdir)

pt=pwd()*"/route.txt"
pwd()
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
#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py`
#run(cmd2)
run(`./routecmd.sh`)
ds = kge_df("qout")
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

r"rge"|>dfp
r"global_radiation"|>dfp


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


dfs = readall(".")
for x in dfs
    x=x[!,Not(Cols(r"^Column"))]
    n=DataFrames.metadata(x)|>only|>last|>basename
    println(n," done!")
end
getnames(dfs)
filterplot("rge",dfs)   
filterplot("sb",dfs)
filterplot!("wind",dfs)


df = filterdf("sb",dfs)
s = Symbol.(filter(x->!occursin(r"date|year",x),names(df)))
@df df Plots.plot(:date,cols(s),legend = :outertopright, title=
"titel")

legend(:outertopright)

#"/mnt/d/temp/saale/out_smf200/v5"|>cd

"tsoil"|>facets
r"ts_av"|>dfr
df = r"ts_loc"|>glob|>only|>dfr

x=r"ts"
u=first(glob(x))

#df2 = dfr(r"ts_av")
df2 = readdf(r"ts_av")
dfp(df2)
dfp(df)
plotlyjs()
dfpjs(df)

#filter(x->endswith(x,"qwlf"),readdir())|>first|>linp
df = filter(x->endswith(x,"qwlf"),readdir())|>first|>readdf

dfp(r"qwlf")
df = dfread(r"qwlf")
linp(df)

dfs=rglob(r"SCN"i)|>


dfs=readall(".")


pt="/mnt/d/temp/saale/output/v6/flows.2017"
dfp(pt)

cd("/mnt/d/temp/saale/saale_25/thulba/")
plotlyjs()
"/mnt/d/temp/saale/saale_25/thulba/saale.wit"|>rplot
#r=Raster("/mnt/d/temp/saale/saale_25/thulba/saale.wit")
r=readras("/mnt/d/temp/saale/saale_25/thulba/saale.wit")
#describe(r)
max(r.data)

function stats(r::Raster)
    m = mean(r) # get the mean for each band
    n = minimum(r) # get the minimum for each band
    x = maximum(r) # get the maximum for each band
    d = median(r) # get the median for each band
    s = std(r) # get the standard deviation for each band
    # get the number of missing values for each band
    
    c = try
        parse(Float64, Rasters.missingval(r)) 
        catch
        @warn "No missval in metadat! -set to 0.0"
        c = 0.0
        end
    
    
    arr=[m,n,x,d,s,c]'
    #println("$nm\n",arr)
    df = DataFrame(arr,:auto)
    nm=["mean", "min", "max", "median", "sd", "missval"]
    rename!(df,nm)

    # Matrix(arr) # convert the adjoint to a matrix
    # m = collect(arr) # convert the adjoint to a matrix
    # xc=[
    #     "mean", "min", "max", "median", "sd", "missval",
    # m,n,x,d,s,c]
    
    return(df)
end

r=Raster("/mnt/d/temp/saale/saale_25/thulba/saale.wit",missingval="-9999")

pt="/mnt/d/temp/saale/output/thulba/v0/AnnualTemperature.nc"
r=Raster(pt)
stats(r)    # get summary statistics
describe(r) #bei ncs gehts.

cpl(r)

x=Rasters.rebuild(r;missingval=0)
Plots.contourf(x; c=cgrad(:matter),xlabel="",ylabel="")
Plots.plot(x; c=cgrad(:matter),xlabel="",ylabel="")
plot(r)

df=readdf("/mnt/d/temp/saale/output/thulba/v0/obth-qout")
theplot(df)


"/mnt/d/Wasim/regio/out/rc200/v5/station3"|>cd

#########windows tst. ########################

wrkdir="D:/Wasim/regio/out/rc200/v5/station3/"
cd(wrkdir)
ctlpath=abspath("../../../../control")
cd(ctlpath)
pt=split(wrkdir,"/")|>lastbefore
filter(file -> occursin(Regex(pt,"i"),file), readdir())
v = filter(file -> occursin(Regex(pt,"i"),file), readdir())|>first

#corresponding ctlfile
ctlfile=abspath(v)
cd(wrkdir)
#cp(abspath(v),"route.txt")
cpcmd = `cp -v $ctlfile .`
run(cpcmd)

cpcmd = `routeg $ctlfile > route.txt`
run(cpcmd)
# umlauts route.txt
pt=pwd()*"/route.txt"
pwd()
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
#cmd2 = `python /mnt/c/Users/Public/Documents/Python_Scripts/qgk.py`
#run(cmd2)
run(`./routecmd.sh`)
ds = kge_df("qout")
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)

"bruecken"|>ftp
"sinn"|>ftp
"wolf"|>ftp
"arnst"|>ftp

r"rge"|>dfp
r"global_radiation"|>dfp


r"ett"|>glob
df = "ett"|>glob|>first|>dfr
vgjl("By")

##subset DF by value (all positive vals..)
df2 = filter(1 => x -> x > 0, df) 
ftp(df2)


##--> gofbatch 
gofbatch()
using DelimitedFiles
using Grep
grep_KGE(".")

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


"D:/Wasim/regio/out/rc200/v5/station3/"|>cd
ds=kge_df("qout")
sort(ds, :KGE,rev=true)
sort(ds, :NSE,rev=true)


"Salz"|>ftp

"ts"|>facets
"tsoil"|>facets
"gwst"|>facets
r"gwst"|>dfp
xd = r"gwst"|>globdf
dfr(only(xd))


dx = r"gwst"|>dfr


xd = readdf(r"scnrad"i)
xd = readdf(r"rad"i)
describe(xd)
names(xd)
dfp(xd)
"Bischbrunn" |>ftp
r"Bischbrunn" |>bardf

df = r"Bischbrunn" |>readdf

dfp(df)
ndf = filter(1 => x -> x > 0, df)
dfp(ndf)
dfl(ndf)

"D:/Wasim/regio/out/rc200/v5/station4/"|>cd

"qi"|>facets

r=readras(r"qi")
r=readras(r"qi")
r=Raster("qiflrcm.2012.nc")
rp(r)
describe(r)

"Schweinhof"|>ftp

r"scnrad"i|>dfp
r"flows"i|>bardfm



f="D:/temp/saale/saale_25/thulba_v2/th.dhk"

using GeoArrays
rp3(f)


x="D:/temp/saale/output/thulba/v2/o-qout"
ftp(x)

v="D:/temp/saale/output/thulba/v2/All_HydrologicResponseUnits.nc"
rp(v)
v=nothing

"D:/temp/saale/output/thulba/v2/sobw.v2"|>
"D:/temp/saale/output/thulba/v2/sobw.v2"|>dfp
"D:/temp/saale/output/thulba/v2/sobw.v2"|>bardf

ras=Raster(v)
Plots.surface(
    lookup(ras, X),     
    lookup(ras, Y), 
    ras.data) 


x="D:/temp/saale/output/thulba/v2/vapoth.2017.nc"
rpr(x)
rp(x)


"D:/Wasim/regio/coarse-pest/v5/jun/winout/"|>cd
ncs,dfs = loadall()
dfs = readfall("a")

k = mall(dfs)
dfp(k)
bardf(k)

#waba()
rpm(r"hge";msk=250.0,gt=true)
rpm(r"tsoi";msk=2.0,gt=false)

@vv "brend"