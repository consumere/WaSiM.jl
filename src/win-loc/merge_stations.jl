
"D:\\remo\\qm\\prec"|>cd
ls()
dfs = readall(r"tsv$")
getnames(dfs)|>grep("M")
xm = mall(dfs)
xm = hcat(xm[!,:date],xm[!,Not(:date)])

@vv "hdr"
v = glob(r"tsv$")
hs=[]
for x in v
    h = CSV.read(x,DataFrame,
        header=1,ntasks=1,transpose=false,limit=4,
        ignorerepeated = true,
        select = [5],
        delim="\t",
        types = Float64,
        silencewarnings=false,
        normalizenames=true)
    push!(hs,h)
end

hs
dfs

hz = reduce(hcat,hs)
data = [
    ["YY", "MM", "DD", "HH"],
    ["YY", "MM", "DD", "HH"],
    ["YY", "MM", "DD", "HH"],
    ["YY", "MM", "DD", "HH"]
]
df = DataFrame(data, :auto)|>permutedims
# Rename the columns
rename!(df, ["YY", "MM", "DD", "HH"])

hz = hcat(df,hz)
propertynames(xm)
rename!(xm, 1=>:date)
dout = copy(xm)
dout.YY = map(x ->year(x),dout.date)
dout.MM = map(x ->month(x),dout.date)
dout.DD = map(x ->day(x),dout.date)
dout[!, "HH"] .= 0

dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]

#hse = hz[!,Not(Cols(r"Meiningen"i))]
#dout = vcat(hse,dout)
dout = vcat(hz,dout)
outfile="pre_rcm85.wa"
CSV.write(outfile, dout, 
transform = (col, val) -> something(val, missing), delim="\t")  


ncol(dout)


#############now to temperature#############
"D:\\remo\\qm\\tas"|>cd
ls()
rgx = r"tsv$"
outfile="tas_rcm85.wa"

dfs = readall(rgx)
#getnames(dfs)|>grep("M")
xm = mall(dfs)
xm = hcat(xm[!,:date],xm[!,Not(:date)]) #reorder

v = glob(rgx)
hs=[]
for x in v
    h = CSV.read(x,DataFrame,
        header=1,ntasks=1,transpose=false,limit=4,
        ignorerepeated = true,
        select = [5],
        delim="\t",
        types = Float64,
        silencewarnings=false,
        normalizenames=true)
    push!(hs,h)
end

hs
hz = reduce(hcat,hs)
data = [
    ["YY", "MM", "DD", "HH"],
    ["YY", "MM", "DD", "HH"],
    ["YY", "MM", "DD", "HH"],
    ["YY", "MM", "DD", "HH"]
]
df = DataFrame(data, :auto)|>permutedims
# Rename the columns
rename!(df, ["YY", "MM", "DD", "HH"])
hz = hcat(df,hz)
rename!(xm, :x1=>:date)
propertynames(xm)
dout = copy(xm)
dout.YY = map(x ->year(x),dout.date)
dout.MM = map(x ->month(x),dout.date)
dout.DD = map(x ->day(x),dout.date)
dout[!, "HH"] .= 0

dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]

#any.(propertynames(hz) in propertynames(dout))
#hse = hz[!,Not(Cols(r"Meiningen"i))]
dout = vcat(hz,dout)

CSV.write(outfile, dout, 
transform = (col, val) -> something(val, missing), delim="\t")  

ncol(dout)


#C:/Users/Public/Documents/Python_Scripts\rsds_rcm.py
vgpy("rsds")
#block ######################################
wrkdir="D:\\remo\\qm\\rsds"
cd(wrkdir)
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
outfile=tv*"_rcm85.wa"
op()
begin
    dfs = readall(rgx)
    #getnames(dfs)|>grep("M")
    xm = mall(dfs)
    xm = hcat(xm[!,:date],xm[!,Not(:date)]) #reorder

    v = glob(rgx)
    hs=[]
    for x in v
        h = CSV.read(x,DataFrame,
            header=1,ntasks=1,transpose=false,limit=4,
            ignorerepeated = true,
            select = [5],
            delim="\t",
            types = Float64,
            silencewarnings=false,
            normalizenames=true)
        push!(hs,h)
    end

    hz = reduce(hcat,hs)
    data = [
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"]
    ]
    df = DataFrame(data, :auto)|>permutedims

    # Rename the columns
    rename!(df, ["YY", "MM", "DD", "HH"])
    hz = hcat(df,hz)
    rename!(xm, :x1=>:date)
    #propertynames(xm)
    dout = copy(xm)
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 0

    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    dout = vcat(hz,dout)

    CSV.write(outfile, dout, 
    transform = (col, val) -> something(val, missing), delim="\t")  

    cls = ncol(dout)
    println("$cls Stations merged and written to $outfile !")
end

##########ps##################################
wrkdir="D:\\remo\\qm\\ps"
cd(wrkdir)
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
outfile=tv*"_rcm85.wa"

#see block above

wrkdir="D:\\remo\\qm\\rh"
cd(wrkdir)
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
outfile=tv*"_rcm85.wa"

wrkdir="D:\\remo\\qm\\wind"
cd(wrkdir)
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
outfile=tv*"_rcm85.wa"

wrkdir="D:\\remo\\qm\\humid"
cd(wrkdir)
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
outfile=tv*"_rcm85.wa"

begin

    dfs = readall(rgx)
    #getnames(dfs)|>grep("M")
    xm = mall(dfs)
    xm = hcat(xm[!,:date],xm[!,Not(:date)]) #reorder

    v = glob(rgx)
    hs=[]
    for x in v
        h = CSV.read(x,DataFrame,
            header=1,ntasks=1,transpose=false,limit=4,
            ignorerepeated = true,
            select = [5],
            delim="\t",
            types = Float64,
            silencewarnings=false,
            normalizenames=true)
        push!(hs,h)
    end

    hz = reduce(hcat,hs)
    data = [
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"]
    ]
    df = DataFrame(data, :auto)|>permutedims

    # Rename the columns
    rename!(df, ["YY", "MM", "DD", "HH"])
    hz = hcat(df,hz)
    rename!(xm, :x1=>:date)
    #propertynames(xm)
    dout = copy(xm)
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 0

    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    dout = vcat(hz,dout)

    CSV.write(outfile, dout, 
    transform = (col, val) -> something(val, missing), delim="\t")  

    cls = ncol(dout)
    println("$cls Stations merged and written to $outfile !")
end

#cmd="bash -c \"cd /mnt/d/remo/qm/prec; ls -1 *.tsv | wc -l\\"
#run(`wsl -d Ubuntu-20.04 -e /mnt/c/Users/Public/Documents/Python_Scripts/heat-codon $outfile`)



#block ######################################
wrkdir="/mnt/d/remo/qm/corgrids/rsds"
cd(wrkdir)
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
outfile=tv*"_cor.wa2"
op()
begin
    dfs = readall(rgx)
    #getnames(dfs)|>grep("M")
    xm = innerjoin(unique.(dfs, :date)..., on = :date, makeunique=true)
    xm = hcat(xm[!,:date],xm[!,Not(:date)]) #reorder

    v = glob(rgx)
    hs=[]
    for x in v
        h = CSV.read(x,DataFrame,
            header=1,ntasks=1,transpose=false,limit=4,
            ignorerepeated = true,
            select = [5],
            delim="\t",
            types = Float64,
            silencewarnings=false,
            normalizenames=true)
        push!(hs,h)
    end

    hz = reduce(hcat,hs)
    data = [
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"]
    ]
    df = DataFrame(data, :auto)|>permutedims

    # Rename the columns
    rename!(df, ["YY", "MM", "DD", "HH"])
    hz = hcat(df,hz)
    rename!(xm, :x1=>:date)
    #propertynames(xm)
    dout = copy(xm)
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 24

    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    dout = vcat(hz,dout)

    CSV.write(outfile, dout, 
    transform = (col, val) -> something(val, missing), delim="\t")  

    cls = ncol(dout)
    println("$cls Stations merged and written to $outfile !")
end



#preci wsl ######################################
wrkdir="/mnt/d/remo/qm/corgrids/pre"
cd(wrkdir)
ssup()
tv = last(splitdir(wrkdir))
rgx = r"tsv$"
#map(rm,glob(rgx)) #rm all
outfile=tv*"_cor.wa"
function mall2(files)
    innerjoin(files..., on = :date, makeunique=true)
end

function mall3(files::Vector{DataFrame})
    df = files[1] 
    for i in 2:length(files)
        @show i
        global df = innerjoin(df, files[i], on = :date, makeunique=true)
    end
    return df
end


    #dfs = readall(rgx)
    files = map(waread,glob(rgx))
    #xm = mall3(dfs)
df = files[1] 
for i in 2:length(files)
    @show i
    #global df = innerjoin(df, files[i], on = :date, makeunique=true)
    global df = hcat(df, select(files[i],1); makeunique=true)
end

xm = df

begin
    xm = hcat(xm[!,:date],xm[!,Not(:date)]) #reorder

    v = glob(rgx)
    hs=[]
    for x in v
        h = CSV.read(x,DataFrame,
            header=1,ntasks=1,
            transpose=false,
            limit=4,
            ignorerepeated = true,
            select = [5],
            delim="\t",
            types = Float64,
            silencewarnings=false,
            normalizenames=true)
        push!(hs,h)
    end

    hz = reduce(hcat,hs)
    data = [
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"]
    ]
    df = DataFrame(data, :auto)|>permutedims

    # Rename the columns
    rename!(df, ["YY", "MM", "DD", "HH"])
    hz = hcat(df,hz)
    rename!(xm, :x1=>:date)
    #propertynames(xm)
    dout = copy(xm)
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 24

    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    dout = vcat(hz,dout)

    CSV.write(outfile, dout, 
    transform = (col, val) -> something(val, missing), delim="\t")  

    cls = ncol(dout)
    println("$cls Stations merged and written to $outfile !")
end

dx = "pre_cor.wa"|>waread
setup()

dfm(dx;fun=yrmean,ann=false,mode=:scatter)
wa.stplot("pre_cor.wa")

###wind ###
pt="/mnt/d/remo/qm/corgrids/wind"
cd(pt)

tv = last(splitdir(pt))
rgx = r"tsv$"
#map(rm,glob(rgx)) #rm all
outfile=tv*"_cor.wa"

files = map(waread,glob(rgx))
df = files[1] 
for i in 2:length(files)
    @show i
    #global df = innerjoin(df, files[i], on = :date, makeunique=true)
    global df = hcat(df, select(files[i],1); makeunique=true)
end

xm = df

begin
    xm = hcat(xm[!,:date],xm[!,Not(:date)]) #reorder

    v = glob(rgx)
    hs=[]
    for x in v
        h = CSV.read(x,DataFrame,
            header=1,ntasks=1,
            transpose=false,
            limit=4,
            ignorerepeated = true,
            select = [5],
            delim="\t",
            types = Float64,
            silencewarnings=false,
            normalizenames=true)
        push!(hs,h)
    end

    hz = reduce(hcat,hs)
    data = [
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"]
    ]
    df = DataFrame(data, :auto)|>permutedims

    # Rename the columns
    rename!(df, ["YY", "MM", "DD", "HH"])
    hz = hcat(df,hz)
    rename!(xm, :x1=>:date)
    #propertynames(xm)
    dout = copy(xm)
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 24

    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    dout = vcat(hz,dout)

    CSV.write(outfile, dout, 
    transform = (col, val) -> something(val, missing), delim="\t")  

    cls = ncol(dout)
    println("$cls Stations merged and written to $outfile !")
end


wa.dfm(outfile;fun=yrmean,ann=false,mode=:scatter)
wa.stplot(outfile)

#@time innerjoin(unique.(files, :date)..., on = :date, makeunique=true);
#first(xm)

size(xm)
xu = unique(xm, keep=:noduplicates)
size(xu)

xu == xm



###tas ###
pt="/mnt/d/remo/qm/corgrids/tas"
cd(pt)
ls()
tv = last(splitdir(pt))
rgx = r"tsv$"
#map(rm,glob(rgx)) #rm all
outfile=tv*"_cor.wa"

files = map(waread,glob(rgx))
df = files[1] 
for i in 2:length(files)
    @show i
    #global df = innerjoin(df, files[i], on = :date, makeunique=true)
    global df = hcat(df, select(files[i],1); makeunique=true)
end

xm = df

begin
    xm = hcat(xm[!,:date],xm[!,Not(:date)]) #reorder

    v = glob(rgx)
    hs=[]
    for x in v
        h = CSV.read(x,DataFrame,
            header=1,ntasks=1,
            transpose=false,
            limit=4,
            ignorerepeated = true,
            select = [5],
            delim="\t",
            types = Float64,
            silencewarnings=false,
            normalizenames=true)
        push!(hs,h)
    end

    hz = reduce(hcat,hs)
    data = [
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"],
        ["YY", "MM", "DD", "HH"]
    ]
    df = DataFrame(data, :auto)|>permutedims

    # Rename the columns
    rename!(df, ["YY", "MM", "DD", "HH"])
    hz = hcat(df,hz)
    rename!(xm, :x1=>:date)
    #propertynames(xm)
    dout = copy(xm)
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= 24

    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    dout = vcat(hz,dout)

    CSV.write(outfile, dout, 
    transform = (col, val) -> something(val, missing), delim="\t")  

    cls = ncol(dout)
    println("$cls Stations merged and written to $outfile !")
end


wa.dfm(outfile;fun=yrmean,ann=false,mode=:scatter)
wa.stplot(outfile)



a = zeros(3)
map!(x -> x * 2, a, [1, 2, 3])
zeros(3,5)


pt="/mnt/d/remo/qm/corgrids/wind/Hofhe_i_sfcWind_cor.tsv"
df = waread(pt)

unique(nonunique(df,:date))

v = nonunique(df,:date)
all(x -> x == v[1], v)

dubs = df[findall(x -> count(==(x), df.date) > 1, df.date),:]

dfp(dubs)
grp = groupby(dubs,:date)

grp2 = groupby(dubs,:date)
v = [size(yx,1) for yx in grp2]|>unique
println("unique counts: $(v)")



@df dubs boxplot(:date, cols(propertynames(dubs)[1]), 
title="Boxplot for Each Group", ylabel="Value", 
xlabel="Date", legend=false)


grp = groupby(dubs,:date)


@df grp[1] violin(:date, cols(propertynames(dubs)[1]), 
title="Boxplot for Each Group", ylabel="Value", 
xlabel="Date", legend=false)


# Create a customized plot for each group
plot()
for group in grp
    @df group violin!(cols(propertynames(dubs)[1]), 
    ylabel="Wind", 
    xlabel="Group", legend=false)
end
title!("Boxplot for Each Group")

function dubplot(pt::String)
    #df::DataFrame
    df = waread2(pt)
    #v = nonunique(df,:date)
    yl=collect(DataFrames.metadata(df))[1][2]|>dirname|>splitpath|>last
    dubs = df[findall(x -> count(==(x), df.date) > 1, df.date),:]
    grp = groupby(dubs,:date)
    v = [size(yx,1) for yx in grp]|>unique
    println("unique counts: $(v)")
    if length(v) == 0
        @info("Warning: no unique count!")
        return
    else
        println("unique counts: $(v)")
    end
    
    plot()
    #i = size(grp,1)
    i = 0;
    for group in grp
        cntz = size(group,1)
        i =+ 1
        # xi = Symbol.(propertynames(dubs)[1])
        # yi = Symbol.(propertynames(dubs)[2])
        # anns = Plots.text(0.5,0.5,
        # "Count $cntz", 
        # halign=:center, valign=:center, pointsize=6)
        # ann = map(x->string.(round(x;sigdigits=3)),dfi.NSE)
        @df group violin!(
            cols(propertynames(dubs)[1]), 
            # annotations = Plots.text("no: $cntz"),
            # annotations = (i,
            # maximum(group[!,1]),
            # "no: $cntz", :top),
            ylabel=yl, 
            xlabel="Group", 
            legend=false)
    end
    for i in 1:size(grp, 1)
        Plots.annotate!(i,
        maximum(grp[i][!,1]),    
        (size(grp[i],1),9,:center,:top,:black))
        # println("R² "*ann[i]*" of Basin 
        # "*(df_A)[i]*" added")
    end
    title!("Bandwith of duplicate values")
end

pt="/mnt/d/remo/qm/corgrids/pre/Hofhe_i_pre_cor.tsv"
pt="/mnt/d/remo/qm/corgrids/rsds/HOF_rsds_cor.tsv"
pt="/mnt/d/remo/qm/prec/Hofheim_pre+rcp85.tsv"
dubplot(pt)

hydromon(pt)

i = 0;
for group in grp; 
    i += 1 
    cntz = size(group,1); b=i;println("$cntz,$b");end


"/mnt/d/Wasim/regio/out/rc200/x22/spin"|>cd
wa.isroute()

lnk="/mnt/d/remo/qm/corgrids/sfcWind-cor.nc_list"


fn=raw"D:/remo/qm/wind/wind_cor.wa"
wa.dubplot(fn)
fn=raw"E:\tmp\rsds_cor.wa"
df = waread(fn)
dfm(df)

fn=raw"E:/tmp/pre_cor24.wa"
df = waread(fn)
sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
obs = waread(sfpt)

ko = selt(obs, r"kiss"i)
ks = selt(df, r"Kiss"i)
qplot(mall(ko,ks))

dz=mall(ko,ks)
dfm(dz)

st = rst.stread(sfpt)
@rsubset st occursin(r"Kiss",:name )

result = @chain df begin 
	@subset :2 .> 20 
	@transform :yr = year.(:date)
	@groupby(:yr) 
	@transform :whmean = mean(:2)
    @select(:date, :whmean, :2, :yr)
end

dfp(result[!,Not(:yr)])
dfp(result[!,Not(3:4)])

ksg = filter(w->occursin(r"Kiss",string.(w.name)),st)
ksg.geometry

phr=dfread("D:/remo/genRE/genRE_daily.wa_ctlp")
dfp(phr)
obsprec = dfread("D:/Wasim/Tanalys/DEM/Input_V2/meteo/prec_ce_bad_kissingen.wa")
xm = mall(obsprec, phr|>dropmissing)
qplot(xm)
xm=reorder_df(xm)
kgs = byear(xm)
wa.pxm(kgs;col=:kge)
wa.pxm(kgs;all=true)
DataFrames.metadata!(xm, "fn" , "prec_ce_bad_kissingen.wa", style=:note);
wa.pxm(kgs;col=:ve)

res = @chain kgs begin
    @subset :2 .> 0.5
    #@select(:2, :3)
    @groupby(:year) 
    @transform :x_mean = mean(:ve)*mean(:kge)
    #@groupby(:ve_mean) 
end

sort!(res, :x_mean; rev=true)