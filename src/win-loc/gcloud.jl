raw"L:\04-Phil\Geo1data\prj-Efre-Daten\CSchaefer_TUllman\saale\gcloud\g2"|>cd
# #station 
@time setup()
###########
infile = wa.ctl2()
infile = replace(infile,r"ctl.*"=>"ctl",
"control"=>"D:/Wasim/regio/control")
grep_with_context("routing_model", infile, 1) #1
##non - observed runoff! 
grep_with_context("[SpinUp", infile, 1)
grep_with_context("[groundwater_f", infile, 1)
xf = @gl "xml"
hrs,min = wa.extract_duration_from_xml(xf)
message = "run took $hrs hrs and $min min..."
println(message)
clipboard(message)
#run took 0 hrs and 6 min...
#npp(infile)
qbb()
cdb()
findlog()
qgk()
cd("g1")
#stats to clipboard
begin
    kd=qba()
    strs = sprint(io -> pretty_table(io, kd, header=uppercasefirst.(names(kd)), backend = Val(:text)))
    clipboard(strs)
end
┌───────┬───────────────┬───────────────┬───────────────┬───────────────┐
│ Basin │ LIN. R-SQUARE │ LOG. R-SQUARE │ COEFF.VAR.LIN │ COEFF.VAR.LOG │
├───────┼───────────────┼───────────────┼───────────────┼───────────────┤
│     8 │        0.2918 │       -1.3964 │        0.2919 │       -0.9402 │
│     9 │      -14.8651 │       -3.2109 │      -12.0120 │       -2.9631 │
│    10 │       -0.1555 │       -3.3859 │       -0.1555 │       -2.2044 │
│    11 │       -0.3656 │       -7.3988 │       -0.3568 │       -4.7247 │
│    12 │        0.5733 │       -1.0967 │        0.6683 │       -0.0251 │
│    13 │        0.0859 │       -3.5318 │        0.0910 │       -2.2358 │
│    14 │        0.2387 │        0.0164 │        0.2649 │        0.0212 │
│    15 │       -2.4262 │       -0.8180 │       -1.8829 │       -0.7718 │
│    16 │        0.3961 │       -2.2205 │        0.3973 │       -1.6497 │
└───────┴───────────────┴───────────────┴───────────────┴───────────────┘
dfp(r"qbas")

#######qgk####################################
#hydroeval.objective_functions.kge(simulations, evaluation)
#10.1016/j.jhydrol.2009.08.003

begin 
    ofl = "route.txt"
    routeg(infile, ofl)
    sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
    sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/"
    specfile=joinpath(sfpt,sfn)
    obs = readdf(specfile)

    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
    df.name=map(x->replace(x,r"_>.*" => ""),df.name)
    sort!(df, :sim)
    sim = r"qges"|>glob|>first|>waread
    @info "innerjoin has to be sim -> obs!"
    outd = []
    for i in eachrow(df)
        println(i[1],"->",i[3])
        try
            dm =innerjoin(
                sim[!,Cols("C"*string(i[1]),end)],
                obs[!,Cols(Regex(i[3],"i"),end)],    
                on=:date)
            onam = i[3]*"-qoutjl"
            wawrite(dm,onam)
            println("$onam saved!")
            DataFrames.metadata!(dm, "filename", onam, style=:note);
            push!(outd,dm)
            println(names(dm)," on $onam pushed!")
        catch
            onam = i[3]*"-qoutjl"
            @warn "merge is empty for $onam ! ..."
            continue
        end
    end
    vz = wa.cntcolv("outjl")
    vz = filter(xx->xx[2]==5,vz) #|>x->x[1][1]
end

qb = dfr(r"qbas")
qb2 = @rsubset qb year(:date) .> 2010
baryrsum(qb)
baryrsum(qb2)
dfm(qb2)
dfp(qb2)
ad = dfr(r"qges")
wa.rename_columns!(ad,df)
dfp(ad)
dfp(qb)

broadcast(m->display(wa.ftplin(m)),outd)
broadcast(z->display(wa.dpr(z)),outd)

ks = map(byear,outd)
map(x->Matrix(x[!,Not(:year)])|>maximum,ks)
ks|>cb

typeof(ks)
plot_grouped_metrics(ks;col=:kge)
plot_grouped_metrics(ks;col=:ve)
plot_grouped_metrics(ks;all=true)

#wajs.dfpjs(r"Wolf")

dm = pwd()|>splitpath|>last
wa.kgeval()
@doc wa.pers()
wa.pers(;tofile=true)
npplat()
#run(`powershell.exe -command hyd Wolf`)
@pwrs_str "ls"
pw()

wa.kgeval()
savefig("kgeplot-$dm.svg")
wa.nsevalraw()
savefig("nse-$dm-raw.svg")
ggofjl_nosvg()
#ggofjl()
nsegrep()
kgegrep()
kgewrite()
# pyplot()
# gr()
wa.waba()
lx=lat()
irfan(lx)
dfm(r"Capi")
dfm(r"Capi";fun=false)
dfm(r"so_snow_st";fun=false)
dfr(r"Wolf")|>qplot
#@cmk
facets("sb0")
facets("sb1")

r = readras("sb1_rcm_1400.mit.nc")

@cmk
cmk.mkrheat(r"sb1";missingval=-9999)
cmk.mkrheat(r"sb1";mskval=3000)
cmk.mkrheat(r"sb05";mskval=.4)
cmk.mkrheat(r"stack";mskval=200,layer=2)
Main.cmk.mkrheat(r"^tem")
Main.cmk.mkrheat(r"^tem";mskval=-.12)

dfp("prec")
glob("pre")
@ncrm
Main.cmk.mkrheat(r"^pre";mskval=1)
cl="D:/Wasim/regio/out/rc200/x22/cl4/precrcm_1000.sum.nc"
Main.cmk.mkrheat(cl;mskval=1)

nconly("pre")
Main.cmk.mkrheat(r"^pre+.*1200";mskval=1)
Main.cmk.mkrheat(wa.regand("pre","1100");mskval=1)
cl="D:/Wasim/regio/out/rc200/x22/cl4/precrcm_1100.sum.nc"
cl="D:/Wasim/regio/out/rc200/x22/cl4/precrcm_1200.sum.nc"

Main.cmk.mkrheat(cl;mskval=1)


setup()
wa.rcont("sb1_rcm_1400.mit.nc")
r"rad"|>wa.rcont
r"gws"|>wa.rcont
"rad"|>facets
facets("gws")
r = readras(r"gws")
z = r"gws"|>glob|>last
agheat(z;msk=-20,umask=-0.1)
cmk.mkrheat(z;mskval=-20)
dfm(r"gws")
dfm(r"^gwn";fun=false)
wa.tline(dfr(r"^gwn"),:date)
#theme(:dao)
wa.tline(outd[end],:date)
xm = mall(outd)
wawrite(xm,"mall-jl.txt")
xm =  wa.mall(glob(r"qges|qbas|qd"))
wawrite(xm,"flows-jl.txt")
dropmissing!(xm)
dfm(xm;ann=false)

ds = mall(outd[1:4])
wa.heat(ds)
wa.heat(r"mall-")
wa.heat(r"qbas")

pwc()
dfm(r"^sb0";fun=yrmean)
my = dfr(r"^sb0")|>yrmean
wa.tline!(my;date_col=:year)
wa.tline(my, :year)
dfm(r"^sb1";fun=yrmean)
getm(rpr)
td=tdiff()
te(td)
wawrite(td,"tdiff-jl.txt")
baryrsum(td)
wa.dfm(td;fun=false)
latx()
(r"gws")|>dfp
glob("sf")
rmeq()
#wa.rmdub()
@ncrm
#methods(dd)[1].sig
#@edit wa.dd(;msg=true)
wa.dd(;msg=true)
# 612 files
# D:\Wasim\regio\out\rc200\x34: 122.48 MB

cd("s1")
infile = replace(infile,"x34"=>"x34-s1")
npp(infile)
grep_with_context("routing_model", infile, 1) #1
##non - observed runoff! 
grep_with_context("[SpinUp", infile, 1)
grep_with_context("[groundwater_f", infile, 1)
xf = @gl "xml"
hrs,min = wa.extract_duration_from_xml(xf)
message = "run took $hrs hrs and $min min..."
println(message)
clipboard(message)
#run took 0 hrs and 3 min...
qbb()
#stats to clipboard
begin
    kd=qba()
    strs = sprint(io -> pretty_table(io, kd, header=uppercasefirst.(names(kd)), backend = Val(:text)))
    clipboard(strs)
end
┌───────┬───────────────┬───────────────┬───────────────┬───────────────┐
│ Basin │ LIN. R-SQUARE │ LOG. R-SQUARE │ COEFF.VAR.LIN │ COEFF.VAR.LOG │
├───────┼───────────────┼───────────────┼───────────────┼───────────────┤
│    13 │        0.0105 │       -0.9085 │        0.0744 │       -0.8646 │
│    14 │        0.5108 │       -0.4257 │        0.5138 │       -0.2922 │
│    15 │        0.5934 │       -0.7529 │        0.6075 │       -0.1722 │
│    16 │        0.6097 │       -1.6042 │        0.7237 │       -0.2630 │
│    17 │        0.6766 │       -0.6879 │        0.7179 │       -0.0192 │
│    18 │        0.0532 │       -0.1222 │        0.1465 │       -0.1100 │
│    19 │        0.4528 │       -0.5364 │        0.4605 │       -0.4332 │
│    21 │        0.6061 │        0.2929 │        0.6095 │        0.3506 │
│    22 │        0.7221 │        0.4040 │        0.7318 │        0.4178 │
└───────┴───────────────┴───────────────┴───────────────┴───────────────┘
dfp(r"qbas")
begin 
    ofl = "route.txt"
    routeg(infile, ofl)
    sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
    sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/"
    specfile=joinpath(sfpt,sfn)
    obs = readdf(specfile)

    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
    df.name=map(x->replace(x,r"_>.*" => ""),df.name)
    sort!(df, :sim)
    sim = r"qges"|>glob|>first|>waread
    @info "innerjoin has to be sim -> obs!"
    outd = []
    for i in eachrow(df)
        println(i[1],"->",i[3])
        try
            dm =innerjoin(
                sim[!,Cols("C"*string(i[1]),end)],
                obs[!,Cols(Regex(i[3],"i"),end)],    
                on=:date)
            onam = i[3]*"-qoutjl"
            wawrite(dm,onam)
            println("$onam saved!")
            DataFrames.metadata!(dm, "filename", onam, style=:note);
            push!(outd,dm)
            println(names(dm)," on $onam pushed!")
        catch
            onam = i[3]*"-qoutjl"
            @warn "merge is empty for $onam ! ..."
            continue
        end
    end
    vz = wa.cntcolv("outjl")
    vz = filter(xx->xx[2]==5,vz) #|>x->x[1][1]
end

broadcast(m->display(wa.ftplin(m)),outd)
broadcast(z->display(wa.dpr(z)),outd)
ks = map(byear,outd)
plot_grouped_metrics(ks;col=:kge)
plot_grouped_metrics(ks;col=:ve)
plot_grouped_metrics(ks;all=true)

#wajs.dfpjs(r"Wolf")
dm = pwd()|>splitpath|>last
wa.kgeval()
@doc wa.pers()
wa.pers(;tofile=true)

rmeq()
@ncrm

cd("D:/Wasim/regio/out/rc200/x34/s2/")
#infile = replace(infile,"x34"=>"x34-s2")
#npp(infile)
grep_with_context("routing_model", infile, 1) #1
##non - observed runoff! 
grep_with_context("[SpinUp", infile, 1)
grep_with_context("[groundwater_f", infile, 1)
xf = @gl "xml"
hrs,min = wa.extract_duration_from_xml(xf)
message = "run took $hrs hrs and $min min..."
println(message)
clipboard(message)
#run took 0 hrs and 3 min...
qbb()
#stats to clipboard
begin
    kd=qba()
    strs = sprint(io -> pretty_table(io, kd, header=uppercasefirst.(names(kd)), backend = Val(:text)))
    clipboard(strs)
end
┌───────┬───────────────┬───────────────┬───────────────┬───────────────┐
│ Basin │ LIN. R-SQUARE │ LOG. R-SQUARE │ COEFF.VAR.LIN │ COEFF.VAR.LOG │
├───────┼───────────────┼───────────────┼───────────────┼───────────────┤
│    13 │       -0.3778 │    -4114.7092 │        0.0222 │     -506.8220 │
│    14 │       -0.4897 │    -2579.2164 │        0.0201 │     -264.8959 │
│    15 │       -0.5348 │    -2685.8513 │        0.0283 │     -357.0017 │
│    16 │       -0.7230 │    -2692.5164 │        0.0218 │     -255.6325 │
│    17 │       -0.6244 │    -2184.2005 │        0.0351 │     -272.2024 │
│    18 │       -0.3133 │    -2597.7090 │        0.0095 │     -338.4138 │
│    19 │       -0.3326 │    -3019.1251 │        0.0025 │     -350.6235 │
│    21 │       -0.0566 │       -2.6086 │        0.3807 │        0.8552 │
│    22 │        0.3013 │       -1.4466 │        0.5733 │        0.8379 │
└───────┴───────────────┴───────────────┴───────────────┴───────────────┘

dfp(r"qbas")
begin 
    ofl = "route.txt"
    routeg(infile, ofl)
    sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
    sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/"
    specfile=joinpath(sfpt,sfn)
    obs = readdf(specfile)

    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
    df.name=map(x->replace(x,r"_>.*" => ""),df.name)
    sort!(df, :sim)
    sim = r"qges"|>glob|>first|>waread
    @info "innerjoin has to be sim -> obs!"
    outd = []
    for i in eachrow(df)
        println(i[1],"->",i[3])
        try
            dm =innerjoin(
                sim[!,Cols("C"*string(i[1]),end)],
                obs[!,Cols(Regex(i[3],"i"),end)],    
                on=:date)
            onam = i[3]*"-qoutjl"
            wawrite(dm,onam)
            println("$onam saved!")
            DataFrames.metadata!(dm, "filename", onam, style=:note);
            push!(outd,dm)
            println(names(dm)," on $onam pushed!")
        catch
            onam = i[3]*"-qoutjl"
            @warn "merge is empty for $onam ! ..."
            continue
        end
    end
    vz = wa.cntcolv("outjl")
    vz = filter(xx->xx[2]==5,vz) #|>x->x[1][1]
end
@edit wa.ftplin(df)
broadcast(m->display(wa.ftplin(m)),outd)
broadcast(m->display(wa.ftp(m)),outd)
broadcast(z->display(wa.dpr(z)),outd)
ks = map(byear,outd)
plot_grouped_metrics(ks;col=:kge)
plot_grouped_metrics(ks;col=:ve)
@edit plot_grouped_metrics(ks;all=true)

#wajs.dfpjs(r"Wolf")
dm = pwd()|>splitpath|>last
wa.kgeval()
@doc wa.pers()
wa.pers(;tofile=true)

rmeq()
@ncrm



wa.kgeval()
savefig("kgeplot-$dm.svg")
wa.nsevalraw()
savefig("nse-$dm-raw.svg")
ggofjl_nosvg()
#ggofjl()
nsegrep()
kgegrep()
kgewrite()
# pyplot()
# gr()
wa.waba()
lx=lat()
irfan(lx)
dfm(r"Capi")
dfm(r"Capi";fun=false)
dfm(r"so_snow_st";fun=false)
dfr(r"Wolf")|>qplot
#@cmk
facets("sb0")
facets("sb1")

r = readras("sb1_rcm_1400.mit.nc")

@cmk
cmk.mkrheat(r"sb1";missingval=-9999)
cmk.mkrheat(r"sb1";mskval=3000)
cmk.mkrheat(r"sb05";mskval=.4)
cmk.mkrheat(r"stack";mskval=200,layer=2)
Main.cmk.mkrheat(r"^tem")
Main.cmk.mkrheat(r"^tem";mskval=-.12)

dfp("prec")
glob("pre")
@ncrm
Main.cmk.mkrheat(r"^pre";mskval=1)
cl="D:/Wasim/regio/out/rc200/x22/cl4/precrcm_1000.sum.nc"
Main.cmk.mkrheat(cl;mskval=1)

Main.cmk.mkrheat(r"gwst";msk=false)
agheat("gwstrcm_1200.mit.nc";msk=-50,umask=-0.1)
"D:/Wasim/regio/rcm200/v11/rcm.dhk"|>cmk.mkrheat
fn="D:/Wasim/regio/rcm200/v11/rcm.s01"
agheat(fn;msk=0.00005,umask=0.0045,roundto=6)

cd("D:/Wasim/regio/rcm200/v11/")
op()

#wit calcs.