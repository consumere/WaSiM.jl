raw"D:\Wasim\regio\out\vt\s0"|>cd
#eobs local
@time setup()
mvwasim2()
###########
qgk()
xf = @gl "xml"
infile = wa.ctl2()
# infile = replace(infile,r"ctl.*"=>"ctl"
# ,"control"=>"D:/Wasim/regio/control")
## observed runoff! 
grep_with_context("routing_model", infile, 1)
##non - observed runoff! 
grep_with_context("[SpinUp", infile, 1) #yes
grep_with_context("[groundwater_f", infile, 1)
hrs,min = wa.extract_duration_from_xml(xf)
message = "run took $hrs hrs and $min min..."
println(message)
clipboard(message)
qbb()
#findlog()
s = @gl "wq"
wqplot(s)
@time facets(r"ts")
dfp(r"qbas")
#######qgk####################################
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
    @info "taking prefix of sim and colnumber -> more robust merge with regex of obs!"
    outd = []
    for i in eachrow(df)
        println(i[1],"->",i[3])
        try
            dm =innerjoin(
                obs[!,Cols(Regex(i[3],"i"),end)],    
                sim[!,Cols("C"*string(i[1]),end)],
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

map(wcl,glob(r"qoutjl$")) 
#map(rm,glob(r"qoutjl$")) 
qb = dfr(r"qbas")
# df2 = unique(df,:sim)
# wa.rename_columns(qb,df2)
qb2 = @rsubset qb year(:date) .> 2015
baryrsum(qb)
baryrsum(qb2)
dfm(qb2)
dfp(qb2)
####anscheinend brauchts min. 1 Jahr zur Initialize
dfp("qgk")
#map(x->rm(x),glob(r"qoutjl$")) 
ad = dfr(r"qges")
wa.rename_columns(ad,df)
dfp(ad)
dfp(r"qoutjl")
dfp(r"qgko")
wa.dfl(r"qgko")
ro = dfr(r"qgko";)
wa.rename_columns(ro,df)
dfm(ro;fun=yrsum,ann=false)
wa.tline(r"qges"|>fread,:date)
wa.heat(r"qbas")
wa.heat(r"Wolf")
wa.qplot(r"Wolf"|>fread)
baryrsum(ro)
#broadcast(z->display(ftp(z)),r"qoutjl$"|>glob)
broadcast(m->display(ftplin(m)),outd)
broadcast(z->display(wa.dpr(z)),outd)
dm = pwd()|>splitpath|>last
@time wa.kgeval()
savefig("kgeplot-$dm.svg")
wa.nsevalraw()
savefig("nse-$dm-raw.svg")
ggofjl_nosvg()
#ggofjl()
kgegrep()
kgewrite()
wa.waba()
#wa.waba2()
lx=lat()
irfan(lx)
dfm(r"Capi")
dfm(r"Capi";fun=false)
dfm(r"so_snow_st";fun=false)
pwc()
dm = pwd()|>splitpath|>last
dfr(r"Wolf")|>qplot
facets("sb0")
facets("sb1")
"rad"|>facets
facets("gws")
dfm(r"gws")
dfm(r"^gwn";fun=false)
wa.tline(dfr(r"^gwn"),:date)
#theme(:dao)
wa.tline(outd[end],:date)
xm = mall(outd)
wawrite(xm,"mall-jl.wa")
heat(xm)

xm =  wa.mall(glob(r"qges|qbas|qd"))
dropmissing!(xm)
wawrite(xm,"flows-jl.wa")
wa.heat(r"Wolf")
wa.heat(r"qbas")
wa.heat(selt(xm,r"C14"))

pwc()

dfm(r"^sb0";fun=yrmean)
wa.dfrib(r"^sb0";col=:tot_average)
wa.dfrib(r"^qgk")
#wa.dfrib(xm)

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
mx = (r"Wolf")|>dfr
findindf(df,"Wolf")
ftplin(mx)
(r"gws")|>dfm
(r"gws")|>dfp
glob("sf")
rmeq()
#wa.rmdub()
@ncrm

fzplot(;dirpath="D:/Relief_DGMs/FABDEM/wasim/merit")
ezplot(;dirpath="D:/Relief_DGMs/FABDEM/wasim/merit")
@doc fzplot2
ezplot(;dirpath="D:/Relief_DGMs/FABDEM/wasim/merit")

