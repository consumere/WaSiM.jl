k="/mnt/d/Wasim/regio/coarse-pest/v5/eobs/"
k=towin(k)
cd(k)
#@time setup()
#wa.mvwasim2()
###########
qgk()
xf = @gl "xml"
infile = ctl2()
##observed runoff! 
grep_with_context("routing_model", infile, 1)
grep_with_context("[SpinUp", infile, 1)
grep_with_context("[groundwater_f", infile, 1)
wa.extract_duration_from_xml(xf)
# run took 0 hrs and 32 min...
#npp(infile)
qbb()
#findlog()
s = @gl "wq"
wqplot(s)
@time facets(r"vapo")
@time facets(r"ts")
dfp(r"qbas")

wa.routeg(infile, "route.txt")
df = wa.dfroute()
outd = wa.pout(infile)
# ad = dfr(r"qges")
# wa.rename_columns!(ad,df)
# dfp(ad)
# dfp(select(ad,Cols(Not(r"C"))))
ds = wa.kge_df3()
wa.nsx(sort(ds,2))
outd = map(waread2,r"qoutjl$"|>glob)
#broadcast(z->display(ftp(z)),r"qoutjl$"|>glob)
broadcast(m->display(ftplin(m)),outd)
#broadcast(z->display(wa.dpr(z)),outd)
dm = pwd()|>splitpath|>last
@time wa.kgeval()
savefig("kgeplot-$dm.svg")
wa.nsevalraw()
savefig("nse-$dm-raw.svg")
#ggofjl_nosvg()
#ggofjl()
using RCall
R"version"
# R"library(ggplot2)"
# @rimport defaultPackages as dfpkg
@rgof

kgegrep()
#kgewrite()
wa.waba()
lx=lat()
irfan(lx)
dfm(r"Capi")
dfm(r"Capi";fun=false)
dfm(r"so_snow_st";fun=false)
pwc()
#bsgs sosubs
dfp(r"rsi"i)
dfp(r"lai"i)

dfr(r"W")|>qplot
facets("sb0")
facets("sb1")
"rad"|>facets
facets("gws")
theme(:dao)
dfm(r"gws")
dfm(r"^gwn";fun=false)

s = @nco "gwst"
rpr(s)
td=tdiff()
te(td)
#wawrite(td,"tdiff-jl.txt")
baryrsum(td)
wa.dfm(td;fun=false)
glob("sf")
rmeq()
wa.rmdub()
@ncrm
#methods(dd)[1].sig
wa.dd(;msg=true)
# 195 files
# D:\Wasim\regio\coarse-pest\v5\eobs: 23.99 MB
wslpath()|>clipboard
pwd()|>clipboard
op()
wa.corrbar(r"qge",r"qgk")
grep_with_context("routing_model", infile, 1)

