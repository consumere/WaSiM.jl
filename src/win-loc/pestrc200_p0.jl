raw"D:/Wasim/regio/out/rc200/p0"|>cd
#setup()
####################1 lyr
qgk()
#xf = @gl "xml"
#hrs,min = wa.extract_duration_from_xml(xf)
#println("run took $hrs hrs and $min min...")        
# min

xc = ctl()|>first
xc = split(xc,"\"")|>first
xc = replace(split(xc," ")|>last, "control"=>"D:/Wasim/regio/control")
#xc = split(xc," ")|>last
infile = String(xc)
ofl = "route.txt"
inf="D:/Wasim/regio/out/rc200/x23/cdx/route.txt"
cp(inf,ofl;force=true)
#routeg(infile, ofl)
#sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last

sfn = "specdis_kmu.txt"
sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/"
specfile=joinpath(sfpt,sfn)
obs = readdf(specfile)
#sim = r"qgko"|>glob|>first|>readdf
sim = r"qges"|>glob|>first|>readdf
df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
rename!(df,1=>"sim",2=>"obs",3=>"name")
df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
df.name=map(x->replace(x,r"_>.*" => ""),df.name)
sort!(df, :sim)
#map(x->rm(x),glob("qoutjl"))     #rm all 
@info "
taking prefix of sim and colnumber -> more robust merge with regex of obs!"
outd = []
for i in eachrow(df)
    println(i[1],"->",i[3])
    try
        dm =innerjoin(
            #sim[!,Cols(Regex(string(i[1]),"i"),end)],    
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

#cntcols("outjl")
vz = wa.cntcolv("outjl")
vz = filter(xx->xx[2]==5,vz) #|>x->x[1][1]
#map(x->rm(x[1]),vz)      #remover
#map(x->rm(x),glob(r"qoutjl$"))      #remover
# cd("../")
# kr = kge_rec()
ds=wa.kge_df3()
#gr()
wa.nsx(sort(ds,2))
dfp(r"qbas")
#broadcast(z->display(ftp(z)),r"qoutjl$"|>glob)
broadcast(z->display(wa.ftp(z)),outd)
broadcast(z->display(wa.dpr(z)),outd)
dm = pwd()|>splitpath|>last
@time wa.kgeval()
savefig("kgeplot-$dm.svg")
# @time nseval()
# savefig("nseplot-$dm.svg")
wa.nsevalraw()
savefig("nse-$dm-raw.svg")
ggofjl_nosvg()
#ggofjl()
kgegrep()
kgewrite()
dpr(r"Gem")
dpr(r"Wolf")
facets("sb0")
facets("sb1")
"rad"|>facets
# "rad_rcm_1200.mit.nc"|>facets
# sd = wa.read_soildata(infile)
# names(sd)
# sd[!,end-1]

td=tdiff()
te(td)
wawrite(td,"tdiff-jl.txt")
baryrsum(td)
dfp(td)
#dfl(td)

#dfs = r"so_"|>dfonly|>loadalldfs
#pall(dfs[2:3])
begin  
    qdf=qba()
    rename!(df,1=>"basin")
    kd = innerjoin(qdf,df,on=:basin,makeunique=true)
    kd = kd[!,Cols(end,1,2,3)]
    sort!(kd, 3;rev=true)
    pretty_table(kd, backend = Val(:text))
    # create a string with the LaTeX code
    latex_str = sprint(io -> pretty_table(io, kd, backend = Val(:latex)))
    dm = (pwd())|>splitpath|>last
    write("eval-table-$dm.tex", latex_str)

    tx_str = sprint(io -> pretty_table(kd, backend = Val(:text)))
    writedf("eval-table-$dm.txt", kd)
end

latx()

wslpath()|>clipboard
pwd()|>clipboard
rmeq()
@ncrm
#D:\Wasim\regio\out\rc200\p0: 142.135 MB
#(2003 files)

raw"D:\Wasim\regio\rcm200\pest\v6\out\txts"|>cd
cnt()
zp("cnt")
dsf = wa.readallx(r"qg")
getnames(dsf)

dfs = wa.readall(r"sb1")
sb0 = wa.readall(r"^sb0")
getnames(dfs)

sb0|>first|>wa.cmk.cloudplot
sb0|>last|>wa.cmk.tsp
sb0|>last|>wa.cmk.tspblack
sb0|>last|>dfp
@wajs

sb0|>last|>wa.wajs.dfpjs

sb0|>first|>wa.cmk.tsbar


zp(dpr)|>println

function_str = "$(getnames)"
function_str = string(getnames)
zp(function_str)
wa.zp(fread)
getnames|>wa.zp


baryrsum(dfs[1])
@cmk
wa.cmk.tsp(dfs[1])

bas = wa.readallx(r"qbas")

byr = map(x->yrsum(x),bas )
byr = map(x->select(x,[1,ncol(x)]),byr)
df = mall(byr;xcol=:year)
df = stack(df)
@df df StatsPlots.violin(:year, :value,#group=:variable, 
ylabel="[mm/year]", 
xlabel="Year", legend=false, 
title="Annual discharge")
#ylims=(0, maximum(df.value)+0.1))

dm = map(x->select(x,[1,ncol(x)]),bas )
dm = mall(dm;)
vio(dm)

dfp(bas|>second)
cmk.tsp(bas|>second)

wa.dfp(bas|>second;leg=false)
wa.dfp(r"bas";leg=:outerbottomleft)

zp(ctlg)
@edit wa.ctlg("D:/Wasim/regio/control","x4")

wa.ctlg("D:/Wasim/regio/control","rc200/v4")
vgr("v4")