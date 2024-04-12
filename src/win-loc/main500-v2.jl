"D:/Wasim/main/out_main/2-mb/"|>cd
mvwasim2()
xc = ctl()|>first
xc = split(xc,"\"")|>first|>y->split(y," ")|>last
#xc = replace(xc, "control"=>"D:/Wasim/regio/control")
infile = String(xc)
ofl = "route.txt"
routeg(infile, ofl)
DelimitedFiles.readdlm(ofl, ' ', String)|>k->k[6]
#specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu_12"
#sim = r"qgko"|>glob|>first|>readdf
sim = r"qges"|>glob|>first|>readdf
obs = readdf(specfile)
begin
    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
    df.name=map(x->replace(x,r"_>.*" => ""),df.name)
    sort!(df, :sim)
end

#map(x->rm(x),glob("qoutjl"))     #rm all 
@info "
taking prefix of sim and colnumber -> more robust
merge with regex of obs!
"
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

ds=kge_df3()
nsx(sort(ds,2))
wa.cntcols("outjl")
dpr(r"outjl")
#map(x->dfp(x*"-qoutjl")|>display,df.name)
map(x->dpr(x*"-qoutjl")|>display,df.name)

#dx = qba()

# dm = pwd()|>splitpath|>last
# @time kgeval()
# savefig("kgeplot-$dm.svg")
# @time nseval()
# ggofjl_nosvg()
# kgegrep()

dpr(r"Mi")
findindf(df,"Mi")
facets("sb0")
waba2()
dfp(r"gwst")
td=tdiff()
wawrite(td,"tdiff-jl.txt")
#bardf(td)
dfp(td)