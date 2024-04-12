"D:/Wasim/sinn/out/spin2/v2"|>cd
#@time setup()
infile = "D:/Wasim/regio/coarse-pest/controls-win/pst_v5_v2.ctl"
###########
qgk()
#wa.mvwasim2()
xf = @gl "xml"
hrs,min = wa.extract_duration_from_xml(xf)
message = "run took $hrs hrs and $min min..."
println(message)
clipboard(message)
#run took 0 hrs and 20 min...
begin 
    #npp(infile)
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
    sim = r"qges"|>glob|>first|>readdf
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

end

# cd("../")
# kr = kge_rec()
ds=wa.kge_df3()
#gr()
wa.nsx(sort(ds,2))
dfp(r"qbas")
dfm(r"qbas";leg=false)
#broadcast(z->display(ftp(z)),r"qoutjl$"|>glob)
#broadcast(z->display(wa.ftp(z)),outd)
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
#wa.waba()
#irfan("waba-jl.png")
#latx()
lat()
#rmlat()
dpr(r"Gem")
facets("sb0")
facets("sb1")
"rad"|>facets
facets("gws")
v=@gl "gws"
rpr(v)
v=@gl "sb1"

td=tdiff()

te(td)
wawrite(td,"tdiff-jl.txt")
baryrsum(td)
dfm(td)
dfm(td;fun=false)

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

lat()

glob("ex")
#facets("sfst")
facets("ex")

dfp(r"qbas")

"r"|>glob

rmeq()
@ncrm
#methods(dd)[1].sig
wa.dd(;msg=true)
# 266 files
#  D:\Wasim\regio\coarse-pest\v5\jun\winout: 30.76 MB

wslpath()|>clipboard
pwd()|>clipboard

ser = "^pre"|>nconly
xs = RasterSeries(ser,:t;missingval=0.0)
plot(xs[2:end];xlabel="",ylabel="") #,legend=false)




fn=raw"E:\radolan\pre-dm-2016.nc"
r = Raster(fn)