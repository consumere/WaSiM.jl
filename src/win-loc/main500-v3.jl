raw"D:\Wasim\main\out_main\m500-v3"|>cd
#cluster
# xc = ctl()|>first
# xc = split(xc,"\"")|>first|>y->split(y," ")|>last
# xc = replace(xc, "control"=>"D:/main/control")
infile = raw"D:\Wasim\main\control\mbcl_v3.ctl"
ofl = "route.txt"
routeg(infile, ofl)
DelimitedFiles.readdlm(ofl, ' ', String)|>k->k[6]
specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
#specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu_12"
#sim = r"qgko"|>glob|>first|>readdf
sim = r"qges"|>glob|>first|>readdf
obs = readdf(specfile)

    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
    df.name=map(x->replace(x,r"_>.*" => ""),df.name)
    sort!(df, :sim)


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
#map(x->dpr(x*"-qoutjl")|>display,df.name)
map(x->dpr(x*"-qoutjl"),df.name)

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
waba()
dfp(r"gwst")
td=tdiff()
wawrite(td,"tdiff-jl.txt")

wajs.baryrjs(td|>yrmean)
dfp(td)
rmeq()
@ncrm
#dfs = dfonly(".")
#filter!(x->!occursin(r"^wq",x),dfs)
latx()
dfs=loadalldfs(".")

dxm = mall(dfs)
wawrite(dxm,"mall-jl.txt")
wslpath()|>clipboard
op()

rmdub()

vgjl("-table")

dpr(r"Kleinheu")
wajs.ftpjs(r"Steinbach")

npp(infile)

qdf=wa.qba()
rename!(df,1=>"basin")
kd = innerjoin(qdf,df,on=:basin,makeunique=true)


kd = kd[!,Cols(end,1,2,3)]
sort!(kd, 1)
pretty_table(kd, backend = Val(:text))
# create a string with the LaTeX code
latex_str = sprint(io -> pretty_table(io, kd, backend = Val(:latex)))
dm = dirname(pwd())|>splitpath|>last
write("eval-table-$dm.tex", latex_str)

tx_str = sprint(io -> pretty_table(kd, backend = Val(:text)))
writedf("eval-table-$dm.txt", kd)
latx()


value_to_filter = 1
selected_indices = findall(x -> x == value_to_filter, values)

s=raw"D:/Wasim/main/main_500/mb.ezg"
agheat(s::AbstractString;step=30,lyr=1,msk=0.001)


step=50
lyr=1
lower=0
upper=1000

function agmask(s::AbstractString;step=50,lyr=1,lower=0,upper=1000)
    """
    AG only.
    plots heatmap from raster with ArchGDAL
    upper and lower bound
    """
    if !isfile(s)
        error("file not found!")
    end
    if !endswith(s,".nc")
        @warn("file seems not to be a NetCDF !")
    end
    r = ArchGDAL.readraster(s)
    dx = ArchGDAL.getband(r,lyr)
    #ArchGDAL.getnodatavalue(dx)
    
    if endswith(s,".nc")
        dx = reverse(dx, dims=2)
        dx = reverse(dx, dims=1)
    end
    #findmin(dx)
    if (minimum(dx) <= -9999.0) # || (maximum(dx) >= 9999.0)
        println("extrema: ",join(extrema(dx)," <-> "))
        dmin=minimum(dx)
        @info "
        $dmin set to NaN!
        lower bound: $lower
        upper bound: $upper
        "
        dx = replace(dx, minimum(dx)=>NaN)
        #replace!(dx, minimum(dx)=>missing) #err
        #replace!(dx, minimum(dx)=>NaN)
    end
    
    #heatmap(dx)
    #extrema(dx)
    
    #println("MIN:",minimum(dx),"\nMAX:",maximum(dx))
    
    bitmat = (dx .> lower) .& (dx .< upper)
    # Filter the dx matrix using bitmat
    dx_filtered = dx .* bitmat
    dx_output = Matrix{Float32}(undef, size(dx, 1),size(dx,2))
    dx_output .= NaN
    dx_output[bitmat] .= dx_filtered[bitmat]

    if allequal(dx_output)
        @error("all values are equal!")
        return
    end

    dx = dx_output
    heatmap_plot = heatmap(dx, c=:matter,
                    title=basename(s), 
                    xaxis="", yaxis="");
    step = step
    for i in 1:step:size(dx, 1)
        for j in 1:step:size(dx, 2)
            value = round(dx[i, j]; digits=2)
            color = isnan(value) ? :white : :black
            annotate!(j, i, Plots.text(string(value), 7, color, :center, 
                halign=:center, rotation=-35.0))
        end
    end
    # Show the plot
    display(heatmap_plot)
end

s=raw"D:/Wasim/main/main_500/mb.fzs"
agmask(s;step=90,lyr=1,lower=100,upper=90000)


# r = ArchGDAL.readraster(s)
# dx = ArchGDAL.getband(r,lyr)
# heatmap(dx)

s = @gl "tsoi"
agmask(s;step=30,lyr=5,lower=0.681,upper=2.51)


using PlotlyJS, CSV, DataFrames

df = waread(r"qg")
#df = df[(df.year .== 2007) .& (df.continent .== "Europe") .& (df.pop .> 2e6), :]



dy = wa.yrmean(df)

dy = waread("sobw.m500-v3")|>yrsum
#names(dy)
layout = Layout(uniformtext_minsize=8, uniformtext_mode="hide")
trace = PlotlyJS.bar(
    dy,
    y=:precipitation_588840_5,
    x=:year,
    text=:precipitation_588840_5,
    texttemplate="%{text:.2s}",
    textposition="outside"
)

PlotlyJS.plot(trace, layout)