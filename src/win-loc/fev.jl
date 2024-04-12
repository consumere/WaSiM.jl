###################wasim output evaluation############################
if length(ARGS) == 0
	println("need args! <file>...")
    exit()
end
#wrkdir = ARGS[1]
#using Base.Threads
#nthreads()
setup()
wrkdir = "D:/temp/saale/out_smf200/lnx/v1/"
wrkdir|>cd
# xc = ctl()
# xc = first(xc)
# xc = split(xc,"\"")|>first|>z->split(z," ")|>last
# infile = joinpath(raw"D:\Wasim\regio",xc)

infile = raw"D:\temp\saale\control\smf180-lnx-v1.ctl"
outfile = "route.txt"
routeg(infile, outfile)
run(`wsl -e perl -i -pe 's/ß/ss/g;s/[\/]/_/g;s/_/-/g;s/[,,]//g;s/\xc4/Ae/g;s/\xd6/Oe/g;s/\xdc/Ue/g;s/\xe4/ae/g;s/\xf6/oe/g;s/\xfc/ue/g;s/\xdf/ss/g' $outfile`)
run(`cat route.txt`)
x = filter(file -> occursin(r"qgko",file), readdir())
cp(first(x),"qfile")
qdf = qall()
using PrettyTables
sort(qdf,4)|>pretty_table
##grep the specdis_kmu
pt=pwd()*"/route.txt"
CSV.read(pt,DataFrame,header=false,skipto=6,delim="\t",footerskip=6)
#specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/saale_spec_ordered.09"
specfile="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
sim = r"qgko"|>glob|>first|>readdf
obs = readdf(specfile)
simcols = names(sim)
obscols = names(obs)
obscols = map(x->replace(x,r"-" => "_"),obscols)
df = CSV.read(pt,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
rename!(df,1=>"sim",2=>"obs",3=>"name")
#df.name=map(x->replace(x*"-qoutjl",r" # " => ""),df[:,3])  
#df.name=map(x->replace(x,r" # " => ""),df[:,3])
df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
sort!(df, :sim)
@warn "merge with regex of obs!"
map(x->rm(x),glob("qoutjl"))     #rm all 
cnter = 0
for i in eachrow(df)
    cnter += 1
    #cnter = obsc[1]
    obsc = i
    try
        dm =innerjoin(
            #sim[!,Cols(Regex(string(cnter)),end)],
            sim[!,Cols(cnter,end)],
            #obs[!,Cols(obsc[2],end)],    
            obs[!,Cols(Regex(obsc[3],"i"),end)],    
            on=:date)
        onam = obsc[3]*"-qoutjl"
        wawrite(dm,onam)
        println("$onam saved!")
    catch
        onam = obsc[3]*"-qoutjl"
        @warn "merge is empty on $onam ! ..."
    end
end

cntcols("qoutjl")
dpr(r"Bad_Brueckenau")
dpr(r"Wolf")
hcat(qdf,df)|>f -> sort(f,4,rev=true)|>PrettyTables.pretty_table
dpr(r"Br")
facets("tsoil")
rplot(r"qd__")
dfp(r"av")
dfp(r"loc")
dfp(r"sb05")
dfp(r"sb1")


"Wec"|>ftp
"Sac"|>ftp
"Sch"|>ftp

#ggofjl()
gofbatch_nosvg()


# file = "Salz-qoutjl_output.txt"
# output = DelimitedFiles.readdlm(file, '\t', String)
# match = Grep.grep(r"KGE", output)
# sort(match, by = x -> parse(Float64, split(x,dlm="\t")[end]), rev=true)
# x = only(match)
# parse(Float64, split(x)[end])

using DataFrames

function kgedf()
    path = pwd()
    files = glob(r"_output.txt|_outputjl") #non-recurse
    results = DataFrame(File = String[], KGE = Float64[])  # Create an empty DataFrame to store the results
    
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"KGE", output)
        
        if !isempty(match)
            fn = first(split(file, "_qout"))
            
            for line in sort(match, by = x -> parse(Float64, split(x)[end]); rev = true)
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  # remove inner whitespaces
                push!(results, (File = fn, KGE = parse(Float64,split(line)[end])))  # Add the result to the DataFrame
            end
        end
    end
    
    sort!(results, :KGE, rev = true)  # Sort the DataFrame by the 'Line' column in descending order
    return results
end

ds = kgedf()
ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.File)    #mind the DOT before Asterisk
ann = map(x->string.(round(x;sigdigits=1)),ds.KGE)
bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
title = splitpath(pwd())|>last, xrotation = 45, fmt = :png, size = (800, 600), 
fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
annotations = (ds.name,ds.KGE, ann, :top),
bar_width = 0.6)

rsqgrep()
"tjl"|>globdf
dfp(r"glob")
ftp("Schwein")
dpr(r"Schwein")
dpr(r"Wolf")
dpr(r"Mittel")
dpr(r"Bruecken")

r"qba"|>dfp
r"qba"|>bardf
r"qges"|>bardf

dpr(r"Schwein")
nd=readdf(r"qges")
nd=nd[!,Cols(r"19|date")]
dfp!(nd)

ds = kge_df3()
sort(ds, :NSE,rev=true)|>PrettyTables.pretty_table
nsx(ds)

dpr(r"Wolf")
nd=readdf(r"qges")
nd=nd[!,Cols(r"_13|date")]
dfp!(nd)


baryrmean(r"qges")

###plots bars of scores...works only with gr()
@time kgeval()
@time nseval()
@time nsevalraw()
#tdiff
td=tdiff()
wawrite(td,"tdiff-jl.txt")
dfp(td)
baryrsum(td)

import NCDatasets
tdn = tdifnc()
waba()
using Images
img=Images.load("waba-jl.png")

nconly("rad")|>x->x[4]|>rp
nconly("rad")|>first|>rplot

facets("rad")
facets("rad_rcm_1600")

#file remover.
@ncrm
#ts remover.
rmeq()

#include(raw"C:\Users\Public\Documents\Python_Scripts\julia\sd2.jl")
include(raw"C:\Users\Public\Documents\Python_Scripts\julia\sd.jl")
include(raw"C:\Users\Public\Documents\Python_Scripts\julia\frecurs.jl")

npp(infile)






"D:\\Wasim\\regio\\out\\rc200\\r6"|>cd

begin
#function floplt()
    # Search for files matching the specified patterns
    #patterns = ["qbas", "qdi", "qout", "qoutjl"]
    patterns = ["qbas", "qges", "qgko", "qoutjl"]
    #patterns = ["qbas", "qges", "qgko"]
    patterns = map(x->Regex(x),patterns)
    files = filter(file -> any(occursin(r, file) for r in patterns), readdir())
    
    #files = filter(x -> !endswith(x,"txt"),files) #this is in loadalldfs
    
    if isempty(files)
        println("files not present!")
        return 1
    end
    
    # Call the Python script with the appropriate arguments
    #cmd = `python /mnt/c/Users/Public/Documents/Python_Scripts/intmerge_wkrs.py flows.$ext $files`
    #run(cmd)

    dfs = loadalldfs(files)
    # for nd in dfs
    #     nd=nd[!,Cols(r"_13|date")]
    # end

    #dm = mall(dfs)
    dm = mall(dfs)
    dfp(dm)
    dfp(dm[!,Cols(r"_13|date")])
    dfp(dm[!,Cols(r"_23|date")])
    dfp(dm[!,Cols(r"_11|Schwein|date")])
    
    #baryrsum(dm)
    baryrsum(dm[!,Cols(r"_13|date")])
    #baryrsum(dm[!,Cols(r"_26|date")])
end


de = (dm[!,Cols(r"_11|Schwein|date")])
names(de)
de = reorder_df(de)
plotlyjs()
using PlotlyJS
dfplotjs(readf(r"qges"),logy=true,fact=0.70)
dfl(de)

@vv "year"

df = readf(r"qges")
filter!(row -> year(row.date) > 2015, df)
fdf(df)
fdf(r"wind"|>readf)
fdf(r"reg"|>readf)
xdf(r"sb"|>readf)
fdf(r"sb1"|>readdf)

raw"D:\Wasim\regio\out\rc200\v5"|>cd
"D:\\Wasim\\regio\\out\\rc200\\v5\\evap"|>cd
ds = kgedf()
nsegrep()
fdd()
@ncrm
rmeq()
cdb()
#joinpath(raw"D:\Wasim\regio\out","/rc200/v5/station5/")|>cd
raw"D:\Wasim\regio\out/rc200/v5/station5/"|>cd

fdf(r"rad"|>readf)

begin
        patterns = ["qbas", "qges", "qgko", "qoutjl"]
        #patterns = ["qbas", "qges", "qgko"]
        patterns = map(x->Regex(x),patterns)
        files = filter(file -> any(occursin(r, file) for r in patterns), readdir())
        
        #files = filter(x -> !endswith(x,"txt"),files) #this is in loadalldfs
        
        if isempty(files)
            println("files not present!")
            return 1
        end
   dfs = loadalldfs(files)
        # for nd in dfs
        #     nd=nd[!,Cols(r"_13|date")]
        # end
    
        #dm = mall(dfs)
        dm = mall(dfs)
        #dfp(dm)
        dfp(dm[!,Cols(r"_13|date")])
        dfp(dm[!,Cols(r"_23|date")])
        dfp(dm[!,Cols(r"_11|Schwein|date")])
        
        #baryrsum(dm)
        baryrsum(dm[!,Cols(r"_13|date")])
        #baryrsum(dm[!,Cols(r"_26|date")])
    end


##with PlotlyJS()
begin 
    ds = kgedf()
    ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.File)
    ann = map(x->string.(round(x;sigdigits=1)),ds.KGE)
    bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", 
    legend = false, 
    title = splitpath(pwd())|>last, xrotation = 45,
    #fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
    annotations = (ds.name,ds.KGE, ann, :top),
    bar_width = 0.6)
end

    ds = kgedf()
    nrows=size(df)[2]-1
    ti = try
        DataFrames.metadata(df)|>only|>last|>basename
    catch
        @warn "No basename in metadata!"
        ti = raw""
    end
    ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.File)

    
    ds = filter(r -> r.KGE .> 0.1,ds)
    fig = PlotlyJS.make_subplots(
        shared_xaxes=true, 
        shared_yaxes=true    
        );
        PlotlyJS.add_trace!(fig, 
        PlotlyJS.bar(x=ds.name, y=ds[!,:KGE],
        #template="seaborn",
        legend = false, 
        name=ds.name,xrotation = -15))
    fact=0.7
        PlotlyJS.relayout!(fig,
        template="seaborn",
        height=600*fact,width=900*fact,
        title_text=title = splitpath(pwd())|>last)
    display(fig)
    


##with gr()
gr()
begin 
    ds = kgedf()
    ds.name=map(x->replace(x,r"-qoutjl.*" => ""),ds.File)
    ann = map(x->string.(round(x;sigdigits=1)),ds.KGE)
    bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
    title = splitpath(pwd())|>last, xrotation = 45,
    fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
    annotations = (ds.name,ds.KGE, ann, :top),
    bar_width = 0.6)
end


∑(x) = sum(x)
μ(x) = sum(x) / length(x)