#read_soildata
@time setup()
lk=raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\tblProfile.csv"
cd(dirname(lk))
#select(df,"GWS")
df = CSV.read(lk,DataFrame,header=true,delim=",")
findindf(df,"Symbol")

sel = findindf(df,"DDc")|>DataFrame
sel = sel[!,Cols(1:5 , 10:ncol(sel))]
view(sel,1:5,1:ncol(sel))

#st = df[df.STATUS=="Leitprofil"]
@vv "ByRo"
DataFrames.subset(df,:STATUS => ByRow(==("Leitprofil")))
filter(x->x.STATUS=="Leitprofil",df)
df[(df.STATUS .==("Leitprofil")), :]
df[(df.STATUS .==("Leitprofil")) .& (df.BOF_NR .>= 1), :]

function dfilter(df::DataFrame, col, val)
    DataFrames.subset(df,Symbol(col)=> ByRow(==(val)))
end

sel = dfilter(df,"STATUS","Leitprofil")
wa.writedf("D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/Leitprofile.tsv",sel)

zp(tovec)

names(sel)|>println


df

df = CSV.read(lk,DataFrame,header=true,delim=",",
    footerskip=0,lazystrings=false,select=1:8)
names(df)|>println
select!(df,1:8)
#println(df)

lk=raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\tblHorizonte.csv"
dh = CSV.read(lk,DataFrame,header=true,delim=",")

lk=raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\tblZuordnungGLE_BLE.csv"
dz = CSV.read(lk,DataFrame,header=true,delim=",")

#lk=raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\tblBlattlegendeneinheit.csv"
# umlaut Blattleg_decoded.csv 

lk="D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/Blattleg_decoded.csv"
dl = CSV.read(lk,DataFrame,header=true)

describe(dl)

# filename = lk
# encoding = "ISO-8859-1"
# # Open the file with the specified encoding
# # Open the file and read its content into an IOBuffer
# file = open(filename, "r")
# buffer = IOBuffer(read(file, String))
# close(file)
# # Convert the buffer's content to a string using the specified encoding
# contents = String(take!(buffer))
# # Now parse the CSV contents using CSV.File with the specified encoding
# data = CSV.File(IOBuffer(contents), header=true, delim=',', String) |> DataFrame


ctlf=raw"D:\Wasim\regio\control\rcm200_x12_loc2.ctl"
dsoil=wa.read_soildata(ctlf)

findindf(dz,"10452")
leg = findindf(dl,"10452")
hor = findindf(dh,"10452")|>DataFrame
pro = findindf(df,"10452")|>DataFrame
cf  = findindf(dsoil,"10452")|>DataFrame
names(cf)|>println
cf[!,Cols("Name")]
hor[!,Cols(r"B")]
pro[!,Cols(r"B")]

#broadcast(x->findmax(x.UTIEF),select(dh,"UTIEF"))

select(dh,"UTIEF")|>sort
select(dh,"UTIEF")|>extrema

#vec(Matrix(select(prec, col)))
function tovec(x::DataFrame, col::Any)
    df = select(x,col)
    println(names(df))
    return vec(Matrix(df))
end
tovec(dh,"UTIEF")|>extrema
tovec(dh,3)|>extrema

#vec(Matrix(select(prec, col)))

findmax(dh.UTIEF)

# 10452 {method = MultipleHorizons;EvapMaxDepth = 0.15;
# Name = Sl2 Sl2 Sl2 Sl3 Sl3 ;horizon = 4 5 6 7 7 ;
# ksat = 3.04444908089364e-06 3.04444908089364e-06 3.04444908089364e-06 2.72947060796587e-06 2.72947060796587e-06 ;
# k_recession = $kr1 $kr2 $kr3 $kr4 $kr5 ;
# theta_res = 0.139387 0.139387 0.139387 0.132053 0.132053 ;
# theta_sat = 0.529878 0.529878 0.529878 0.515036 0.515036 ;
# alpha = 1.3355 1.3355 1.3355 1.1606 1.1606 ;Par_n = 1.224317 1.224317 1.224317 1.245574 1.245574 ;Par_tau = 0.5 0.5 0.5 0.5 0.5 ;
# thickness = 0.15 0.34 0.37 0.24 0.24 ;
# ThicknessScaling = $TS;maxratio = 100 100 100 100 100 ;
# layers = $lyr1 $lyr2 $lyr3 $lyr4 $lyr5 ;}


#############testjoin###############
out=CSV.read("D:/Wasim/main/db_table.csv",DataFrame,header=true)
#,delim=",",    footerskip=1,drop=[1])
bks=CSV.read("D:/Wasim/main/bk_table.csv",DataFrame,header=true)
println(vcat(names(out),"<->",names(bks)))

#bk.join<- merge(x = bkd, y=out, by='TKLE_NR')
dropmissing!(out, :TKLE_NR)       #91213 rows
dropmissing!(bks, :TKLE_NR)       #153 rows
bkj = innerjoin(bks,out,on=:TKLE_NR,makeunique=true)

##WOW##
out_column_names = names(out)
bks_column_names = names(bks)
common_columns = Set(out_column_names) ∩ Set(bks_column_names)

Set(out_column_names) ∩ Set(bks_column_names)        #intersect(s, itrs...)

a,b = map(x->Set(names(x)),[out,bks])
a ∩ b

v = [dl,hor,pro]
a,b,c = map(x->Set(names(x)),v)
a ∩ b ∩ c
b ∩ c
a ∩ b

Set(out_column_names) ∪ Set(bks_column_names)        #union(s, sets...)   
methods(∪)

#crazy:
k = Set(out_column_names)
DataFrame(collect(k)|>permutedims,:auto)
xd = innerjoin(bks,out,on=:TKLE_NR,makeunique=true)
(fnames)|>last                #all funcnames
Set(out_column_names)|>last

# using MacroTools
# macro zp(func)
#     return :(@less $func())
# end
# @zp fdi


@__DIR__ #macro provides the path to the directory containing the currently executing script or module. 

wa.zp("fdi")
wa.zp(dpr)
wa.zp(loadalldfs)



function rbtw(io::IO, start::String, stop::String)
    output = Vector{String}()
    while !eof(io)
        line = readline(io)
        if contains(line, start)
            push!(output, line)
            while !eof(io)
                line = readline(io)
                if contains(line, stop)
                    skip(io, (length(line)))
                    break
                end
                push!(output, line)
            end
            break
        end
    end
    return output
end

pt = "C:\\Users\\Public\\Documents\\Python_Scripts\\julia\\func-win.jl"
readbetween(open(pt),string(fdi), "end")
rbtw(open(pt),string(fdi), "end")



Grep.grep(pt,string(fdi))

output = readlines(pt)
match = 
Grep.grep(string(fdi), output)
regand(string(fdi),"end")


function print_lines_between_patterns(filename::AbstractString, start_pattern::AbstractString, end_pattern::AbstractString)
    in_range = false
    
    for line in eachline(filename)
        if occursin(start_pattern, line)
            in_range = true
        end
        
        if in_range
            println(line)
        end
        
        if occursin(end_pattern, line)
            in_range = false
        end
    end
end

# Usage
filename = ctlf
start_pattern = "theta_res"
end_pattern = "thickness"

print_lines_between_patterns(filename, start_pattern, end_pattern)

print_lines_between_patterns(pt, "fdi", "end")



v = [dl,hor,pro]
a,b,c = map(x->Set(names(x)),v)
a ∩ b ∩ c
b ∩ c
a ∩ b

hp = innerjoin(hor,pro,on=:BF_ID,makeunique=true)
names(hp)|>println
hp=hp[!,Not(Cols(r"Colu"))]
PrettyTables.pretty_table(hp)

names(bkj) ∩ names(hp)
hp.BODSYSTEINH

pwd()
outfile="D:/Wasim/main/saale_jlmerge.tsv"
wa.writedf(outfile,bkj)
zp(wa.writedf)

println(names(bkj))
using RCall
@rimport terra
R"library(terra)"
using GeoDataFrames 
const GDF=GeoDataFrames
GDF.DataFrame(bkj)

gd = GDF.read(outfile)
names(gd)



#outfile|>clipboard
dropmissing!(gd)
#ndf = filter(:BF_ID=> x -> !any(f -> f(x), (ismissing, isnothing)), gd)
gd = filter(:BF_ID=> x -> !any(length(x)==0), gd)
gd.BF_ID= map(x->parse.(Int,x),gd.BF_ID)
#gd[!,"BF_ID"=10659]
select!(gd,Not(Cols(r"Hinweis|geometry|BGL")))
filter(:BF_ID=> x -> any(x==10659), gd)|>  x->x[!,Cols(2,3,7:ncol(x))]|>
PrettyTables.pretty_table 

fl="D:/Wasim/regio/rcm200/v10/nptf.txt"
fl="D:/Wasim/regio/rcm200/v10/dt-ptf.csv"
df = CSV.read(fl,DataFrame,header=true,footerskip=0,lazystrings=false)
df[!,Cols(3:8)]|>x->first(x,12)|>clipboard

using DataFrames

for col in names(df)
    col_data = df[!, col]
    missing_indices = findall(ismissing, col_data)
    
    for idx in missing_indices
        if idx > 1
            col_data[idx] = col_data[idx - 1]
        end
    end
    
    df[!, col] = col_data
end

df

df.ksat isa Number

typeof(df.ksat[5])

col="ksat"
eltype(df[!, col]) == Union{Missing, Float64}

if eltype(df[!, col]) == Union{Missing, Float64}
    col_data = df[!, col]
    missing_indices = findall(ismissing, col_data)
end

using DataFrames

for col in names(df)
    if eltype(df[!, col]) == Union{Missing, Float64}
        col_data = df[!, col]
        missing_indices = findall(ismissing, col_data)
        for idx in missing_indices
            #if idx > 1
                col_data[idx] = col_data[idx - 1]
            #end
        end

        df[!, col] = col_data
        println("nas in $col replaced!")
    end
end

for idx in missing_indices
    println(idx);end

df|>first

findindf(dl,"10452")
findindf(dh,"10452")|>DataFrame
findindf(df,"10452")|>DataFrame
names(dl) ∩ names(sel)
names(dh) ∩ names(sel)

bkj = innerjoin(dh,sel,on=:BF_ID,makeunique=true)

#bkj = innerjoin(bkj,dl,on=:BF_ID,makeunique=true)

names(bkj)|>println
#bkj.SYMBOL

names(bkj) ∩ names(dl) |>println
names(bkj) ∩ "SYMBOL" |>println

v = [dl,sel,dh]
a,b,c = map(x->Set(names(x)),v)

#broadcast(y-> 
hcat(map(x->(names(x)),v)) ∩ "SYMBOL" |>println

#flattened=  join(map(x->(names(x)),v), " ")
flattened = reduce(vcat,(map(x->(names(x)),v), " "))

reduce(vcat,flattened) ∩ "HOR_NR" |>println

f2 = reduce(vcat,flattened)
stringified = join(reduce(hcat,f2), ", ")

join(f2," ") ∩ "TK" |>println
stringified ∩ "TK" |>println
stringified  ∩ "HOR_NR" |>println


k = (map(x->(propertynames(x)),v)) 
flattened = reduce(vcat,k)

join(k," ") ∩ Symbol.("HOR_NR")


"D:/Wasim/Tanalys/DEM/brend_fab/qbfiles/"|>cd
dfs = glob("^qb")
dfs = loadalldfs(dfs)

zp(dfilter)

tm = monmean(dfs[2])
subset(tm,:month=>1)
#tm.nm = getnames(tm)

zp(monsum)

dm = []
for x in dfs
    tm = yrsum(x)
    #tm = DataFrames.subset(tm,:month => ByRow(==(9)))
    tm = DataFrames.subset(tm,:year => ByRow(==(1999)))
    tm[!, "nm"] .= getnames(tm)
    push!(dm,tm)
end

#dm = reduce(hcat, dm)
dm = DataFrame(dm)
sort!(dm,3)
dropmissing!(dm)
DataFrames.subset(dm,:year => ByRow(!=(1999)))
"qbasfab.p1.2016_out-qbas.wa"|>dfp
"qbasfab.v3.1995_out-w3-v4.wa"|>dfp
r"qbasfab.v3.1995_out-qbas.wa_qb"|>dfp


reduced = []
for col in eachcol(dm)
    reduce(vcat, values(col))
    push!(reduced, reduce(vcat, values(col)))
end
dx = DataFrame(reduced,:auto)
dx = DataFrames.subset(dx,:x2 => ByRow(>(10)))
bar(dx[!,6],dx[!,2],legend=false,xrotation = 45)


s = select(DataFrames.subset(dx,:x2 => ByRow(>(50))),:x6)
#s = join(reduce(hcat,s.x6),", ")
collect(s.x6)

typeof(s)
nms = readdir()

ndf = s.x6 ∩ nms
ndf = loadalldfs(ndf)
pall(ndf)

homreg()
nms = rglob("^qbas")
dfs = loadalldfs(nms)

dm = []
for x in dfs
    tm = monmean(x)
    #tm = yrsum(x)
    tm = DataFrames.subset(tm,:month => ByRow(==(9)))
    #tm = DataFrames.subset(tm,:year => ByRow(==(1999)))
    tm[!, "nm"] .= getnames(tm)
    select!(tm,(Cols("month","tot_average","nm")))
    push!(dm,tm)
end

typeof(dm)
#DataFrame(dm|>permutedims,:auto)
dm = reduce(vcat, values(dm))
sort!(dm,2)

findmax(dm.tot_average)

intersect(nms , dm.nm)
k = filter(t->occursin(dm.nm[end],t),nms)
dk = loadalldfs(k)
pall(dk)

".\\x17\\full2\\"|>cd
ds = kge_df3()
nsx(ds)


#qbas für brend
"D:\\Wasim\\Tanalys\\DEM\\brend_fab\\out\\"|>cd
begin
    nms = rglob("^qbas")
    dfs = loadalldfs(nms)
    dm = []
    for x in dfs
        tm = monmean(x)
        #tm = yrsum(x)
        tm = DataFrames.subset(tm,:month => ByRow(==(9)))
        #tm = DataFrames.subset(tm,:year => ByRow(==(1999)))
        tm[!, "nm"] .= getnames(tm)
        select!(tm,(Cols("month","tot_average","nm")))
        push!(dm,tm)
    end
end
typeof(dm)
#DataFrame(dm|>permutedims,:auto)
dm = reduce(vcat, values(dm))
sort!(dm,2)

findmax(dm.tot_average)

intersect(nms , dm.nm)
k = filter(t->occursin(dm.nm[end],t),nms)
dk = loadalldfs(k)
dk = filter(x->nrow(x)>5,dk)
pall(dk)


#qbas für streu
"D:/Wasim/streu/out/"|>cd
nms = rglob("^qbas")
begin
    dfs = loadalldfs(nms)
    dm = []
    for x in dfs
        tm = monmean(x)
        #tm = yrsum(x)
        tm = DataFrames.subset(tm,:month => ByRow(==(9)))
        #tm = DataFrames.subset(tm,:year => ByRow(==(1999)))
        select!(tm,(Cols("month",ncol(tm))))
        rename!(tm,ncol(tm)=>"tot_average")
        tm[!, "nm"] .= getnames(tm)
        push!(dm,tm)
    end
end
typeof(dm)
#DataFrame(dm|>permutedims,:auto)
dm = reduce(vcat, values(dm))
sort!(dm,2)

findmax(dm.tot_average)

#DataFrames.subset(dm,3 => ByRow(>(0.1)))
#filter(t-> t.nm > 0.1,dm)
k = filter(t->occursin(dm.nm[end-1],t),nms)
dk = loadalldfs(k)
dk = filter(x->nrow(x)>5,dk)
pall(dk)

baryrsum(dk[1])

".\\v11\\"|>cd
dfp(r"qout")
dpr(r"qout")
cdb()
ctl()
kgegrepr()

raw"coarse\v1"|>cd
dfp(r"ww-qout")
baryrsum(r"ww-qout")
md = dfr(r"ww-qout")
dropmissing!(md)
ftp(md)

#now saale 
"D:/temp/saale/out_smf200/"|>cd
nms = rglob("^qbas")
nms = (nms[3:end])
nms = nms[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt"),nms)]
nms = nms[broadcast(x->!occursin(r"log|cordex",x),nms)]
println(nms)
dm = []
begin
    for x in nms
        tm = waread(x)
        #tm = wa.xread(x)
        tm = monmean(tm)
        #tm = yrsum(x)
        tm = DataFrames.subset(tm,:month => ByRow(==(9)))
        #tm = DataFrames.subset(tm,:year => ByRow(==(1999)))
        select!(tm,(Cols("month",ncol(tm))))
        rename!(tm,ncol(tm)=>"tot_average")
        tm[!, "nm"] .= getnames(tm)
        push!(dm,tm)
    end
end

typeof(dm)
#DataFrame(dm|>permutedims,:auto)
dm = reduce(vcat, values(dm))
sort!(dm,2)

findmax(dm.tot_average)

dk = DataFrames.subset(dm,2=> ByRow(>(0.01)))
#filter(t-> t.nm > 0.1,dm)
k = filter(t->occursin(dk.nm|>only,t),nms)

# -> qbasrcm.v0.2017
# also really bad...

xd = waread(k|>only)
baryrsum(xd)
dfp(xd)
dpr(r"qout")
ctl()
kgegrepr()

#now goldbach
"D:/Wasim/Goldbach/out"|>cd
nms = rglob("^qbas")
#nms = (nms[3:end])
nms = nms[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt"),nms)]
nms = nms[broadcast(x->!occursin(r"log|cordex",x),nms)]
nms = nms[broadcast(x->!occursin(r"old|yrly",x),nms)]
println(nms)
typeof(nms)|>clipboard

function cntcolread(x::Vector{Any})
    # Get a list of all files in the directory
    files = x

    files = filter(inF->isfile(inF),files)
                        #if isfile(inF)
    file_columns = []
    
    for file in files
        # Open the file for reading
        open(file, "r") do io
            # Read the first line of the file
            line = readline(io)
            # Split the line on tabs
            columns = split(line, '\t')
            # Count the number of columns
            num_columns = length(columns)
            push!(file_columns, (file, num_columns))
        end
    end
    
    # Sort the file_columns array based on the number of columns in descending order
    sorted_files = sort(file_columns, by = x -> x[2], rev = true)
    
    for (file, num_columns) in sorted_files
        printstyled(
            rpad("File: $file",45),
        lpad(" | Columns: $num_columns\n",10),color=:green,bold=true)
    end
    return file_columns
end

cls = cntcolread(nms)
cls = filter(x->x[2]>=3,cls)
#collect(cls)
nmsf = map(x->x[1],cls)

dm = []
begin
    for x in nmsf
        tm = waread(x)
        #tm = wa.xread(x)
        tm = monmean(tm)
        #tm = yrsum(x)
        tm = DataFrames.subset(tm,:month => ByRow(==(9)))
        #tm = DataFrames.subset(tm,:year => ByRow(==(1999)))
        select!(tm,(Cols("month",ncol(tm))))
        rename!(tm,ncol(tm)=>"tot_average")
        tm[!, "nm"] .= getnames(tm)
        push!(dm,tm)
    end
end

typeof(dm)

dm = reduce(vcat, values(dm))
sort!(dm,2)

findmax(dm.tot_average)
dropmissing!(dm)
dk = DataFrames.subset(dm,2=> ByRow(>(0.1)))
k = dk[!,end]
writedf("D:/Wasim/Goldbach/qbas_vals.txt",dk)

#dk = loadalldfs(k)
#pall(dk)

fx = rglob("qbasggb.v21_2011")

dfp(fx[1])
dfp!(fx[2])
dfp!(fx[3])
dfp!(fx[4])

fx = "qbasggb.v46.2021"|>rglob
dfds = loadalldfs(fx)
pall(dfds)
dfp(fx[1])







###########fib bench##########
function fibonacci_julia(n)
    if n <= 0
        return 0
    elseif n == 1
        return 1
    else
        return fibonacci_julia(n-1) + fibonacci_julia(n-2)
    end
end

using RCall
R"
fibonacci_r <- function(n) {
    if (n <= 0) {
        return(0)
    } else if (n == 1) {
        return(1)
    } else {
        return(fibonacci_r(n-1) + fibonacci_r(n-2))
    }
}
"

using PyCall
@pyimport math
function fibonacci_python(n)
    if n <= 0
        return 0
    elseif n == 1
        return 1
    else
        return fibonacci_python(n-1) + fibonacci_python(n-2)
    end
end


# Load the Fibonacci functions

using BenchmarkTools
using BenchmarkPlots

# Set up the benchmarking
n = 30
julia_time = @benchmark fibonacci_julia($n)
r_time = @benchmark R"fibonacci_r($n)"
python_time = @benchmark fibonacci_python($n)

# Compare the benchmarks

plot( julia_time,legend=true,label="Julia")
plot!(r_time,legend=true,label="R",yaxis=:log10)
plot!(python_time,legend=true,label="Python",yaxis=:log10,
title = "Fibonacci Benchmark")


plot( julia_time,legend=true,label="Julia")
#plot!(r_time,legend=true,label="R",yaxis=:log10)
plot!(python_time,legend=true,label="Python",yaxis=:log10,
title = "Fibonacci Benchmark")




D=Dict(
    "Julia" => julia_time,
    "R" => r_time,
    "Python" => python_time
)
values(D)|>first|>plot #or
values(D)|>last|>plot

plot(julia_time,legend=true,label="Julia")
r_time|>plot!
python_time|>plot!
title!("Fibonacci Benchmark")
yaxis!(:log)
#xlabel!(String.(keys(D)))
#xlabel!(join(["Julia","R","Py"]," "))

#xlabel!(String.(keys(D)))
# Get the languages and benchmark times from the dictionary
languages = keys(D)
times = values(D)



###################some selection operations####################
#https://medium.com/@bkamins/dataframes-jl-survey-selecting-columns-of-a-data-frame-based-on-their-values-8c746ce7e389
lk=raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\tblProfile.csv"
cd(dirname(lk))
df = CSV.read(lk,DataFrame,header=true,delim=",")
findindf(df,"Symbol")

select(df, Cols(startswith("a"),
                       any.(ismissing, eachcol(df))))

select(df, Cols(startswith("a"),
                       any.(ismissing, eachcol(df));
                       operator=intersect))

select(df, Cols(startswith("N"),
                       any.(ismissing, eachcol(df));
                       operator=intersect))

select(df, Cols(startswith("B"),
                       any.(!ismissing, eachcol(df));
                       operator=intersect))

pairs(eachcol(df))


filter(row -> any(occursin(Regex("Leitprofil"), 
    string(value)) for value in row), 
    eachrow(df))

wa.findindf(df,"Leitprofil")
wa.findindf(df,r"LeIT")
wa.findindf(df,r"LeIT"i)
wa.findindf(df,r"pseudo"i)
wa.findindf(df,"66")
wa.findindf(df,r"^66$")
#df.STATUS=="Leitprofil"         
#df[df.STATUS=="Leitprofil",:] #wwrong
view(df,1:5,1:ncol(df))
DataFrames.subset(df,:STATUS => ByRow(==("Leitprofil")))
filter(x->x.STATUS==("Leitprofil"),df)

#no --> Select the columns between first and last (including both) from a table.
# filter(x->x.GEN_ID==DataFrames.Between(10,2680),df)
# DataFrames.subset(df,:GEN_ID => ByRow(DataFrames.Between(10,2680)))

select(df, Cols(DataFrames.Between(:GEN_ID,:BODTYP),
                       any.(!ismissing, eachcol(df));
                       operator=intersect))

select(df, Missing .<: eltype.(eachcol(df)))


select(df, Cols(DataFrames.Between(:GEN_ID,:BODTYP),
                Not(Int64 .<: eltype.(eachcol(df))),
                any.(!ismissing, eachcol(df))                
                ;operator=intersect))


select(df, Cols(Float64 .<: eltype.(eachcol(df))))
select(df, Cols(Float16 .<: eltype.(eachcol(df))))
select(df, Cols(Int64 .<: eltype.(eachcol(df))))

typeof.(eachcol(df))


select(df, Cols(Not(String15 .<: eltype.(eachcol(df))),
                any.(!ismissing, eachcol(df))                
                ;operator=intersect))


#occursin("txt"i, names(df),
select(df, Cols(endswith("TXT"),
any.(!ismissing, eachcol(df)); operator=intersect)) 

propertynames(df)|>println
describe(df)
gd = groupby(df, :STATUS);
#gd = groupby(df, :RLFORM);
select(gd, nrow, proprow, groupindices, eachindex)

nd = waread(raw"D:\Wasim\regio\out\rc200\x9\loc7\Bad_Kissingen_Golfplatz-qoutjl")
#hd = groupby(nd, year.(nd.date))
hd = yrsum(nd)
hd = groupby(hd, :year)
select(hd, nrow, proprow, groupindices, eachindex)

dfp(nd)


###wesso###
lk=("C:/Users/chs72fw/Documents/EFRE_GIS/Bodendaten/PTF/wessolek.csv")
wd = CSV.read(lk,DataFrame,header=true,delim=",")
wd|>first

lnk="D:/Wasim/regio/rcm200/v11/rcm.art-bfid"
volcano = readdlm(lnk, ' ';use_mmap=true, skipstart=7,skipblanks = true)
replace!(volcano,"" => NaN64)
v = parse.(Float64,volcano)
#parse.(Int,volcano[1,:]) |>sum


############downscaling results#############
fx = "D:\\remo\\qm\\prec\\wkg.tsw"|>waread
mx = "D:\\remo\\qm\\prec\\wkcdx.tsw"|>waread
simx = "D:\\remo\\qm\\prec\\wkobs.tsw"|>waread
qplot([fx,mx]|>mall)

myr = [fx,mx]|>mall|>yrsum

mse = filter(x->x.year>2010 && x.year < 2050,myr)
mse |> baryrsum
mse |> dfp
mse |> qplot


qplot([fx,simx]|>mall)

[fx,simx]|>mall|>yrsum|>z->filter(u->u.wasserkuppe_gamma>0,z)|>dfp
[fx,simx]|>mall|>z->filter(u->u.wasserkuppe_gamma>0,z)|>baryrsum
[fx,simx]|>mall|>z->filter(u->u.wasserkuppe_gamma>0,z)|>qplot



#nms = DataFrame(CSV.File(pt;header=1,limit=0))
#df = DataFrame(CSV.File(pt;header=false,skipto=2,delim="\t"))
pt = "D:/Wasim/Testdata/Control/wasim_soil_props.txt"
#df = CSV.read(pt, DataFrame;header=false,skipto=2,delim="\t",transpose=true)
df = CSV.read(pt, DataFrame;header=false,skipto=2,delim="\t")
#dropmissing!(df,ncol(df))
dt = permutedims(df)
#nm = dt[2,:]
#nm = readdlm(pt,String;delim='\t')[1:end]
dt = dt[3:end-1,:]
# xm = split(string(nm),"\n")|>last
# xm = split(xm," ")
#|>x->map(y,length(y>0),x)
#filter(s->length(s>0),xm)
new_names = convert(Vector{String}, df[!,2])
rename!(dt,Symbol.(new_names))
cdof(pt)
wa.writedf(dt,"soilprops.txt")
op()

th = dt[28:41,:] #theta sat 
ks = dt[59:end,:] #ksat 
suc = dt[4,:] #suction
hv = dt[42:55,:] #h-values

#innerjoin(hv[!,1],ks[!,1],on=propertynames(ks)[1],makeunique=true)
# sand = hcat(hv[!,1],ks[!,1])
# sand = DataFrame(sand,:auto)
# rename!(sand,Symbol.(["hv","ks"]))
#@df sand plot(:hv,:ks,seriestype=:scatter,title="hv vs ksat",legend=false)
#@df sand corrplot([:hv :ks])

#corrplot(tovec(ks,1),tovec(suc,1))
#tovec(ks,1)|>plot
@df ks plot(cols(propertynames(ks)),xaxis=:identity,
    yticks=tovec(suc,1),thickness_scaling=1.5)
#xaxis!(logscale=true)

@df dt[59:end,:] plot(cols(propertynames(ks)),
    title="k_sat [m/s] for each layer",
    xaxis=:identity,
    #yaxis=:log,
    thickness_scaling=1.1)

"D:/Wasim/regio/rcm200/v11/rcm.art1000"|>agcont2

@rsubset dt  
cb(names(dt))
s = map(x->split(x,"_")|>last,names(dt))
#matches = map(x->match(r"[(].*", x),names(dt))
matches = map(x->match(r"\((.*?)\)", x), s)
inside_brackets = map(m -> m === nothing ? "" : m[1], matches)
string.(inside_brackets)

@df th plot(cols(propertynames(th)),
    title="theta sat for each soiltype",
    xaxis=:identity,
    #yaxis=:log,
    thickness_scaling=1.1)

ksat = dt[3,:]
rename!(ks,Symbol.(inside_brackets))
theme(:tst)
@df ks boxplot(cols(propertynames(ks)),bins=20,
    title="k_sat [m/s] for each layer",
    xaxis=:identity,
    #yaxis=:log,
    thickness_scaling=1.1)

#cl = repeat(1.7909375040117e-05,length(cly))

#clay from ctlfile
cly = Float64.(ks[!,:CL])
value = 1.7909375040117e-05

#sand from ctlfile
cly = Float64.(ks[!,:S])
value = 5.9270139021654e-05
value = 0.00014324652809865

kres = 0.975
vector = Float64[]

push!(vector, value)
for i in 2:21
    value *= kres
    push!(vector, value)
end

qplot(cly,vector)

using GLM
function plot_and_calculate_r2(x::Vector{Float64}, y::Vector{Float64})
    # Fit a linear regression model
    model = lm(@formula(y ~ x), DataFrame(x=x, y=y))
    # Calculate R²
    r2 =  GLM.r2(model)
    # Plot
    p = scatter(x, y, label="data")
    plot!(x, predict(model), label="fit")
    annotate!(p, :bottomright, text("R² = "*string(round(r2, digits=3)), :black))
    return p
end

plot_and_calculate_r2(cly, vector)
plot(cly,label="wasim")
plot!(vector,label="wessolek",thickness_scaling=1.1)

# using Pkg
# Pkg.build("RCall")
vgr("Brook")
using RCall
@rimport LWFBrook90R as lwf

#@rimport data.table as dtt
#ENV|>grep("R_")|>println
#ENV["R_RTOOLS42_PATH"]=nothing
cd("D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8")
ls()

hor=CSV.read("tblHorizonte.csv",DataFrame,
    header=true,delim=",")
hor = dropmissing(hor,:BOART)

bard = unique(hor.BOART)
size(bard)
ts = lwf.hydpar_wessolek_tab(bard)
bard[7]
#"Su4"
lwf.hydpar_wessolek_tab("Su4")
dt
#saturated hydraulic conductivity 
v = 4.44791667663e-05 #ksat in m/s
n = 832.97 #ksat (mm/d)
n_m_s = n * 0.001 / 86400
#n/1000/(24*60*60)
#9.640856503077e-06
#@time setup()
bar([v,n_m_s],
    label=["ctlf","wessolek"],
    title="ksat [m/s] for Su4")

ks
#silty_loam_(SIL)
#stark schluffiger Sand, Su4, 0-8 % T / 40 - 50 % U / 42 - 60 % S
832.97/1000/(24*60*60)
#"fSms"
x = lwf.hydpar_wessolek_tab(bard[2])
y = rcopy(x[6])

#Press $ to activate the R REPL mode and the R prompt will be shown. (Press backspace to leave R REPL mode in case you did not know.)

wa.vgctl("wessolek")
setup()


hometeo()
mkdir("weekly")
cd("weekly")
wa.findctl("f18")
# temperature	# precipitation
# wind_speed    # air_humidity # vapor_pressure
# global_radiation
inpath_meteo = raw"D:\Wasim\Tanalys\DEM\Input_V2\meteo"
pts = ["$inpath_meteo//tmean_ce.txt",
"$inpath_meteo//preci_1970.txt",
"$inpath_meteo//windspeed_ce.txt",
"$inpath_meteo//feuchte_1970.txt",
"$inpath_meteo//dampf_1970.txt",
"$inpath_meteo//radiation_ce24.txt"]

df = dfr(pts[1])
dw = wa.qrtr(df;fun=mean,agg=week)
dfp(dw)


fn=raw"D:\Wasim\Tanalys\DEM\brend_fab\control\fab_c8-loc3.ctl"
dx = wa.read_soildata_4(fn)
#dx = wa.read_soildata_raw(fn)
th = wa.read_soildata_4(fn)|>x->grep(r"thickness",x)
tx = replace.(th," thickness = " => "",";" => "")
tx = map(x->split(x),tx)
tx = map(x->parse.(Float64,x),tx)
thd = DataFrame(reduce(hcat,tx),:auto)
thd = wa.colsums(thd) .* 0.1
nm = wa.read_soildata_4(fn)|>x->grep(r"Name",x)|>x->replace.(x,"Name = " => "",";" => "")|>x->strip.(x)
plot(thd;seriestype=:bar,legend=false,
    xticks=(1:20:length(nm),nm),
    rotation=90)

mx = hcat(nm,thd)
mx = DataFrame(mx,:auto)
sort(mx,2,rev=true)
sort!(mx,2)
rename!(mx,Symbol.(["Feinboden","Horizontmächtigkeit [m]"]))
cdof(fn)

pretty_table(mx;header=(names(mx)), backend = Val(:text))
    
# create a string with the LaTeX code
latex_str = sprint(io -> pretty_table(io, mx, 
    header=uppercasefirst.(names(mx)), backend = Val(:latex)))
dm = (pwd())|>splitpath|>last
write("thickness-table-c8-$dm.tex", latex_str)

open("thickness-c8-$dm.txt", "w") do io
    pretty_table(io, mx; 
    header=uppercasefirst.(names(mx)), 
    backend = Val(:text))
end

"D:/Wasim/Tanalys/DEM/brend_fab/out/c8/loc2/sb05fab.2012.nc"|>rplot
r = Raster("D:/Wasim/Tanalys/DEM/brend_fab/out/c8/loc2/sb05fab.2012.nc")
sl = raw"D:\Wasim\Tanalys\DEM\Input_V2\meteo\preci_1970.txt"
addplot(sl)

rl = wa.project(r;)
plot(rl,xlabel="",ylabel="",title="")
stplot!(sl)



###make new soilmap##########
lk="D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/Leitprofile.tsv"
cd(dirname(lk))
op()
# ]
# generate dbjl
# activate dbjl
# add ODBC,DataFramesMeta
# add DataFrames,CSV

#ODBC doesnt work
begin 
    #"D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/BUEK200DE_Sachdaten_V0.8.accdb"
    using ODBC
    # Define the connection string
    conn_str = """
    Driver={Microsoft Access Driver (*.mdb, *.accdb)};
    Dbq=D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/BUEK200DE_Sachdaten_V0.8.accdb;
    Uid=Admin;
    Pwd=;
    """

    # Connect to the database
    conn = ODBC.Connection(conn_str)
    pt = "D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/BUEK200DE_Sachdaten_V0.8.accdb" ;
    conn = ODBC.Connection(basename(pt))

    # Execute a query
    df = ODBC.query(conn, "SELECT * FROM YourTable")

    # Don't forget to close the connection when you're done
    ODBC.disconnect!(conn)

    # #java conversion
    # mkdir xlsdir
    # sf acc
    # #java -jar xlsdir/client-0.0.5.jar convert --output-format=xlsx BUEK200DE_Sachdaten_V0.8.accdb xlsdir/
    # java -jar xlsdir/client-0.0.5.jar convert --output-format=csv BUEK200DE_Sachdaten_V0.8.accdb xlsdir/
end 

#beziehungen
cd(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8")
@time_imports setup()
file = raw"C:\Users\chs72fw\.julia\logs\repl_history.jl" 
#hd = CSV.read(file,DataFrame,header=false) #very big file...

v=glob("csv")
leg = CSV.read(v[1],DataFrame,header=true,delim=",")
#leg = CSV.read(v[2],DataFrame,header=true)

"tblZuordnungGLE_BLE.csv" :TKLE_NR -> "tblBlattlegendeneinheit.csv"
"tblZuordnungGLE_BLE.csv" :GEN_ID -> "tblProfile.csv" :BF_ID 1 -> ∞ "tblHorizonte.csv"

#tblProfile.csv
#DataFrames.subset(df,:STATUS => ByRow(==("Leitprofil")))
lk="D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/Leitprofile.tsv"
df = CSV.read(lk,DataFrame,header=true,delim="\t")
gle = CSV.read("tblZuordnungGLE_BLE.csv",DataFrame,header=true,delim=",")
hor = CSV.read("tblHorizonte.csv",DataFrame,header=true,delim=",")
names(df)
dm = innerjoin(df,hor,on=:BF_ID,makeunique=true)
dm = DataFrames.subset(dm,:STATUS => ByRow(==("Leitprofil")))
select!(dm,Not(Cols(r"^Column")))

dx = innerjoin(dm,gle,on=:GEN_ID,makeunique=true)
select!(dx,Not(Cols(r"^Column")))

#dz = innerjoin(dm,leg,on=:TKLE_NR,makeunique=true)

select(dx,(Cols(r"^hor|boart|bkid|tief|tk|le_nr|bf_id"i)))
select(dx,(Cols(r"tief"i)))
println(propertynames(dx))

sel = select(dx,(Cols(r"^hor|boart|bkid|tief|tk|le_nr|bf_id"i)))
dropmissing!(sel)
sort!(sel,:LE_NR)
@rsubset leg :1=="231801"
leg.TKLE_NR .= parse.(Int64,leg.TKLE_NR)
dz = leftjoin(sel,leg,on=:TKLE_NR,makeunique=true)
tst = @rsubset dz :TKLE_NR==231801
select(tst,(Cols(r"kurz"i)))

m= map(x->(x=>Int64),[:TKLE_NR,:LE_NR])
m2=map(x->(x=>String),propertynames(leg)[4:5])
hcat(m,m2)
#le2 = CSV.read(v[2],DataFrame,header=true,delim=","
#,    types=m)
describe(le2[:,1])
pwc()

n = CSV.read("ts",DataFrame)


# awk -F, '{print $1}'  $a > ts
# grep -i 'g' ts
# grep -iv 'g' ts
#problem ist, dass in einigen zeilen doch ein string steht
le2 = CSV.read(v[1],DataFrame,header=true,delim=",")
# valst = tryparse.(Int64,le2[!,1])
# leg2 = leftjoin(le2,DataFrame(TKLE_NR=valst),on=:TKLE_NR,makeunique=true)

leg2 = CSV.read(v[1],DataFrame,header=true,delim=",")
leg2.TKLE_NR .= tryparse.(Int64,leg2[!,1])
typeof(leg2.TKLE_NR)
@rsubset leg2 :TKLE_NR==nothing
leg = @rsubset leg2 :TKLE_NR!=nothing
#leg.TKLE_NR .= parse.(Int64,leg.TKLE_NR)

dz = innerjoin(sel,leg,on=:TKLE_NR,makeunique=true)
@rsubset dz :TKLE_NR==231801
dz = groupby(dz,:BF_ID)
#map(x->size(x),dz)
size(dz)
grouped_dz = DataFrames.combine(dz, nrow => :count)
sort!(grouped_dz,:count,rev=true)
pwd()
writedf(grouped_dz,"bfid_counts.txt")
first(dz)
writedf(dz,"bfid_grouped.txt")



r = Raster("D:/Wasim/regio/rcm200/v11/rcm.art-bfid")
v = r.data |> filter(row -> !ismissing.(row))
typeof(v)|>cb
# Assuming v is your Vector{Float32}
v_int = map(x -> convert(Int, x), v)
nt = georef((; bfid = v_int), 
    CartesianGrid(size(r)))
#nt |> Filter(row -> row.moisture > 7000u"mm") |> viewer
nt |> viewer

rd = innerjoin(DataFrame(BF_ID=v_int),sel,on=:BF_ID,makeunique=true)
unique_rd = unique(rd, [1, 2])
@rsubset rd :1>=11202
A = innerjoin(sel,leg,on=:TKLE_NR,makeunique=true)
B = innerjoin(A,unique_rd,on=:BF_ID,makeunique=true)
B = unique(B, [1, 2])
#C = DataFrames.combine(B, (:UTIEF .- :OTIEF) => :thickness)
B.thickness .= B.UTIEF .- B.OTIEF
println(propertynames(B))
select!(B,Not(Cols(r"^Column")))
writedf(B,"Leitprofile_BFID.txt")
pw()


using DataFrames
using Plots

# Group B by BF_ID
grouped_B = groupby(B, :BF_ID)
# Create a scatter plot
p = plot()
for sub_df in grouped_B
    scatter!(p, sub_df.BF_ID, sub_df.thickness, 
        markersize = sub_df.thickness)
end
p

using StatsPlots
s = [:thickness, :OTIEF, :UTIEF]
@df B groupedboxplot(B.BF_ID,cols(s[1]),    #cols(s), 
    legend = :outertopright,
    #xticks = B.BF_ID,
    xrotation = 45,
    xlabel = "", 
    ylabel = "[cm]", 
    title = "")

writedf(B,"soilmap.txt")

@df B boxplot(cols(s[1]), B.BF_ID, 
    orientation = :horizontal,
    legend = false,
    #xticks = :log,
    #yticks = :log,
    xlabel = "", 
    ylabel = "[cm]", 
    title = "")

@df B boxplot(cols(s[1]),B.BF_ID, group = B.thickness,
    orientation = :horizontal,
    legend = false,
    xlabel = "", 
    ylabel = "[cm]", 
    title = "", size = (800*2, 600*2))
savefig("thickness_boxplot.png")
oplat()

valx = B[:,:thickness]
groupx = B.BF_ID
thickness = B.thickness
ax = Mke.boxplot(valx, groupx, color = thickness, 
    horizontal = false)
# Save the figure
Mke.save("boxplot.png", ax)

using RCall
@rimport LWFBrook90R as lwf
v = lwf.hydpar_wessolek_tab("Su4")
#R to julia
M=rcopy(v)
@rimport utils as ut
ut.help("hydpar_wessolek_tab")
ut.help("dataframe")


lk=raw"C:\Users\chs72fw\Documents\EFRE_GIS\Landcover\20181016_Deutschland_LC_clip_for_Ullmann\v2\LULC_DE_2014_nbg_200m.tif"
#r=Raster(lk;crs=EPSG(25832),mappedcrs=EPSG(25832))
r=Raster(lk)
s=gdf.read("D:/Wasim/regio/rcm200/v12/catchment.shp")
r = wa.project(r;src="EPSG:3034",dst="EPSG:25832",imethod=:near)
plot(r)
plot!(s.geometry,fillcolor=false)
rx = mask_trim(r,s.geometry)
plot(rx)


rll = wa.project(rx)
contourf(rll)
wa.rplot(rll)
src="EPSG:25832" 
dst="EPSG:4326"
imethod=:near
flags = Dict(
    "s_srs"=>src,
    "t_srs"=>dst,
    "r"=>imethod)
ts = warp(s.geometry,flags)
#AG method:
rc = wa.reverse_coords(s.geometry)

r

reorder(r, Y)

r = Raster("D:/Wasim/regio/out/rc200/x9/loc7/windrcm.2016.nc")
rev = reverse(r[t=1], dims=(Y,X))
plot(rev)
#surface(r[t=1])
@doc Rasters.reorder




##############dez 23 ######################
struct LandUseData
    id::Int64
    method::String
    rootDistr::Float64
    TReduWet::Float64
    LimitReduWet::Float64
    HReduDry::Float64
    interceptCap::Float64
    julDays::Array{Int64}
    albedo::Array{Float64}
    rsc::Array{Float64}
    rs_interception::Array{Float64}
    rs_evaporation::Array{Float64}
    LAI::Array{Float64}
    Z0::Array{Float64}
    VCF::Array{Float64}
    rootDepth::Array{Float64}
    altDep::Array{Float64}
  end
  
  using CSV

  function read_landuse_data(filename)
    df = CSV.read(filename, DataFrame; header=false,
                  skipto=3,
                  delim=" ", comment="#")
  
    landuse_data = []
    for row in eachrow(df)
      id = parse(Int64, row[1])
      method = row[2]
      rootDistr = parse(Float64, row[3])
      TReduWet = parse(Float64, row[4])
      LimitReduWet = parse(Float64, row[5])
      HReduDry = parse(Float64, row[6])
      interceptCap = parse(Float64, row[7])
      julDays = parse.(Int64, row[8:end])
  
      albedo = zeros(Int64(length(julDays)))
      rsc = zeros(Int64(length(julDays)))
      rs_interception = zeros(Int64(length(julDays)))
      rs_evaporation = zeros(Int64(length(julDays)))
      LAI = zeros(Int64(length(julDays)))
      Z0 = zeros(Int64(length(julDays)))
      VCF = zeros(Int64(length(julDays)))
      rootDepth = zeros(Int64(length(julDays)))
      altDep = zeros(Int64(length(julDays)))
  
      for i in 1:length(julDays)
        albedo[i] = parse(Float64, strip("{}", row[9+(i-1)*12]))
        rsc[i] = parse(Float64, strip("{}", row[10+(i-1)*12]))
        rs_interception[i] = parse(Float64, strip("{}", row[11+(i-1)*12]))
        rs_evaporation[i] = parse(Float64, strip("{}", row[12+(i-1)*12]))
        LAI[i] = parse(Float64, strip("{}", row[13+(i-1)*12]))
        Z0[i] = parse(Float64, strip("{}", row[14+(i-1)*12]))
        VCF[i] = parse(Float64, strip("{}", row[15+(i-1)*12]))
        rootDepth[i] = parse(Float64, strip("{}", row[16+(i-1)*12]))
        altDep[i] = parse(Float64, strip("{}", row[17+(i-1)*12]))
      end
  
      landuse_data = push!(landuse_data, LandUseData(id, method, rootDistr, TReduWet, LimitReduWet, HReduDry, interceptCap, julDays, albedo, rsc, rs_interception, rs_evaporation, LAI, Z0, VCF, rootDepth, altDep))
    end
  
    return landuse_data
  end
  

fn=raw"D:\Wasim\regio\control\lu_table_br-c8.txt"
landuse_data = read_landuse_data(fn)

as = read_soildata(fn)

"""
skips first line after [soil_table] i.e. no of soil types
now returns a DataFrame
"""
function read_soildata(filename::String)
    data = open(filename) do io
        a = readbetween(io, "soil_table", "substance_transport")
        return(a)
    end
    data = broadcast(x -> replace(x,    
    r"^#.*" => "",
    r"^[[].*" => "",
    r"method" => "",
    r"MultipleHorizons" => "",
    r"}" => "",
    r" = " => "=" ), data)


    filter!(s -> !isempty(s), data)
    lines = data[2:end]
    result = []
    # iterate over each line
    for line in lines
        # skip empty lines
        if isempty(line)
            continue
        end
        
        # split the line into number and fields
        number, fields = split(line, " {", limit=2)
        
        # convert the number to an integer
        number = parse(Int, number)
        
        # split the fields into individual fields
        fields = split(fields, ';')
        fields = map(f->strip(f),fields)
        
        # initialize a dictionary to store the data for this line
        data = Dict{String, Any}()
        
        # store the number in the dictionary
        data["number"] = number
        
        # iterate over each field
        try
        for field in fields
            # check if the field contains the " = " substring
            if occursin("=", field)
                # split the field into key and value
                key, value = split(field, "=")
                
                # check if the key is "Name"
                if key == "Name"
                    # keep the value as a string
                else
                  
                    # check if the value is a number
                    if occursin(r"^-?\d+(\.\d+)?$", value)
                        # convert the value to a float
                        value = parse(Float64, value)
                    elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                        # convert the value to a float (scientific notation)
                        value = parse(Float64, value)
                    elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                        # convert the value to a float (scientific notation)
                        value = parse(Float64, value)
                    elseif occursin(r"^\d+ \d+", value)
                        # convert the value to an array of integers
                        value = parse.(Int, split(value))
                        #value = parse.(Float64, split(value))
                    elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)

                        # convert the value to an array of floats
                        value = parse.(Float64, split(value))
                                                
                    end

                end
                
                # store the key-value pair in the dictionary
                data[key] = value
            end
        end
    catch
        @warn "could not parse $value"
        continue
    end
        
        # append the dictionary to the result array
        push!(result, data)
    end
    df = DataFrame(result)
    df.ksat = [parse.(Float64, split(string(x))) for x in df.ksat]
    return df
end  

df = CSV.read(fn, DataFrame; header=false,skipto=3,delim="\t", comment="#")
select!(df,1)
for row in eachrow(df)
    id = parse(Int64, row[1])
end

fn="D:/Wasim/regio/out/rc200/x22/re-0/a.json"
using JSON
m = readlines(fn)
mns = [line for line in m if occursin(r"(?i)timerange|column", line) == false]
JSON.parse()
fn="D:/Wasim/regio/out/rc200/x22/re-0/z.json"
JSON.parsefile(fn)


###jan 29.01.24

soildict = wa.build_soil_dictionary(infile)
soildict[4456]|>permutedims
#soildf=fsoil(infile)
xdf = DataFrame()
# Iterate over the key-value pairs in the dictionary
for (key, df) in soildict
    # Add a column named "key" with the current key value to the DataFrame
    df.key = fill(key, size(df, 1))    
    # Append the current DataFrame to the combined DataFrame
    append!(xdf, df)
end

xdf.ksat = wa.fparse(xdf.ksat)
xdf.ksat[60]|>bar
xdf.ksat[61]|>bar!

dat = select(xdf, [:ksat,:key])
@df dat groupedbar(dat.key,cols(:ksat), legend = :outertopright)
#key to first position
select!(dat, :key, :) #like reorder_df

widen_df = DataFrame([Symbol("ksat_$i") => dat.ksat[i] for i in 1:size(dat, 1)]...)
@doc stack
stack(dat, Not(:key),variable_name=:ksat2)

#dw = transform(dat, :ksat => (x -> DataFrame(x,:auto)) => AsTable)
la=[]
for z in 1:size(dat, 1)
    #dat.$z .= DataFrame(dat.ksat[z],:auto)
    @show DataFrame(dat.ksat[z],:auto)
end
push!(la,DataFrame(dat.ksat[z],:auto))

DataFrame(dat[23,:ksat]|>collect,:auto)


using DataFrames
@rsubset dat begin
    :key .== 4456
end

dat[dat.key .== 4456, :]
firstlayer = [ r.ksat[1] for r in eachrow(dat)]
getproperty(xdf,:key) #is the best way to get the column as vect
#size(df[!,:ksat][1],1)
ok = extract_layers(dat)
#select(xdf, [:Par_m,:key])
nda = select(xdf, [:theta_sat,:key])
nda.theta_sat = wa.fparse(nda.theta_sat)
da = extract_layers(nda,"theta_sat")
#thickness by layers
#dat = select(xdf, [:ksat,:key])
nda = select(xdf, [:thickness,:key])
nda.thickness = wa.fparse(nda.thickness)
thi = extract_layers(nda,"thickness")
npp(infile)
#sums in meter...
thi.sums = sum.(eachrow(thi[!, Not(:key)]))*.1
@view thi[1:3,:]
df=xdf
#map(x->match(r"ksat",x),names(df))
#.captures[1],

Symbol(first(filter(x->occursin(r"ksat",x),names(df))))
propertynames(df)[findfirst(x->occursin(r"ksat",x),names(df))]


"""
Convert DataFrame Column to a Vector
returns only first match, see tovec for multiple matches
or:
m=Symbol.(filter(x->occursin(r"k",x),names(df)))
map(x->getproperty(df, x),m)
DataFrame(map(x->getproperty(df, x),m),:auto)
"""
function vecdf(x::DataFrame, col::Any)
    m = try
            propertynames(df)[findfirst(x->occursin(col,x),names(df))]
        catch
            @error "no match found"
            return
        end
        @info "$m found!"
    return getproperty(x, m)
end

vecdf(df,"sat")
tovec(df,r"sat")


z=dfr(r"win")
wa.reorder_df(z,true)
wa.reorder_df(z)

wa.vecdf

lp = CSV.read("D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/Leitprofile.tsv",DataFrame)
sm = CSV.read("D:/Bodendaten/buek200_2020/BUEK200DE_Sachdaten_V0.8/soildata_wessolek.txt",DataFrame)
println(propertynames(sm))
#soildata
infile = raw"D:\Wasim\regio\control\rcm200_x22-cl5.ctl"
sd = wa.fsoil(infile)
println(propertynames(sd))
thr = wa.extract_layers(sd,"theta_sat")
#select thr for id 4456
@rsubset thr begin
    :key .== 4456
end
@rsubset sm begin
    :BF_ID .== 4456
end

thr[thr.key .== 4456, :]

describe(sm)
import CairoMakie
CairoMakie.lines(sm.thickness)
CairoMakie.boxplot(sm.thickness,1:size(sm,1))
@doc CairoMakie.hist
CairoMakie.hist(sm.thickness;bins=25)

su = rename(sm, :BF_ID => :key)
xm = innerjoin(su,thr,on=:key,makeunique=true)
names(su)|>println
names(xm)|>println

xm.key
@rsubset xm begin
    :key .== 8512
end

#uniqu by row:
un = unique(xm, [:thr, :bart])|>dropmissing
un[!,1:8]|>cb
un[!,Cols(1:8,:BOART_1,r"thick|lay|k_")]|>cb
un|>cb
describe(un)
describe(un[!,Cols(r"alpha|^t|lay|ks")])
describe(sd[!,Cols(r"alpha|^t|lay|ks")])
nv = propertynames(sd[!,Cols(r"alpha|^t|ks")])
lst=[]
for i in nv
    tm=wa.extract_layers(sd,string(i))
    # Create a dictionary that maps the old names to the new names
    rename_dict = Dict(string("layer", j) => string(i, j) for j in 1:5)
    # Rename the columns
    rename!(tm, rename_dict)
    push!(lst,tm)
end
ndf = mall(lst;xcol=:key)
names(ndf)|>println
@df ndf plot(1:size(ndf,1),cols(propertynames(ndf)[1:5]))
@df ndf plot(1:size(ndf,1),cols(propertynames(ndf[!,Cols(r"thi")])))
#@df ndf[ndf.key .== 8512,:] plot(cols(2))
#ndf[ndf.key .== 8512,:]|>permutedims
select!(ndf,:key,:) #key to front...
wa.correlogram(ndf[!,Cols(r"ksat")])
wa.qplot(ndf[!,Cols(r"ksat")],ndf[!,Cols(r"thi")])
writedf(ndf,raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\wes_select.tsv")

bk = tsread(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\wes_select.tsv")
println(names(bk))

a=bk.alpha1
b=bk.alpha2
qplot(a,b)
kge2(a,b)
@pyjl
@doc pyjl.hekge(a,b)
@edit pyjl.hekge(a,b)
"D:/Wasim/regio/out/rc200/x12/spin/"|>cd
@pyimport hydroeval as he

pyjl.kgevec(bk.ksat1,bk.ksat5)
qplot(bk.ksat1,bk.ksat5)
scatter(bk.ksat1,bk.ksat5)
println(names(bk))
wa.correlogram(select(bk,Cols(r"ksat3|theta3")))
wa.correlogram(select(bk,
    Cols(r"^(thet.*!?[2-3])")))
  
select(bk,Cols(r"ksat[1-3]|^(thet.*!?[0-3])"))|>
names|>println


bk = tsread(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\wes_select.tsv")
println(names(bk))
cdof(bk)
glob("Leit")
nk = tsread(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\Leitprofile_BFID.txt")
println(names(nk))
unique(nk.LE_KURZ)|>length
#select(nk, Cols(r"^(BF_ID|LE_KURZ)"))|>cb
unique(nk.BF_ID)

infile=raw"D:\Wasim\regio\control\rcm200_x12_loc10.ctl"
wdf = wa.fsoil(infile)
names(wdf)|>println
wdf.key in nk.BF_ID
#any key in nk.BF_ID?
any(x->x in nk.BF_ID, wdf.key)
#list matching keys.
import DataFrames: subset
#subset(nk,filter(x->x in nk.BF_ID, wdf.key))
#sublist=filter(x->x in nk.BF_ID, wdf.key)
@rsubset nk :BF_ID in wdf.key

unique(nk.LE_TXT)
names(nk)
sel = select(nk, Cols(r"^(BF_ID|LE_TXT|BOART$)"))
#write to tex

latex_str = sprint(io -> 
pretty_table(io, unique(sel,1), header=uppercasefirst.(names(sel)), backend = Val(:latex))
)
write("BF_ID_LE_TXT_BOART.tex", latex_str)
#maby TexTables
#CSV.write("file.tex", unique(sel), format=:tex)

latex_str = sprint(io -> 
    pretty_table(io, sort(nk,:LE_TXT), header=uppercasefirst.(names(nk)), 
    backend = Val(:latex)))
write("Leitprofile_BFID.tex", latex_str)


##todo : Übersicht über die Bodenschichtungen im UG
ug = agread("D:/Wasim/regio/rcm200/v13/rcm.art-bfid")
ug = replace_missing(ug, 0)
parse.(Int64,ug)
#@edit stats(ug)
dct = ug|>skipmissing|>collect|>unique
vgjl("fsoil")
fn="D:/Wasim/regio/control/rcm200_n1-loc1.ctl"
sd = wa.fsoil(fn)
#sd = wa.extract_layers(sd,"theta_sat")	
ugd = filter(x->x.key in dct, sd)       #93soiltypes in UG by art-bfid
tgd = wa.extract_layers(ugd,"thickness")
tgd.sums = sum.(eachrow(tgd[!, Not(:key)]))*.1
@view tgd[1:3,:]
@df tgd scatter(1:size(tgd,1),
    cols(propertynames(tgd)[ncol(tgd)]),
    legend=false)
for i in 1:nrow(tgd)
    val = string(round(tgd.sums[i],
        digits=2))*" m\n"*string(tgd.key[i])
    #*string(tgd.layer1[i])
    
    annotate!(i, tgd.sums[i], 
    Plots.text(val, #"$val m", 
    10, :black, :right,rotation=-35;
    family="Computer Modern"))
end
plot!(size=(800, 600))

nms = select(ugd,Cols(r"Name|key"))
df = innerjoin(nms,tgd,on=:key,makeunique=true)
#select(groupby(df, :Name, sort=true), :Name)
sort!(df,:sums,rev=false)
scatter(df.sums,
    annotations=(1:nrow(df),
    df.sums,
    df.key))
    # Plots.text([first(x) for x in df.Name], 
    # 10, :black, :right,rotation=-35;
    # family="Computer Modern")))
df.l1 = [split(x," ")[1] for x in df.Name]
df.l2 = [split(x," ")[2] for x in df.Name]
df.l3 = [split(x," ")[end-1] for x in df.Name]
df.l4 = [split(x," ")[end] for x in df.Name]
scatter(df.sums, annotations=(1:nrow(df),
    df.sums,
    df.l4),label="last layer")

df

@chain df begin
    @rsubset :sums .< 0.75
    @df bar(:sums, 
#        xlabel = :key,
        xlabel = :Name,
        #xlabel = :l4,
        labelfontsize = 7,
        xlabelrotation = -35,
        xticks = false,
        label="lower 75 cm")
end

dsub = @rsubset df :sums .< 0.95
@df dsub bar(:sums, 
    bottom_margin = 8mm,
    labelfontsize = 7,
    xrotation=-75,
    xticks = (1:nrow(dsub), dsub.Name),
    label="lower 95 cm")

cnms = Symbol.(filter(x->occursin(r"layer",x),names(dsub)))
cnms = cnms|>reverse
lnms = ["layer $i" for i in 1:length(cnms)]
@df dsub bar(
    cols(cnms),
    #yflip = true,
    #yticks = (1:length(cnms), dsub.sums),
    yticks = false,
    yrange = (0, 5),
    bottom_margin = 8mm,
    labelfontsize = 8,
    xrotation=-70,
    xticks = (1:nrow(dsub), dsub.Name),
    label = true)




############## soil 11.03.24 #################
@time setup() #28sec
ug = agread("D:/Wasim/regio/rcm200/v14/rcm.art-bfid")
ug = replace_missing(ug, 0)
#@edit stats(ug)
dct = ug|>skipmissing|>collect|>unique
fn="D:/Wasim/regio/control/rcm200_n2-loc3.ctl"
sd = wa.fsoil(fn)
first(sd.key)
#ardf = DataFrame(key=parse.(Int64, string.(dct)))
ardf = DataFrame(key=Int64.(dct))
sd.key
df = innerjoin(sd,ardf,on=:key,makeunique=true)
#extract_layers(df::DataFrame, colkey::String="ksat")
sel = extract_layers(df,"ksat")
ths = extract_layers(df,"thickness")
#sel = innerjoin(sel,select(df,Cols(r"key|Name|thick")),on=:key,makeunique=true)
sel = innerjoin(sel,select(df,Cols(r"key|Name")),on=:key,makeunique=true)
xm = mall([sel,ths],xcol=:key)
rename!(xm, [Symbol("layer$i") => Symbol("ksat$i") for i in 1:5]...)
rename!(xm, [Symbol("layer$i"*"_1") => Symbol("thickness_$i") for i in 1:5]...)

names(sd)
Tu = select(sd,Cols(r"Name|ksat|theta|alpha|Par_n"))
#list types of cols
propertynames(Tu)|>println
map(typeof,eachcol(Tu))
#ta = [i[2] for i in eachcol(Tu[!,Not(1)])]
nm = [ split(i)[2] for i in (Tu[!,1]) ]
se = [getindex.(Tu[!, col], 2) for col in names(Tu)[Not(1)]]
sed = DataFrame(se,:auto)
sed.Name = nm
#rename!(sed, [i => names(Tu)[Not(1)] for i in 1:4]...)
rename!(sed, 1:5 .=>names(Tu)[Not(1)] )

@rimport LWFBrook90R as lwf
#Alpha parameter of the van Genuchten 
#water retention function (1/m)
@doc lwf.hydpar_wessolek_tab
R"""
require(LWFBrook90R)
help(hydpar_wessolek_tab)
""" 

hp = convert(DataFrame,lwf.hydpar_wessolek_tab(sed.Name))
hp = hcat(sed.Name,hp)
first(hp)
rename!(hp,1=>"Name")
hsub = unique(hp,:Name)
pretty_table
(hsub)|>cb

fi=readlines(raw"D:\Wasim\regio\control\temp")[3:end]
#cat temp | tr '\n' ' ' |sed 's/}/}\n/g' > temp2
write("D:\\Wasim\\regio\\control\\jtemp",join(fi))
#dann 
#sed -i 's/}/}\n/g' jtemp
write("D:\\Wasim\\regio\\control\\A_regex",join(sel.key,"|"))
#dann
#egrep -f A_regex jtemp > new_soil
#wc -l new_soil #186 new_soil 
#perl -lane 'print $F[0]' jtemp 


s=raw"D:\Wasim\regio\control\rcm200_x22-genre-loc.ctl"
@edit wa.fsoil(s)
nsoil= wa.read_soildata(s)
nsoil= wa.fsoil(s)
#wa.read_soildata_4(s)
#DataFrame(nsoil.thickness)
#fparse(nsoil.thickness)
# soiltable = open(s) do io
#     a = readbetween(io, "soil_table", "substance_transport")
#     return(join(a[3:end-1]," ")) #rem first 2 and last lines
# end
# sdict = wa.build_soil_dictionary(soiltable)


sel = extract_layers(nsoil,"ksat")
ths = extract_layers(nsoil,"thickness")
names(nsoil)
boxplot(nsoil.sums)
#bar(nsoil.key,nsoil.sums)
nsoil.sums|>sort|>plot
#findmax(nsoil[!,Not(Cols(3))]|>Matrix|>collect)
sort!(nsoil,:sums)
nd = select(nsoil,Cols(r"key|Name|sums"))
nd.sums = nd.sums*.1
@rsubset nd :sums .< 1
ndm = @rsubset nd :sums .< 1

ann = map(x->Plots.text(
    x,Plots.font(
        family="Computer Modern",:center,8,:black,    
        halign = :left,
    rotation=90.0)),ndm.Name)
@df ndm bar(:sums, 
    bar_width = 0.7,
    #xlabel = :Name,
    #annotations = (1:nrow(ndm), ndm.sums,ann),
    #annotations = (1:nrow(ndm),mean(ndm.sums),ann),
    annotations = (1:nrow(ndm),0.01,ann),
    labelfontsize = 7,
    labelrotation = -35,
    xticks = false,
    label="lower 1 m")

#get the keys.
write("D:\\Wasim\\regio\\control\\A_regex_lower1m",join(ndm.key,"|"))
egrep -f A_regex_lower1m jtemp > new_soil_lower1m
#invert match
egrep -f -v A_regex_lower1m new_soil > new_soil_upper1m


pt="D:/remo/qm/corgrids/rsds/pout"
pt="D:/remo/qm/corgrids/pre/pre_stations2.wa" #<-fehlerhaft. wird gelöscht
pt="D:/remo/qm/corgrids/pre/pre_stations.wa" #
rm(pt)
pt="D:/remo/qm/corgrids/pre/pre_ext.wa"
df = fread(pt)
@chain selt(df,r"Kiss") begin
    @rsubset year.(:date) == 2010 || year.(:date) == 2099
    @rename :Kissingen = $1
    @select(:Kissingen, :date)
    hydro
end

qplot(selt(df,r"Kiss|Kitz"))
#nk = map(byear,df)
heat(df[:,1:10])
heat(df[:,11:12])
#@doc hydro
hydromon(df[:,Cols(1:10,:date)])
#baryrsum(df[:,Cols(1:10,:date)])
@chain selt(df,r"Kiss|Wue") begin
    @rsubset year.(:date) == 2010 || year.(:date) == 2099
    qplot
end

tmp = @rsubset df year.(:date) >= 2020 && year.(:date) <= 2029
baryrsum(tmp[:,Cols(1:5,11,15,:date)])


fn=raw"C:\Users\chs72fw\AppData\Local\Temp\scp06943\home\ext\schaefer\wasim\rcm\control\rcm_adj_v3.ctl"
fn=raw"D:\Wasim\regio\control\new_soil"

df = wa.fsoil(fn)
nk = tsread(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\Leitprofile_BFID.txt")
nksu = filter(x-> x.BF_ID in df.key, nk)
#nk.BF_ID in df.key
println(names(nksu))
rename!(df, :key => :BF_ID)
df = innerjoin(df,nksu,on=:BF_ID,makeunique=true)
nksu = select(df, Cols(r"Name|HORIZ|BOART|BF_ID|LE_KURZ"))
nksu = unique(nksu,:LE_KURZ)
@view nksu[1:12,[1,3,6,5]]
unique(nksu,:Name)
dfu = select(unique(df,:Name),Not(1:2))
select(dfu, Cols(r"Name|key|thickness|ksat"))
getproperty(dfu,:sums) .*.1 |>plot
#size(df[!,:ksat][1],1)
rename!(dfu,:BF_ID => :key)
ok = extract_layers(dfu,"ksat")
ok = stack(ok, Not(:key),variable_name=:lyr)
#groupby=:lyr,
@df ok groupedbar(group = :lyr,cols(:value),
    legend = :outertopright)

sort!(ok,:value)
sort!(ok,:lyr)
scatter(1:nrow(ok),ok.value,group=ok.lyr,yaxis=:log)
plot(ok.value,1:nrow(ok),group=ok.lyr,yaxis=:log,yflip=true)
plot(ok.value,1:nrow(ok),group=ok.lyr,xflip=true,grid=false) #,yflip=true
#again thickenss vs ksat
ths = extract_layers(dfu,"thickness")
rename!(ths, [Symbol("layer$i") => Symbol("thickness_$i") for i in 1:5]...)
ths = stack(ths, Not(:key),variable_name=:lyr, value_name=:thick)
ndf = innerjoin(ok,ths[!,Not(2)],on=:key,makeunique=true)
rename!(ndf,:value => :ksat)
@df ndf scatter(:thick,:ksat,group=:lyr)
@df ndf corrplot([:thick :ksat])

qplot(ndf[!,Cols(r"thick|ksat")])


ks = extract_layers(dfu,"ksat")
dfxm = mall([ths,ks],xcol=:key)
rename!(dfxm, [Symbol("layer$i"*"_1") => Symbol("thickness_$i") for i in 1:5]...)
#correlation between layer and thickness?
#wa.correlogram(dfu[!,Cols(r"thickness")])
@doc wa.corrplot
#theme(:tst)
theme(:dao)
#@df dfu corrplot([:layer1 :thickness_1 :layer5 :thickness_5])
@df dfxm corrplot([:layer1 :thickness_1 ])
ok = stack(dfxm, Not(:key),variable_name=:lyr)
plot(ok.value,1:nrow(ok),group=ok.lyr,yflip=true,grid=false)


# ##using HDF5 ,terra works.
# using Rasters, ArchGDAL, RasterDataSources
# #s = "E:\\glass\\2014\\GLASS01D01.V60.A2014001.h01v08.2022012.hdf"
# s = "E:/glass/lai2014.vrt"
# r = Raster(s)


fn=raw"D:\Wasim\regio\control\rcm200_y22-f18-loc3.ctl"
@time setup() #40sec
df = wa.fsoil(fn)
nk = tsread(raw"D:\Bodendaten\buek200_2020\BUEK200DE_Sachdaten_V0.8\Leitprofile_BFID.txt")
select!(nk,Not(Cols(8:10)))
select(nk,Cols(r"_ID|IEF"))
select!(nk,Not(Cols(r"_1")))
select!(nk,Not(Cols(r"_2")))


z="D:/Wasim/regio/out/rc200/x22/re4/precrcm_1400.sum.nc"
z=readras(z)
plot(z)


