###global df plots

plotlyjs()
"/mnt/d/Wasim/Tanalys/DEM/brend_fab/out/m3"|>cd

dfs=dfs=loadalldfs(p())
filterplot("qg",dfs)
filterplot!("wind",dfs)
aplot(dfs[3])

lplot("qg",dfs)
lplot("vap",dfs)
denseplot("vap",dfs)
denselog("qges",dfs)
vio("vap",dfs)

dfpall(dfs)

#or from current dir 
("tst2")|>lplot
md = getf("gw")
getdf("gwn",md) |>dfp


md = getdf("vap",dfs)

lplot(md)

denselog("sb",dfs)
getdf("sb",dfs)|>lplot
filterplot!("temper",dfs)

getdf("sb",dfs)|>dfp

##even in pipe...
getf("sb")|>last|>dfplot

#with inex of :) 
getindex(getf("sb05"), 2)|>dfp

v=getf("te")

push!(v,0)  #add 0 to vector
v[1:end-1] # Wählen Sie alle Elemente außer dem letzten aus
pop!(v)    #Entfernen Sie das letzte Element aus v und geben Sie es zurück
m=pop!(v)
push!(v,m)  #add m to vector

length(v)

###and yes, this works with rasters as well...
getf("sb05")[3]|>rplot
v = regand(getf("sb"),"sb1","nc")
#contourplot in js dauern..
gr()
getindex(v,2)|>rp


a=readdf(first(getf("qges")))
b=readdf(first(getf("vapo")))

v = filter(file -> !endswith(file, "nc"), readdir())
#v = regand(v,"sb1","q")
vv = []
# for x in v
#       if (filter(file -> occursin("q",file),x))&&(filter(file -> startswith(file,"so"),x))
#             push!(vv,x)   
# end
# end   
a=(filter(file -> occursin("q",file),v))
b=(filter(file -> startswith(file,"so"),v))
#join(a,b)

filter(file -> startswith(file,"so"),v)
occursin(a,b)

filter(n->occursin(Regex(regex,"i"),n),
    map(x->basename(only(DataFrames.metadata(x))[2]),dfs)
    )

ad = loadalldfs(a)
xx=getdf(ad,2)

function dfpall(dfs::Vector{DataFrame})
      p1 = dfp(first(dfs))
      res = []
      push!(res,p1)
      sel = dfs[2:length(dfs)]
      for xx in sel
          println("adding ",basename(only(DataFrames.metadata(xx))[2]))
          y = filter(x->!occursin("date",x),
          names(xx))
          s = map(y -> Symbol(y),y)
          #ti = basename(only(DataFrames.metadata(xx))[2]);
          px=@df xx Plots.plot(:date,cols(s),
          #yaxis = :log,
          #title = ti,
          legend = :topright)
          push!(res,px)
          return(res)
      end
  end
#  DataFrame(x = 1:10, y = rand(10), group = "A")
#df = vcat(ad[2], ad[end])
df = vcat(ad[2], ad[1])

df = innerjoin(ad[7], ad[end],on = :date,makeunique=true)

dfs=ad
df = reduce((left, right) -> 
      innerjoin(left, right, on = :date,makeunique=true), 
      dfs)
y = filter(x->!occursin("date",x), names(df))
s = map(y -> Symbol(y),y)
@df df Plots.plot(:date,cols(s),#yaxis = :log,
         legend = :topright)


plotlyjs()
a=getf("so_t")
#ad = loadalldfs(a)          
ad = loadalldfs(r"q")
dfpall(ad)
dfpall(ad[2:5])



v[1]
v[2]

#########maths.....
A = [1 2 3]'
B = [1 2 3]

using LinearAlgebra
# Berechnen Sie das Skalarprodukt
dot(A, B)
A ⋅ B  # Eine andere Möglichkeit ist die Verwendung des Symbols ⋅, 
      # das durch die Eingabe von \cdot gefolgt von der TAB-Taste 
      #erzeugt werden kann. Dieses Symbol ist ein Operator, der das Skalarprodukt zwischen seinen Operanden berechnet.


#A .* B




dfs=dfs=loadalldfs(p())
filterplot("gw",dfs)
filterplot!("wind",dfs)
aplot(dfs[3])

dfp("gwn",dfs)

dfs=[]
files = filter(x -> !occursin(r"grd|xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg|yr|log",x),v) 
for file in files
	if isfile(file) && (!endswith(file,"png|grd|xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|png|svg"))
	           file_path = file
	       println("reading ",file_path,"...")
	       p1 = readdf(file_path)
	       push!(dfs, p1)
	        end
	    end

out = reduce((left, right) -> 
          innerjoin(left, right, on = :date,makeunique=true), 
          ddf)

"/mnt/d/Wasim/regio/out/c8"|>cd

dfs = loadalldfs(pwd())
dfp("qout",dfs)

########barplots##############
#using StatsPlots, DataFrames
#yearsums = DataFrame(year=2010:2020, sales=rand(11), profit=rand(11))

@df yearsums plot(:year, [:sales :profit], seriestype=:bar)

df = readdf(r"qges")
df|>first
y = filter(x->!occursin("date",x),names(df))
s = map(y -> Symbol(y),y)
o = DataFrames.metadata(df)|>collect
ti = basename(o[1][2])
@df df Plots.plot(:date,cols(s),legend = :topright, title=ti,seriestype=:bar)

yearsums = toyrsum(df)
yearsums|>first

#ad = cattoyrsum(df)
#df_yearsum = by(df, :year, names(df)[2:end] .=> sum .=> names(df)[2:end])
#
df = readdf(r"qges")
y = filter(x->!occursin("date",x),names(df))
s = map(y -> Symbol(y),y)
df[!, :year] = year.(df[!,:date]);
df_yearsum = combine(groupby(df, :year), y .=> sum .=> y) 
@df df_yearsum Plots.plot(:year,cols(s),
legend = :topright, title=ti,seriestype=:bar)

bardf(r"wind")
ve=rglob("wind")
ve[3]|>typeof
bardfm(r"wind")
bardfm(ve[3])

yrsum(ve|>first)
readdf(ve|>last)|>yrsum
cnt()

readdir(abspath("."), join=true)