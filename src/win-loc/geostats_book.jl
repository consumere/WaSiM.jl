
#https://juliaearth.github.io/geospatial-data-science-with-julia/04-geoproc.html   
using GeoStats
import CairoMakie as Mke

N = 100*100
a = [2randn(N÷2) .+ 6; randn(N÷2)]
b = [3randn(N÷2); 2randn(N÷2)]
c = randn(N)
d = c .+ 0.6randn(N)
table = (; a, b, c, d)
gt = georef(table, CartesianGrid(100, 100))

gt|>viewer

p = Point(1, 2)
s = Segment((0, 2), (1, 3))
t = Triangle((0, 0), (1, 0), (1, 1))
b = Ball((2, 2), 1)

geoms = [p, s, t, b]
gset = GeometrySet(geoms)
length(s), area(t), area(b)
[measure(g) for g in gset] # same as #measure.(gset)
#using Unitful is integrated
[1.0, 2.0, 3.0]u"m"
@doc Select # =! @doc select  

box = Box((0, 0, 0), (1, 1, 1))
ball = Ball((0, 0, 2), 0.5)
viz([box, ball], color = ["teal", "slategray3"])

Torus((1, 1, 0), (-1, -1, 0), (1, -1, 0), 0.5)|>viz
Triangle((0, 0), (1, 0), (1, 1))|>Mke.plot
Triangle((0, 2), (1, 0), (1, 1))|>viz
@doc viz

subtypes(Primitive)
grid = CartesianGrid(2, 3)
centroid.(grid)
topo = topology(grid)
A = Adjacency{2}(topo)
A(1)
#We can also query the vertices that are on the 
#boundary of a given quadrangle with the Boundary topological relation. In this case, we specify the parametric dimension of the input and output geometries:
B = Boundary{2,0}(topo)
import Pkg; Pkg.add("GeoIO")
using GeoIO
GeoIO.formats()

fr="D:/Wasim/regio/out/rc200/gout/rad_rcm_1200.mit.nc" #breaks on GeoIO
fn="D:/Wasim/regio/rcm200/v12/catchment.shp"
a = GeoIO.load(fn)
a|>first|>view
a|>first


r=wa.readras(fr)
x = wa.project(fr)
x|>plot


#Interpolation (also used in ClimateTools.jl)
data = georef((z=[1.0, 0.0, 1.0],), [(25, 25), (50, 75), (75, 50)])
grid = CartesianGrid(100, 100)
interp = data |> Interpolate(grid, IDW())
interp|>viewer

seg = Segment((25, 25), (50, 75))
z = interp[seg, "z"]
fig = Mke.Figure()
Mke.Axis(fig[1,1])
for β in [1,2,3,4,5]
  interp = data |> Interpolate(grid, IDW(β))
  Mke.lines!(interp[seg, "z"], label = "β=$β")
end
Mke.axislegend(position = :lb)
Mke.current_figure()
#The larger is the exponent, the more abrupt is the transition of values between the two locations. In addition, the IDW solution will converge to the nearest neighbor solution as β → ∞.

data |> Interpolate(grid, IDW(1, Chebyshev())) |> viewer
#The Chebyshev distance is the maximum of the absolute differences
# between the coordinates of the two points. 
# It is a good choice when the data is sparse and the locations are irregularly distributed.

#Kriging
γ = GaussianVariogram(range=30.0)
data |> Interpolate(grid, Kriging(γ)) |> viewer

fig = Mke.Figure()
Mke.Axis(fig[1,1])
for r in [10,20,30,40,50]
  γ = GaussianVariogram(range=r)
  interp = data |> Interpolate(grid, Kriging(γ))
  Mke.lines!(interp[seg, "z"], label = "range=$r")
end
Mke.axislegend(position = :lb)
Mke.current_figure()


# using GeoArtifacts
# ximg = GeoArtifacts.image("WalkerLakeTruth")

##raster to Geotable +  correlation
ras = Raster("D:/Wasim/regio/rcm200/v12/rcm.slp")
#Rasters.write("D:/Wasim/regio/rcm200/v12/rcm.geojson", img)
#Source driver is raster-only whereas output driver is vector-only
typeof(ximg)

gri = CartesianGrid(size(ras))
val = ras.data
#data = Rasters.points(img)
#table = (; a, b, c, d)
vd = DataFrame(val,:auto)
#table = (; vd.x1, vd.x10, vd.x21, vd.x31)
#gt = georef(table, gri)
Z = ras.data|>vec
table = (; Z)
img = georef(table, gri)
# data = georef((z=[1.0, 0.0, 1.0],), [(25, 25), (50, 75), (75, 50)])
# grid = CartesianGrid(100, 100)
# interp = data |> Interpolate(grid, IDW())
using Random
samples = img |> Sample(10000, replace=false, rng=MersenneTwister(123))
samples |> viewer
#g = EmpiricalVariogram(samples, "x21", maxlag = 100.0)
#g = EmpiricalVariogram(samples, "Z", maxlag = 100.0)
#robust estimator:cressie
g = EmpiricalVariogram(samples, "Z", maxlag = 100.0, estimator = :cressie)
Mke.plot(g)
γ = fit(SphericalVariogram, g)
interp = samples |> InterpolateNeighbors(img.geometry, Kriging(γ))

interp |> viewer
img |> viewer
samples |> viewer

pwd()

#Raster -> Geotable
cd("D:/Wasim/regio/rcm200/v12")
nconly("")
ras = Raster("rcm.intern.sl.nc")
gri = CartesianGrid(size(ras))
Z = (; z=ras.data|>vec)
img = georef(Z, gri)
img |> viewer

img |> Sample(10, replace=false, rng=MersenneTwister(123))
#add units
v=ras.data|>vec
vu = [v]u"°"
slr  = georef((; deg = only(vu)), gri)
slr |> viewer

#with transformation
v = ras[:, :, 1]|>vec #|>permutedims
vu = [v]u"°"

g2 = CartesianGrid(size(ras)[1:2]) #drop last dim i.e. time/layer
slr  = georef((; deg = only(vu)), g2)
slr |> viewer
slr

fn="D:/Wasim/regio/out/rc200/gout/sb1_rcm_1300.mit.nc"
ras = Raster(fn;missingval=0)|>x->x[:, :, 1]

#ras|>plot

dat = ras.data|>vec
slr = georef((; moisture = only(
        [ras.data|>vec|>reverse]u"mm")), 
    CartesianGrid(size(ras)))
slr |> viewer

gt = GeoTable(slr)
#DropMissing(gt, :moisture) 
#DropMissing(gt, 1)

ts = slr |> Filter(row -> !ismissing.(row.moisture))
ts |> Filter(row -> row.moisture != 0.0u"mm") |> viewer
ts |> Filter(row -> row.moisture > 7000u"mm") |> viewer

##Raster -> Geotable this method leads to DimensionMismatch
v = ras.data |> filter(row -> !ismissing.(row))
#v = [ras.data|>skipmissing]
nt = georef((; moisture = only([v]u"mm")), 
    CartesianGrid(size(ras)))
nt |> Filter(row -> row.moisture > 7000u"mm") |> viewer
nt |> viewer


sb = georef((; moisture = only(
        [ras.data|>vec]u"mm")), 
    CartesianGrid(size(ras)))|>
        Filter(row -> !ismissing.(row.moisture) && row.moisture != 0.0u"mm")
sb |> viewer
plot(ras)

#without use of Rasters.jl
using NCDatasets
ds = NCDatasets.Dataset(fn)
dsd = collect(ds.dim |> Dict |> values)
#values_vector = collect(values(ds))
var = Dict(ds)|>keys|>collect|>first
dat = ds[var]|>collect

gri = CartesianGrid(dsd[3], dsd[2])

sb2 = georef((; moisture = only(
        [dat[:, :, 1]|>vec]u"mm")), gri)|>
        Filter(row -> !ismissing.(row.moisture) && row.moisture != 0.0u"mm")
sb2 |> viewer

##w reverse
# gri = CartesianGrid(dsd[3], dsd[2])
# sb2 = georef((; moisture = only(
#         [dat[:, :, 1]|>vec|>reverse]u"mm")), 
#         gri)|>
#         Filter(row -> !ismissing.(row.moisture) && 
#         row.moisture != 0.0u"mm")
# sb2 |> viewer

function rastogeoplt(fn)
    ds = NCDatasets.Dataset(fn)
    dsd = collect(ds.dim |> Dict |> values)
    var = Dict(ds)|>keys|>collect|>first
    dat = ds[var]|>collect
    gri = CartesianGrid(dsd[3], dsd[2])
    sb2 = georef((; moisture = only(
            [dat[:, :, 1]|>vec]u"mm")), gri)|>
            Filter(row -> !ismissing.(row.moisture) && 
            row.moisture != 0.0u"mm")
    sb2 |> viewer
end

fn="D:/Wasim/regio/out/rc200/gout/sb1_rcm_1200.mit.nc"
rastogeoplt(fn)

m = dat[:, :, 1]
typeof(m)

function rastbl(fn)
    if !isfile(fn) 
        error("File not found")
        return
    end
    ds = NCDatasets.Dataset(fn)
    dsd = collect(ds.dim |> Dict |> values)
    var = Dict(ds)|>keys|>collect
    var = filter(z->!occursin(r"^x|^y|^t$",z),var)|>first
    dat = ds[var]|>collect
    gri = CartesianGrid(dsd[3], dsd[2])
    tbl = georef((; Z = dat[:, :, 1]|>vec), gri)|>
            Filter(
                row -> !ismissing.(row.Z) && 
                row.Z != 0)
    return tbl
end

fn="D:/Wasim/regio/out/rc200/gout/qiflrcm_1200.sum.nc"
rplot(fn)

fn="D:/Wasim/regio/out/rc200/gout/vaporcm.mit.nc"
b = rastbl(fn)
b |> viewer
sample = b |> Sample(1000, replace=false)
hscatter(sample, :Z, :Z, lag=0.0)
hscatter(sample, :Z, :Z, lag=5.0)
hscatter(sample, :Z, :Z, lag=15.0)


function readmoist(fn)
    ds = NCDatasets.Dataset(fn)
    dsd = collect(ds.dim |> Dict |> values)
    var = Dict(ds)|>keys|>collect
    var = filter(z->!occursin(r"^x|^y|^t$",z),var)|>first
    dat = ds[var]|>collect
    gri = CartesianGrid(dsd[3], dsd[2])
    sb2 = georef((; Z = only(
            [dat[:, :, 1]|>vec])), gri)|>
            Filter(row -> !ismissing.(row.Z) && 
            row.Z != 0.0)
end
fn = "D:/Wasim/regio/out/rc200/x22/f19/sb1_rcm.2017.nc"
fn = "D:/Wasim/regio/out/rc200/x22/f19/sb05rcm.2017.nc"
img = readmoist(fn)
gₕ = DirectionalVariogram((1.0, 0.0), img, :Z, maxlag = 10.0)
gᵥ = DirectionalVariogram((0.0, 1.0), img, :Z, maxlag = 10.0)
Mke.plot(gₕ, hshow = false, vcolor = "maroon")
Mke.plot!(gᵥ, hshow = false, vcolor = "slategray")
Mke.current_figure()