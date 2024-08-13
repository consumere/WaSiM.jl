##grid adjust
pw()
inpath = "D:/Wasim/regio/rcm200/v4/"
cd(inpath)


wa.ezplot(; dirpath=inpath)

using Rasters
inras = "rcm.kx1"
outras = "rcm.kxn"
newval = 1e-03
#isfile(inras)
if isfile(outras)
    @info "file $outras exists, overwrite ..."
end

begin
    raster = Raster(inras)
    raster.data .= round(newval, digits=4)
    write("tmp.asc", raster; force=true)
    run(`wsl ./rwperl.sh tmp.asc $outras`)
    rm("tmp.asc")
end
cp(outras, "rcm.kyn")


inras = "rcm.aq0"
outras = "rcm.aqx"
newval = 50
#isfile(inras)
if isfile(outras)
    @info "file $outras exists, overwrite ..."
end

begin
    raster = Raster(inras)
    raster.data .= round(newval, digits=4)
    write("tmp.asc", raster; force=true)
    run(`wsl ./rwperl.sh tmp.asc $outras`)
    rm("tmp.asc")
end




heatmap(agread("rcm.ezg");
    title="", xlab="", ylab="")

wa.ezplot()

r = agread("rcm.dhk")
heatmap(r; title="", xlab="", ylab="")
heatmap(r; clims=(111, 980), xlab="", ylab="", axis=(false, false))
surface(r; clims=(300, 980), xlab="", ylab="", axis=(false, true))
#clims=extrema(z|>skipmissing),

ctlg("D:/Wasim/regio/control", "rc200/v4/")
@doc ctlg("D:/Wasim/regio/control", "v4")


# import CairoMakie
#CairoMakie.heatmap(Raster("rcm.dhk"))
#cmk.mkrheat(r;msk=false)
# begin
#     raster = Raster("rcm.kx1")
#     raster.data .= round(1e-03,digits=4)
#     write("rcm.asc", raster;force=true)
#     run(`wsl ./rwperl.sh rcm.asc rcm.kxn`)
# end
# cp("rcm.kxn","rcm.kyn")

#Main.cmk.mkrheat(r"gwst";msk=false)
cd("D:/Wasim/regio/rcm200/v11/")
"D:/Wasim/regio/rcm200/v11/rcm.dhk" |> cmk.mkrheat
fn = "D:/Wasim/regio/rcm200/v11/rcm.s01"
agheat(fn; msk=0.00005, umask=0.0045, roundto=6)
op()
#wit calcs.
w = agread("rcm.wit")
rcont(w)
CairoMakie.heatmap(w)

d = agread("rcm.dep")
CairoMakie.heatmap(d)
vx = d.data |> skipmissing |> collect  #|>vec        #|>dropmissing
vz = w.data |> skipmissing |> collect  #|>vec        #|>dropmissing

CairoMakie.hist(vx)
boxplot(vx, orientation=:horizontal, yaxis=false,
    label="Gerinnetiefe", xmirror=:true)
boxplot!(vz, orientation=:horizontal, yaxis=false,
    label="Gerinnebreite", xaxis=:log10)
#xaxis!(vz)

qplot(vx, vz)



n = @gl ".kol_e04"
cd("D:/Wasim/regio/rcm200/v12/")
cp(n, "D:/Wasim/regio/rcm200/v12/rcm.kol_e04")
n = @gl ".aqx"
cp(n, "D:/Wasim/regio/rcm200/v12/rcm.aqx")
n = @gl ".kyn"
cp(n, "D:/Wasim/regio/rcm200/v12/rcm.kyn")
n = @gl ".kxn"
cp(n, "D:/Wasim/regio/rcm200/v12/rcm.kxn")
ctlg("D:/Wasim/regio/control", "rc200/v12/")


ds = dfr("D:/Wasim/regio/out/rc200/x22/eobs/loc7/dout")
ds = dfr("D:/Wasim/regio/out/rc200/x22/eobs/loc7/pytdiff.txt")
heat(ds)
wa.heat(ds)
wa.heat(ds[!, Not(15:end)])
wa.heat(ds[!, Not(10:15)])
wa.heat(ds[!, Not(7:15)])


#landuse check
lu = Raster("C:\\Users\\chs72fw\\Documents\\EFRE_GIS\\Landcover/20181016_Deutschland_LC_clip_for_Ullmann/v2/LULC_DE_2014_nbg_200m.tif")
import GeoDataFrames as gdf
ez = gdf.read("D:/Wasim/regio/rcm200/ezg.shp")
ez.geometry |> plot
@vv ":transparent"
lx = wa.project(lu; src="EPSG:3034", dst="EPSG:25832")
@doc wa.project(ez.geometry[1])

plot(lx)
plot!(ez.geometry, fillcolor=:transparent; cbar=false)

lc = mask_trim(lx, ez.geometry)

plot(lc)
plot!(ez.geometry, fillcolor=:transparent; title="Landuse")

p = plot(lc; size=(1600, 1200), cbar=true, clims=(1, 10), xlabel="", ylabel="");
plot!(ez.geometry, fillcolor=:transparent; title="Landuse");
savefig(p, raw"D:\Wasim\regio\rcm200\landuse.png")
op()

cd(raw"D:\Wasim\regio\rcm200")
cd(raw"D:\Wasim\regio\rcm200\v12")
rglob("shp")
sa = gdf.read("D:/Wasim/regio/rcm200/v12/catchment.shp")
sa.geometry |> plot

sc = mask_trim(lx, sa.geometry)
p = plot(sc; size=(1600, 1200), cbar=true, clims=(1, 10), xlabel="", ylabel="");
plot!(sa.geometry, fillcolor=:transparent; title="Landuse");
savefig(p, raw"v12-landuse.png")
op()

vgjl("project polygon") ##geht nicht mit ArchGDAL, aber mit ogr2ogr
g = sa
pol = ArchGDAL.IGeometry[]
for i in 1:size(g, 1)
    polygon = g.geometry[i]  # assuming the first geometry is a polygon
    src = EPSG(25832)
    dst = EPSG(4326)
    reprojected_polygon = ArchGDAL.reproject(
        polygon, src, dst) #;order=:compliant
    push!(pol, reprojected_polygon)
end
# gp = g
# gp.geometry .= pol
# plot(gp.geometry, fillcolor=false)
plot(pol)
rll = wa.project(sc; dst="EPSG:4326", src="EPSG:25832")
plot(rll)
plot!(pol, fillcolor=false)

##pwsh cmd:
#ogr2ogr -f "ESRI Shapefile" -s_srs EPSG:25832 -t_srs EPSG:4326 "saale_ll.shp" "catchment.shp"
lshp = gdf.read("D:/Wasim/regio/rcm200/v12/saale_ll.shp")
plot!(lshp.geometry, fillcolor=false)

p = plot(rll; size=(1600, 1200), cbar=false, xlabel="", ylabel="");
plot!(lshp.geometry, fillcolor=:transparent; title="Landuse-v12");
savefig(p, "landuse-v12.png")
oplat()



"D:/Wasim/regio/rcm200/rcm.ezg"
@edit ezplot(; dirpath="D:/Wasim/regio/rcm200")
@edit ezplot()

"D:/Wasim/regio/rcm200/v0/" |> cd
using Rasters
inras = "rcm.kx1"
outras = "rcm.kxn"
newval = 1e-03
#isfile(inras)
if isfile(outras)
    @info "file $outras exists, overwrite ..."
end
begin
    raster = Raster(inras)
    raster.data .= round(newval, digits=4)
    write("tmp.asc", raster; force=true)
    run(`wsl ./rwperl.sh tmp.asc $outras`)
    rm("tmp.asc")
end
cp(outras, "rcm.kyn")



##gcloud
#/home/christian_schaefer_1985/saale/coarse

# #args = "perl -pi -e 's/-9999/10832/g' rcm.art_bfid")
# perl -pe 's/0.0010000000474974513/0.01/g' rcm.ky1 > rcm.ky01
# perl -pe 's/0.0010000000474974513/0.01/g' rcm.kx1 > rcm.kx01
# cp rcm.kx1 rcm.kol
# uni kol


using Rasters

function process_rasters(input_folder::AbstractString, inras::AbstractString, outras::AbstractString, newval::Real)
    cd(input_folder)

    if isfile(outras)
        @info "File $outras exists, overwriting ..."
    end

    begin
        raster = Raster(inras)
        raster.data .= round(newval, digits=4)
        write("tmp.asc", raster; force=true)
        run(`wsl ./rwperl.sh tmp.asc $outras`)
        rm("tmp.asc")
    end

end


# "D:/Wasim/regio/rcm" |> cd
# glob("kx")
# # Example usage:
# input_folder = "D:/Wasim/regio/rcm/"
# inras = "rcm.kx1"
# outras = "rcm.kxn"
# newval = 1e-03
# process_rasters(input_folder, inras, outras, newval)
# cp(outras, "rcm.kyn")


# glob("kol")
# input_folder = "D:/Wasim/regio/rcm/"
# inras = glob("kol")[1]
# rst.stats(inras) #1e06
# outras = "rcm.kol_e04"
# newval = 1e-04
# process_rasters(input_folder, inras, outras, newval)


# input_folder = raw"D:/Wasim/Tanalys/DEM/brend_fab/in3"
# cd(input_folder)
# inras = glob("kol")[1]
# rst.stats(inras) #1e06
# outras = "fab.kol_e04"
# newval = 1e-04
# process_rasters(input_folder, inras, outras, newval)

# inras = "fab.kx1"
# outras = "fab.kxn"
# newval = 1e-03
# process_rasters(input_folder, inras, outras, newval)
# cp(outras, "fab.kyn")
# pw()


# inras = "fab.aq1"
# outras = "fab.aqn"
# newval = 100
# process_rasters(input_folder, inras, outras, newval)

# vgjl("gw_level")
# vgjl("krec")
# dat = """
# gw_level[-]	1	5
# aq1[m]	5	100
# kx		1e-7	1e-4
# ky		1e-7	1e-7
# kol		1e-7	5e-5
# wit[m]	3	10
# dep[m]	3	9
# dr[m-1]	3	50
# krec[-]	0.1	0.9
# kd[h]	10	750
# ki[h]	10	750
# """

# m = readlines(IOBuffer(dat))
# m = [split(x) for x in m]
# df = DataFrame(m, :auto) |> permutedims
# rename!(df, 1 => "name", 2 => "min", 3 => "max")
# #plot(df.name,df.min,df.max,legend=:topleft)
# df.min = parse.(Float64, df.min)
# df.max = parse.(Float64, df.max)
# plot(df.name, df.min, ribbon=(df.min, df.max),
#     fillalpha=0.2, xlabel="", ylabel="Value",
#     label="Min-Max Range", legend=:topleft)
