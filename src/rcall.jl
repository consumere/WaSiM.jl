# SPDX-License-Identifier: MIT


module rr
using RCall
try
    import Grep: grep
    rversion = grep("string", convert(Dict, R"version")) |> values |> collect |> first
    printstyled("\n$rversion loaded !\n", color=:green, bold=true, underline=true)
catch e
    @error "R not loaded! \n$e"
    return
end
using DataFrames
using Dates
import DelimitedFiles: readdlm
import ArchGDAL #for wa_terra
using GeoFormatTypes: EPSG, ProjString
tr = rimport("terra") #accessable via rr.tr.
bR = rimport("base") #accessable via rr.bR.
dt = rimport("data.table") #accessable via rr.dt.

export tr      #terra
export rfread  #read time series
export nct     #read netcdf
export tval    #raster values as freq-table
export repro_tr
export reproj

"""
R terra reprojection
"""
function reproj(
    input_file::AbstractString,
    output_file::AbstractString,
    target_srs::Int64; mask=true)
    #@rput mask=T

    if mask
        rmask = "TRUE"
    else
        rmask = "FALSE"
    end

    R"""
    rin = terra::rast($input_file)
    terra::crs(rin) <- 'EPSG:25832'
    terra::NAflag(rin) <- -9999

    rin <- terra::flip(terra::trans(
        terra::rev(rin)), direction = 'vertical')
    #rin[rin <= 0] <- NA #masking
    mask = as.logical($rmask)
    ifelse(mask, yes = {
        rin[rin <= 0] <- NA
    }, no = rin)

    out = terra::project(rin, paste0("EPSG:",$target_srs))
    terra::writeCDF(out,$output_file,overwrite=T)
    """
    return Raster(output_file)
    println("$output_file done!")
end

"""
reprojection with R's terra package
```
repro_tr(
    input_file::AbstractString,
    output_file::AbstractString,
    target_srs::Int64)
```
"""
function repro_tr(
    input_file::AbstractString,
    output_file::AbstractString,
    target_srs::Int64)
    #tr = rimport("terra")
    rin = tr.rast(input_file, drivers="NETCDF")
    # tr.NAflag(rin)=-9999
    # rin = tr.flip(tr.trans(
    #     tr.rev(rin)),
    #     direction = "vertical")

    @rput rin

    R"""
    terra::crs(rin) <- 'EPSG:25832'
    terra::NAflag(rin) <- -9999

    rin <- terra::flip(terra::trans(
    terra::rev(rin)), direction = 'vertical')

    """
    # @rput rin
    # R"""terra::crs(rin) <- 'EPSG:25832'"""

    n = rcopy(R"rin")
    out = tr.project(n, "EPSG:" * string(target_srs))
    tr.writeCDF(out, output_file, overwrite=true)
    #return out
    println("$output_file done!")
end

function dfonly(x1::AbstractString; recursive=false)
    if recursive == false
        v = filter(file -> occursin(Regex(x1, "i"), file), readdir())
    else
        v = rglob(x1)
    end
    z = v[broadcast(x -> !endswith(x, r"nc|png|svg|jpg|txt|log"i), v)]
    return (z)
end

function dfonly(x1::Regex)
    z = filter(file -> occursin(x1, file),
        readdir()[broadcast(x -> !occursin(r"^wq|xml|nc|png|zip|7z|svg|jpg", x), readdir())])
    return (z)
end

function rfread_old(x::String)
    #dt = rimport("data.table")
    #rdf = dt.fread(x,nThread=8,skip="Y",colClasses="numeric")        #verbose=true,
    rdf = dt.fread(x, nThread=8, skip="Y")
    #dt.setnafill(rdf, fill = 0 , nan=-9999)
    df = rcopy(rdf)
    filter!(x -> x.YY != "YY", df)
    for i in 5:ncol(df)
        df[!, i] .= replace(df[!, i], -9999.0 => missing)
    end
    for i in 1:4
        df[!, i] .= parse.(Int, df[!, i])
    end
    # dropmissing!(df)
    dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
    df.date = dt2
    df = select(df, Not(1:4))
    DataFrames.metadata!(df, "filename", x, style=:note)
    for x in names(df)
        if startswith(x, "_")
            newname = replace(x, "_" => "C", count=1)
            rename!(df, Dict(x => newname))
        end
    end

    return df
end

function rfread(x::Union{Regex,String}; rename=true)
    #dt = rimport("data.table")
    if isa(x, Regex)
        x = dfonly(x) |> first
    end
    #rdf = dt.fread(x,nThread=8,skip="Y",colClasses="numeric")        #verbose=true,
    rdf = dt.fread(x, nThread=8, skip="Y", check=true, strip=true) #na="-9999")
    dt.setnames(rdf, 1, "YY")
    df = rcopy(rdf)
    filter!(x -> x.YY != "YY", df)
    for i in 5:ncol(df)
        df[!, i] .= replace(df[!, i], -9999.0 => missing)
    end
    for i in 1:4
        if isa(df[!, i], Vector{String})
            df[!, i] .= tryparse.(Int, df[!, i])
        end
    end
    if rename
        #reads the first line of the file again
        nn = readdlm(x, '\t', String)[1, 5:end]
        for (i, n) in enumerate(nn)
            rename!(df, Dict(i + 4 => n))
        end
    end

    # dropmissing!(df)
    dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
    df.date = dt2
    df = select(df, Not(1:4))
    DataFrames.metadata!(df, "filename", x, style=:note)
    # for x in names(df)
    #     if startswith(x,"X")
    #         newname=replace(x,"X"=>"C", count=1)
    #         rename!(df,Dict(x=>newname))
    #     end
    # end


    return df
end

function rfread_dt(x::Union{Regex,String}; rename=true)
    #dt = rimport("data.table")
    if isa(x, Regex)
        x = dfonly(x) |> first
    end
    #rdf = dt.fread(x,nThread=8,skip="Y",colClasses="numeric")        #verbose=true,
    # x = dt.fread(x,
    #     nThread=8,
    #     skip="Y",
    #     check=true,strip=true,na="-9999")
    # #dt.setnames(rdf,1,"YY")
    #@rput x
    R"""
    require(data.table)
    DT = fread(
        file = $x,
        #nThread = 8, #normal:4 getDTthreads()
        header = T,
        check.names	= T,
        skip = 'YY',
        na.strings = '-9999'
      )
    DT = na.omit(DT[, lapply(.SD, function(z) as.numeric(as.character(z)))], cols = 1:4)
    DT[, t := do.call(paste, c(.SD, sep = '-')), .SDcols = 1:4]
    DT = DT[, date := as.IDate(t)][, -c(1:4)]
    setkey(x = DT, 'date')
    setcolorder(DT, 'date')
    """
    # x = na.omit(x[, lapply(.SD, function(z)
    # as.numeric(as.character(z)))], cols = 1:4)
    #rdf[, date := bR.do_call(bR.paste, c(dt._SD, sep = "-")),SDcols = 1:4]
    #rdf = rdf[, date := dt.as.IDate(t)][, -c(1:4)]
    #data.table::setkey(x = x, 't')
    #setcolorder(x, 't')

    rdf = rcopy(R"DT")
    #filter!(x -> x.YY != "YY",df)
    rdf.date .= DateTime.(rdf.t, dateformat"yyyy-m-d-HH")
    select!(rdf, Not(:t)) #remove t
    return rdf
end

function nconly()
    v = readdir()
    z = v[broadcast(x -> endswith(x, "nc"), v)]
    return (z)
end

function nconly(x1::AbstractString)
    v::Vector{String} = readdir()
    v = v[broadcast(x -> endswith(x, "nc"), v)]
    z = v[(broadcast(x -> occursin(Regex(x1), x), v))]
    return (z)
end

function nconly(x1::Regex)
    v = filter(x -> endswith(x, "nc"), readdir())
    z = filter(file -> occursin(x1, file),
        v[broadcast(x -> !occursin(r"^wq|xml|png|zip|7z|svg|jpg", x),
            v)])
    return (z)
end

function nct(x::Union{Regex,AbstractString})
    #tr = rimport("terra")
    if isa(x, Regex)
        x = nconly(x) |> first
    end

    tmp = tr.rast(x)
    println("transforming $x ...w masking")
    #tr.NAflag(tmp)=-9999
    #tr.subst(tmp, -9999:-0.1, 0) #replace vals from, to
    tr.subst(tmp, -9999, 0) #replace vals from, to
    #rin[rin <= 0] <- NA #masking in R
    #tr.crs(tmp) = "EPSG:25832"
    tr.set_crs(tmp, "EPSG:25832") #inline
    #rin = tr.flip(tmp)
    rin = tr.flip(
        tr.trans(
            tr.rev(tmp)),
        direction="vertical")

    return rin
end

# r = nct(r"^sb1")
# tr.plet(r)
# tr.plot(r,type="classes")

"""
deprecated see rgof()
"""
function gof6()
    sc = "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof6.R"
    R"""source($sc)"""
end #deprecated

#moved to 4.3.2 -> add old path to .libPaths()
lpt = raw"C:/Users/chs72fw/AppData/Local/R/win-library/4.2"
R""".libPaths(new=$lpt)"""

"""
´´´
rgof(fn::Union{String,Regex,DataFrame};returnDF=false)
´´´
"""
function rgof(fn::Union{String,Regex,DataFrame}; returnDF=false)
    if isa(fn, DataFrame)
        x = fn
    else
        x = rfread(fn)
    end
    @rput x
    R"""
    x = xts::as.xts(x)
    hydroGOF::ggof(
        sim = x[, 1], obs = x[, NCOL(x)],
        legend = names(x),
        FUN=mean,
        col = c('#d3bfbf', '#6464cc'),
        gofs = c('MAE', 'PBIAS', 'R2', 'KGE', 'VE'),
        na.rm = T,
        ylab = '[mm]',
        leg.cex = 1.2,
        pch = 21,
        ftype='ma')
    """
    if returnDF
        return x
    end

end

"""
´´´
tval()
like Rs table(terra::values(r,mat=T,na.rm=T))
´´´
"""
function tval(r::RObject{S4Sxp})
    #import base R
    #@rimport "base" as bR #see above at module
    tab = bR.table(tr.values(r, mat=true, na_rm=true))
    vals = parse.(Float64, convert(Array, bR.names(tab)))
    freq = convert(Array, tab) #parse.(Int,convert(Array,tab))
    return DataFrame(value=vals, frequency=freq)
end

"""
´´´
wa_terra(file::String, turn::Bool = false, crsval::Int = 25832, encoding::String = "UTF-8")
usage:
s="D:/Wasim/Tanalys/DEM/Input_V2/meteo/specdis_kmu.txt"
st=wa_terra(s,turn=false,crsval=25832,encoding="UTF-8") #-> SpatVector
a=convert(Matrix,tr.as_data_frame(st)) #-> DataFrame of Names and EZG
b=convert(Matrix,tr.crds(st)) #-> Matrix of crds
hcat(a,DataFrame(x=b[:,1],y=b[:,2])) #-> DataFrame of Names, EZG, x, y
´´´
"""
function wa_terra(file::String, turn::Bool=false, crsval::Int=25832, encoding::String="UTF-8")
    # Read the data file using R's fread function from data.table
    if !isempty(file)
        @rput file encoding
        R"st <- data.table::fread(input = file, skip = 0, header = TRUE, nrows = 4, encoding = encoding)"
        @rget st
    else
        println("Error: please provide a valid file path to read in.")
        return nothing
    end

    # Define how to handle the data based on the 'turn' parameter
    if turn
        # Extract longitude and latitude when `turn` is true
        @rput st
        R"""
        cords <- data.frame(
            lon = as.numeric(t(st[3, 5:ncol(st)])),
            lat = as.numeric(t(st[2, 5:ncol(st)]))
        )
        ez <- st[1, 5:ncol(st)]
        st <- terra::vect(cords, crs = paste0("epsg:", $crsval))
        st$rn <- names(ez)
        st$ezg <- as.numeric(ez)
        """
        @rget st
    else
        # Extract longitude and latitude when `turn` is false
        @rput st
        R"""
        cords <- data.frame(
            lon = as.numeric(t(st[2, 5:ncol(st)])),
            lat = as.numeric(t(st[3, 5:ncol(st)]))
        )
        ez <- st[1, 5:ncol(st)]
        st <- terra::vect(cords, crs = paste0("epsg:", $crsval))
        st$rn <- names(ez)
        st$ezg <- as.numeric(ez)
        """
        @rget st
    end

    #return st
    a = convert(Matrix, tr.as_data_frame(st)) #-> DataFrame of Names and EZG
    b = convert(Matrix, tr.crds(st)) #-> Matrix of crds
    aa = hcat(a, DataFrame(x=b[:, 1], y=b[:, 2])) #-> DataFrame of Names, EZG, x, y
    pts = ArchGDAL.IGeometry[]
    for i in axes(aa)[1]
        pt = ArchGDAL.createpoint([aa.x[i], aa.y[i]])
        pt = ArchGDAL.reproject(pt, EPSG(25832), ProjString("+proj=longlat +datum=WGS84 +no_defs"))
        push!(pts, pt)
    end
    od = DataFrame(geometry=pts, name=aa.rn, ezg=aa.ezg, xc=aa.x, yc=aa.y)
    return od
end

"""
this sources from /rfile/setup.R
not /mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions/setup.R
"""
function rsetup()
    #script_path="D:/Fernerkundungsdaten/Klassifikation/R-Sessions/setup.R"
    script_path = joinpath(dirname(src_path), "rfile", "setup.R")
    @rput script_path
    R"source($script_path)"
    # #moved to 4.4.1 -> add old path to .libPaths() <-already in module
    # lpt = raw"C:/Users/chs72fw/AppData/Local/R/win-library/4.2"
    # R""".libPaths(new=$lpt)"""
end

end # endof module

# using RCall
# function rrtoMain()
#     fnames = names(Main.rr, all=true)
#     for submodule in fnames
#         @eval import Main.rr.$submodule
#     end
# end
println("rr Module loaded!")
#println("use rrtoMain for loading all submodules")
#println("using .rr for loading all submodules to Main...")

#v = dfonly(r"qoutjl$")
# @time dfs2 = map(dfr,v)
# @time dfs = map(rr.rfread,v)
# @time dfs3 = map(waread2,v) #fastest

# x = (dfs[2])
# @rput x
# R"""ggpairs(x)"""
# rr.rgof("Wolfsmuenster-qoutjl")
# rr.rgof(x)

# x = (dfs[3])
# @rput x
# R"""str(x)"""
# R"""
# setDT(x)
# x[,yr:=year(x$date)]
# ggally_facetdensitystrip(x, aes(y = x$C320, x = x$yr))
# """

# R"version" #4.2.2

# l="D:/Wasim/regio/out/rc200/z4/windrcm.2017.nc"
