module rr

# using RCall
# begin    
#     ## source my R setup now!
#     script_path="D:/Fernerkundungsdaten/Klassifikation/R-Sessions/setup.R"
#     @rput script_path
#     R"source($script_path)"
# end


using RCall

export rfread
export repro_tr
export nctc

function nctc(   
    input_file::AbstractString, 
    output_file::AbstractString, 
    target_srs::Int64;mask = true)
    """
    R terra reprojection
    """
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

function repro_tr(   
    input_file::AbstractString, 
    output_file::AbstractString, 
    target_srs::Int64)
    tr = rimport("terra")
    rin = tr.rast(input_file,drivers="NETCDF")
    """
    reprojection with R's terra package
    """
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
    out = tr.project(n, "EPSG:"*string(target_srs))
    tr.writeCDF(out,output_file,overwrite=true)
    #return out
    println("$output_file done!")
end

function dfonly(x1::AbstractString;recursive=false)
    if recursive==false
        v = filter(file -> occursin(Regex(x1,"i"),file), readdir());
    else
        v = rglob(x1)
    end
    z = v[broadcast(x->!endswith(x,r"nc|png|svg|jpg|txt|log"i),v)];
    return(z)
end

function dfonly(x1::Regex)
    z = filter(file -> occursin(x1,file), 
    readdir()[broadcast(x->!occursin(r"^wq|xml|nc|png|zip|7z|svg|jpg",x),readdir())]);       
    return(z)
end

function rfread_old(x::String)
    dt = rimport("data.table")
    #rdf = dt.fread(x,nThread=8,skip="Y",colClasses="numeric")        #verbose=true,
    rdf = dt.fread(x,nThread=8,skip="Y")
    #dt.setnafill(rdf, fill = 0 , nan=-9999)
    df = rcopy(rdf)
    filter!(x -> x.YY != "YY",df)
    for i in 5:ncol(df)
        df[!,i] .= replace(df[!,i],-9999.0 => missing)
    end
    for i in 1:4
        df[!,i] .= parse.(Int,df[!,i])
    end
    # dropmissing!(df)
    dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
    df.date = dt2
    df = select(df, Not(1:4))
    DataFrames.metadata!(df, "filename", x, style=:note)
    for x in names(df)
        if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
        end
    end
   
    return df     
end

function rfread(x::Union{Regex,String})
    dt = rimport("data.table")

    if isa(x,Regex)
        x = dfonly(x)|>first
    end

    #rdf = dt.fread(x,nThread=8,skip="Y",colClasses="numeric")        #verbose=true,
    rdf = dt.fread(x,nThread=8,skip="Y",check=true,strip=true)
        #na="-9999")
    dt.setnames(rdf,1,"YY")
    df = rcopy(rdf)
    filter!(x -> x.YY != "YY",df)
    for i in 5:ncol(df)
        df[!,i] .= replace(df[!,i],-9999.0 => missing)
    end
    for i in 1:4
        if isa(df[!,i],Vector{String})
            df[!,i] .= tryparse.(Int,df[!,i])
        end
    end
    # dropmissing!(df)
    dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
    df.date = dt2
    df = select(df, Not(1:4))
    DataFrames.metadata!(df, "filename", x, style=:note)
    for x in names(df)
        if startswith(x,"X")
            newname=replace(x,"X"=>"C", count=1)
            rename!(df,Dict(x=>newname))
        end
    end

    return df     
end

function nconly()
    v = readdir();
    z = v[broadcast(x->endswith(x,"nc"),v)];
    return(z)
end

function nconly(x1::AbstractString)
    v::Vector{String} = readdir();
    v = v[broadcast(x->endswith(x,"nc"),v)];
    z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
    return(z)
end

function nconly(x1::Regex)  
    v = filter(x->endswith(x,"nc"),readdir())
    z = filter(file -> occursin(x1,file), 
    v[broadcast(x->!occursin(r"^wq|xml|png|zip|7z|svg|jpg",x),
    v)]);       
    return(z)
end

function rastr(x::Union{Regex,AbstractString})
    tr = rimport("terra")
    if isa(x,Regex)
        x = nconly(x)|>first
    end
    
    tmp = tr.rast(x)
    println("transforming $x ...w masking")
    #tr.NAflag(tmp)=-9999
    #tr.subst(tmp, -9999:-0.1, 0) #replace vals from, to
    tr.subst(tmp, -9999, 0) #replace vals from, to
    #rin[rin <= 0] <- NA #masking in R
    #tr.crs(tmp) = "EPSG:25832" 
    tr.set_crs(tmp, "EPSG:25832" ) #inline
    #rin = tr.flip(tmp)
    rin = tr.flip(
      tr.trans(
      tr.rev(tmp)), 
      direction = "vertical")
      
    return rin
end

# r = rastr(r"^sb1")
# tr.plet(r)
# tr.plot(r,type="classes")

# #RCalls
    # using RCall
    # jmtcars = reval("mtcars");
    # gg = rimport("ggplot2")
    # gg.ggplot(jmtcars,gg.aes(:mpg, :wt))+gg.geom_point()

    #ga = rimport("GGally")
    #ga.ggpairs(jmtcars)
    #ga.ggpairs(df[!,1:2])
    #ga.ggpairs(select(df,Not(:date)))

    # @rlibrary "data.table"
    # # mode: r
    # wa.dd <- function(x,flag=T) {
    #     x <-
    #       data.table::fread(
    #         file = x,
            
    #         header = T,
    #         check.names = T,
    #         skip = ifelse(test=flag, yes="YY",no=0),
    #         na.strings = "-9999"
    #       )
    #     x = na.omit(x[, lapply(.SD, function(z)
    #       as.numeric(as.character(z)))], cols = 1:4)
    #     x[, t := do.call(paste, c(.SD, sep = "-")), .SDcols = 1:4]
    #       x = x[, t := data.table::as.IDate(t)][, -c(1:4)]
    #     data.table::setkey(x = x, 't')
    #     }
    # # time: 2023-03-30 11:32:45 W. Europe Summer Time
    # # mode: r
    #   oo=wa.dd("preci_1970.txt")
    #   md = rcopy(R"oo")
    #   describe(md)

end # module