using RCall
using DataFrames
using Dates

#julia --startup-file=no --threads auto -q "C:\Users\Public\Documents\Python_Scripts\julia\RCall_gof.jl"
#macro rgof() pt=src_path*"/RCall_gof.jl";include(pt);end
#@rgof

###wa_dd version geht auch...
# function wa_dd(fp::AbstractString, flag::Bool=true)
#     #RCall.@rlibrary data.table

#     R"""
#     library(data.table)
#     x <- fread(
#         file = $fp,
#         header = TRUE,
#         check.names = TRUE,
#         skip = ifelse(test=$flag, yes="YY", no=0),
#         na.strings = "-9999"
#     )
    
#     x = na.omit(x[, lapply(.SD, function(z) as.numeric(as.character(z)))], cols = 1:4)
#     x[, t := do.call(paste, c(.SD, sep = "-")), .SDcols = 1:4]
#     x = x[, t := as.IDate(t)][, -c(1:4)]
#     setkey(x = x, 't')
#     setcolorder(x, 't')
#     """
#     #return rcopy(x)
# end

# function rfread(x::Union{Regex,String})
#     dt = rimport("data.table")

#     if isa(x,Regex)
#         x = dfonly(x)|>first
#     end

#     #rdf = dt.fread(x,nThread=8,skip="Y",colClasses="numeric")        #verbose=true,
#     rdf = dt.fread(x,nThread=8,skip="Y",check=true,strip=true)
#         #na="-9999")
#     dt.setnames(rdf,1,"YY")
#     df = rcopy(rdf)
#     filter!(x -> x.YY != "YY",df)
#     for i in 5:ncol(df)
#         df[!,i] .= replace(df[!,i],-9999.0 => missing)
#     end
#     for i in 1:4
#         if isa(df[!,i],Vector{String})
#             df[!,i] .= tryparse.(Int,df[!,i])
#         end
#     end
#     # dropmissing!(df)
#     dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
#     df.date = dt2
#     df = select(df, Not(1:4))
#     DataFrames.metadata!(df, "filename", x, style=:note)
#     for x in names(df)
#         if startswith(x,"X")
#             newname=replace(x,"X"=>"C", count=1)
#             rename!(df,Dict(x=>newname))
#         end
#     end

#     return df     
# end

@rimport hydroGOF
@rimport data.table as dt

function process_files()
    nms = filter(file -> occursin(r"qoutjl$", file), readdir())
    # Move @rlibrary outside the local scope
    #@rlibrary data.table

    for fcnt in 1:length(nms)
        #x = wa_dd(nms[fcnt])
        #xj = convert(DataFrame, x) #geht
        
        #xj = rfread(nms[fcnt]) #geht
        fn = nms[fcnt]
        xj = convert(DataFrame, dt.fread(fn))
        xj.date = Date.(string.(xj.YY,"-",xj.MM,"-",xj.DD),"yyyy-mm-dd")
        select!(xj,Not(1:4))

        if size(xj)[1] <= 5
            @error("abort: input file has less than 5 lines")
            return
        end
        #Rgof = hydroGOF.gof(sim = xj[:, 2], obs = xj[:, end])
        #println("GOF of Basin $(names(xj)[2]) (sim) to $(names(xj)[end])...")
        Rgof = hydroGOF.gof(sim = xj[:, 1], obs = xj[:, 2])
        println("GOF of Basin $(names(xj)[1]) (sim) to $(names(xj)[2])...")
        printstyled(Rgof, bold = true, color = :blue)
        output_file = nms[fcnt] * "_output.txt"
        open(output_file, "w") do fy    #f fails!
            redirect_stdout(fy) do
                try
                    println("GOF of Basin $(names(xj)[2]) (sim) to $(names(xj)[end])...")
                    println("Performance indices obs|sim:")
                    println(Rgof)
                catch e
                    @error("$e go to next")
                end
            end
        end
    end
    #println("consider ..\n afnse", nms, "hyd \noplat")
    println("..the end")
end

process_files()



# @pyimport pandas as pd
# fn=glob("qoutjl")|>last
# fn=dfonly("sb")|>last
# ydf = pd.read_csv(fn, delim_whitespace=true, na_values=-9999)

# @rimport data.table as dt
# dt.fread(fn)

# df = convert(DataFrame, dt.fread(fn))
# df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd")
# select!(df,Not(1:4))
#DataFrames.metadata!(df, "filename", fn, style=:note)





