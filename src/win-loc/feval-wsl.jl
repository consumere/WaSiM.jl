##No julia gofplot
# no ncrm. run on ubu
#pt=/mnt/c/Users/Public/Documents/Python_Scripts/julia/win/feval-wsl.jl 
#julia --startup-file=no $pt 

#raw"D:\Wasim\regio\out\rc200\x22\reloc2"|>cd
if Sys.isapple()
    platform = "osx"
    const homejl = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
    const mybash = "/Users/apfel/.bash_aliases"
    src_path = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
elseif Sys.iswindows()
    platform = "windows"
    src_path = "C:\\Users\\Public\\Documents\\Python_Scripts\\julia"
    macro wasim() pt="C:\\Users\\chs72fw\\.julia\\dev\\WaSiM\\src\\wa.jl";include(pt);end
else
    platform = "unix"
    winpt = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
    pcld = "/home/cris/pCloud Drive/Stuff/Python_Scripts/julia"
    src_path = isdir(winpt) ? winpt : pcld
    println("sourcepath is $src_path")
    if isdir(winpt)
        macro wasim() pt="/mnt/c/Users/chs72fw/.julia/dev/WaSiM/src/wa.jl";include(pt);end
    end
end 

println("eval on: ")
printstyled(pwd()*"\n",color=:green)

@time using DataFrames,CSV,Dates,Statistics
printstyled("loading StatsPlots\n",color=:blue)
@time using StatsPlots, Plots
using Plots.PlotMeasures
import Grep: grep
#import DelimitedFiles: readdlm


#@time import PlotThemes
Plots.theme(:dao)

##########funcs##########
begin
    function mvwasim2(;ta=pwd(),pt="C:/Users/chs72fw/Documents/EFRE_GIS/Hydrologie/WaSiM-ETH/wasimvzo64_10.06.05",kw...) 
            
        println("\nmoves all wq, xml and log files to from
        $pt  to current pwd")
        
        println("target dir is $ta");

        af = filter(x -> occursin(r"wq", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i));kw...) #force=false
            println(basename(i)," --> ", ta)
        end

        af = filter(x -> occursin(r"xml", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i));kw...)
            println(basename(i)," --> ", ta)
        end
        af = filter(x -> occursin("modell", x), readdir(pt,join=true))
        for i in af
            mv(i,joinpath(ta,basename(i));kw...)
            println(basename(i)," --> ", ta)
        end

    end

    function kge1(simulations, evaluation)
        r = cor(simulations, evaluation)
        α = std(simulations) / std(evaluation)
        β = mean(simulations) / mean(evaluation)
        return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
    end

    """
    sim, obs
    hydroeval.objective_functions.nse(simulations, evaluation)
    """
    function nse(predictions::Vector{Float64}, targets::Vector{Float64})
        return (1 - (sum((predictions .- targets).^2) / sum((targets .- mean(targets)).^2)))
    end

    """
    simulated, observed
    """
    function vef(sim, obs)
        return (1 - ( sum( map(x->abs(x),(obs - sim) ) ) / sum( obs ) ))
    end

    function ftplin(df::DataFrame)
        if size(df)[2]!=3
            #throw(@warn "wrong number of columns - using dfp!")
            @warn "wrong number of columns - using dfp!"
            @warn "need :Simulated,:Observed,:Date !"
            display(dfp(df))
            return
        end
        rename!(df,map(x->replace(x,r"_"=>" "),names(df)))
        ndf = copy(df)
        ndf = hcat(ndf[!,Not(Cols(r"date"i))],ndf[:,Cols(r"date"i)])
        #ndf = filter(:date=> x -> !any(f -> f(x), (ismissing, isnothing)), ndf)
        rename!(ndf, [:Simulated,:Observed,:Date])
        dropmissing!(ndf) ##hmm sketchy..
        overall_pearson_r = cor(ndf[!, :Simulated],ndf[!, :Observed])
        r2 = overall_pearson_r^2
        #nse(simulations, evaluation)
        nse_score = nse(ndf[!, :Simulated],ndf[!, :Observed])
        kge_score = kge1(ndf[!, :Simulated],ndf[!, :Observed])
        ve = round(vef(ndf[!, :Simulated],ndf[!, :Observed]), digits=2)
        subs = "RSQ: $(round(r2, digits=2))\nNSE: $(round(nse_score, digits=2))\nKGE: $(round(kge_score, digits=2))\nVE: $ve"
        ti = try
            split(basename(last(collect(DataFrames.metadata(ndf)))[2]),"-")[1]
        catch
        @warn "No basename in metadata!"
            raw""
        end 
        
        p = plot(
            ndf[!, :Date], ndf[!, :Simulated], 
        title=ti, 
        line=:dash, 
        color=:blue, 
        label=names(df)[1],
        ylabel="[mm/Tag]", 
        xlabel="", 
        size=(1200,800) .* 0.7,
        right_margin = 10mm,
        left_margin = 5mm,
        top_margin = 5mm,
        bottom_margin = 15mm,         
        legend=:topleft)
        plot!(p, ndf[!, :Date], 
            ndf[!, :Observed], 
            color=:red, 
            label=names(df)[ncol(df)])
        annotate!(
        :topright,
        Plots.text("$subs", 12, :black, 
        :right;
        family="Computer Modern"))
        return p
    end

    """
    newer version with copy df and switched func positions
    """
    function wawrite(df::DataFrame,file::AbstractString)
        dout = copy(df)
        #dout[!,Cols(r"date")]
        #in("date",names(dout))
        if in("year",names(dout))
            @warn "yearcol found!"
            CSV.write(file, dout, 
            transform = (col, val) -> something(val, missing), delim="\t")  
            return
        end
        dout.YY = map(x ->year(x),dout.date)
        dout.MM = map(x ->month(x),dout.date)
        dout.DD = map(x ->day(x),dout.date)
        dout[!, "HH"] .= 0
        dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
        CSV.write(file, dout, 
        transform = (col, val) -> something(val, missing), delim="\t")  
        nothing
    end

    """
    dfr = waread
    """
    function waread(x::Regex)
        """
        Read the text file, preserve line 1 as header column
        """
        x = glob(x)|>first
        ms = ["-9999","lin","log","--"]
        df = CSV.read(x, DataFrame; delim="\t", 
            header=1, missingstring=ms, 
            silencewarnings=true,
            normalizenames=true, types=Float64)
        df = dropmissing(df, 1)
        dt2 = map(row -> Date(Int(row[1]), Int(row[2]), Int(row[3])), eachrow(df))
        df.date = dt2
        df = select(df, Not(1:4))
        metadata!(df, "filename", x, style=:note)
        #renamer
        for x in names(df)
            if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    """
    readdf for loop 
    --- main reader ---
    delim: if no argument is provided, 
    parsing will try to detect the most consistent delimiter on the 
        first 10 rows of the file
    """
    function readdf(x::AbstractString)
        ms=["-9999","lin","log"]
        df::DataFrame = CSV.read(x,DataFrame,
        missingstring=ms,
        silencewarnings=true,
        types = Float64,
        normalizenames=true,
        drop=(i, nm) -> i == 4) #|> dropmissing
        dropmissing!(df,1)
        df.YY=map(x ->Int(x),df.YY);
        df.MM=map(x ->Int(x),df.MM);
        df.DD=map(x ->Int(x),df.DD);
        df.date = Date.(string.(df.YY,"-",df.MM,"-",df.DD),"yyyy-mm-dd");
        df=df[:,Not(1:3)]
        DataFrames.metadata!(df, "filename", x, style=:note);
        for x in names(df)
            if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
            end
        end
        return df 
    end

    """
    looks for control file in all xmls, takes last one, and splits
    """
    function ctl2()
        # Loop through the current directory and its subdirectories
        matches::Vector{Any} = []
        for (root, dirs, files) in walkdir(".")
        # Loop through each file name
        for file in files
            # If the file name ends with .xml
            if endswith(file, ".xml")
            path = joinpath(root, file)
            # Open the file for reading
            open(path) do f
                # Loop through each line of the file
                for line in eachline(f)
                # If the line contains 'compiling symbols in control file '
                if occursin("compiling symbols in control file ", line)
                    # Split the line by whitespace and get the fields from index 9 to 15
                    fields = split(line)[8:end] #," "
                    # Join the fields by space and print them
                    println(join(fields, " "))
                    out=join(fields, " ")
                    push!(matches,out)
                end
                end
            end
            end
        end
        end
        fl = last(matches)
        fl = split(fl)|>last
        fl = split(fl,"\"")|>first  #[2]
        return(string(fl))
    end

    """
    greps from current dir iRegex
    """
    function glob(x::AbstractString)
        filter(file -> occursin(Regex(x,"i"),file), readdir())
    end

    """
    greps from current dir iRegex
    """
    function glob(x::Regex)
        filter(file -> occursin(x,file), readdir())
    end

    function routeg(input_file::String, output_file::String)
        open(output_file, "w") do output
            line_num = 0
            in_range = false

            for line in eachline(input_file)
                line_num += 1

                if line_num > 50 && contains(line, "routing_model")
                    in_range = true
                    line = replace(line, r"ß" => "ss")
                    # line = replace(line, r"[\/]" => "_")
                    # line = replace(line, r"_" => "-")
                    line = replace(line, r"[,,]" => "")
                    line = replace(line, r"\xc4" => "Ae")
                    line = replace(line, r"\xd6" => "Oe")
                    line = replace(line, r"\xdc" => "Ue")
                    line = replace(line, r"\xe4" => "ae")
                    line = replace(line, r"\xf6" => "oe")
                    line = replace(line, r"\xfc" => "ue")
                    line = replace(line, r"\xdf" => "ss")
                    println(output, line)
                elseif in_range && contains(line, "timeoffset")
                    in_range = false
                    line = replace(line, r"ß" => "ss")
                    # line = replace(line, r"[\/]" => "_")
                    # line = replace(line, r"_" => "-")
                    line = replace(line, r"[,,]" => "")
                    line = replace(line, r"\xc4" => "Ae")
                    line = replace(line, r"\xd6" => "Oe")
                    line = replace(line, r"\xdc" => "Ue")
                    line = replace(line, r"\xe4" => "ae")
                    line = replace(line, r"\xf6" => "oe")
                    line = replace(line, r"\xfc" => "ue")
                    line = replace(line, r"\xdf" => "ss")
                    println(output, line)
                elseif in_range
                    line = replace(line, r"ß" => "ss")
                    line = replace(line, r"[,,]" => "")
                    line = replace(line, r"\xc4" => "Ae")
                    line = replace(line, r"\xd6" => "Oe")
                    line = replace(line, r"\xdc" => "Ue")
                    line = replace(line, r"\xe4" => "ae")
                    line = replace(line, r"\xf6" => "oe")
                    line = replace(line, r"\xfc" => "ue")
                    line = replace(line, r"\xdf" => "ss")
                    println(output, line)
                end
            end
        end
    end

    """
    with lapply on all qoutjl files...
    """
    function ggofjl()
        println("batch R Script for GOF")
            run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof5.R"`)
    end

    function ggofjl_nosvg()
        """
        with lapply on all qoutjl files...
        """
        println("batch R Script for GOF")
            run(`cmd.exe /c Rscript "D:/Fernerkundungsdaten/Klassifikation/R-Sessions/gof6.R"`)
    end

    """
    topdown dirsize
    """
    function du(;cwd=pwd())
        osize = 0
        n = 0
        for (root, dirs, files) in walkdir(cwd)
        for file in files
            osize += stat(joinpath(root, file)).size
            n += 1
        end
        for dir in dirs
            printstyled("check dir: $dir\n",color=:light_red)
        end
        end 
        println("$(n) files in directory")
        #@printf("%-40s %15.2f MB\n","$(cwd):",osize/1024^2)
        # printstyled("$(cwd):", " "^(40-length(cwd)), " ",
        #     round(osize/1024^2, digits=2), " MB\n",color=:light_green)
        printstyled(lpad("$(cwd):", 40), " ",
            round(osize/1024^2, digits=2), " MB\n", 
            color=:light_green)
    end
    
    function reorder_df(df::DataFrame)
        """
        date to last position
        """
        df = hcat(df[!,Not(Cols(r"date"))],df[:,Cols(r"date")])
        return(df)
    end

    function yrsum(x::DataFrame)
        df = copy(x)
        y = filter(x->!occursin("date",x),names(df))
        s = map(y -> Symbol(y),y)
            ti = try
            DataFrames.metadata(df)|>only|>last|>basename
        catch
            @warn "No basename in metadata!"
            ti = raw""
        end
        df[!, :year] = year.(df[!,:date]);
        df_yearsum = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        return(df_yearsum)
    end

    function te(x::String)
        df = dfr(x)
        dy = yrsum(df)
        for row in eachrow(dy)
            pot=row[2]
            real=row[3]
            year=row[1]
            ref=row[end]
            ref=round((ref/365)*100; digits=2)
        
        printstyled(rpad("Tdiff for 100 days in $year is $ref [mm]", 45)        
        ,color=:magenta)
        if ref <= 0
            println("\n")
        elseif ref <= 5
            println("conditions are very moist")
        elseif ref <= 10
            println("conditions are moist")
        elseif ref <= 15
            println("conditions are rather moist")
        elseif ref <= 20
            println("conditions are quite moist")
        elseif ref <= 30
            println("conditions are quite dry")
        elseif ref <= 40
            println("conditions are rather dry")
        elseif ref <= 50
            println("conditions are dry")
        elseif ref <= 70
            println("conditions are very dry")
        elseif ref <= 7e10
            println("conditions are exceptionally dry")
        end
    end
    end

    """
    new tdiff
    """
    function tdiff()
        npot = waread(r"^so_pot_trans")
        nreal =  waread(r"^so_real_trans")
        td = innerjoin(npot,nreal,on=:date)
        td = reorder_df(td)
        td.Tdiff = td[!,1] .- td[!,2]
        df = copy(td)
        y = filter(x->!occursin("date",x),names(df))
        df[!, :year] = year.(df[!,:date]);
        dy = DataFrames.combine(groupby(df, :year), y .=> sum .=> y);
        for row in eachrow(dy)
            pot=row[2]
            real=row[3]
            year=row[1]
            ref=row[end]
            ref=round((ref/365)*100; digits=2)
            printstyled(rpad("Tdiff for 100 days in $year is $ref [mm]", 45)        
            ,color=:magenta)
            if ref <= 0
                println("\n")
            elseif ref <= 5
                println("conditions are very moist")
            elseif ref <= 10
                println("conditions are moist")
            elseif ref <= 15
                println("conditions are rather moist")
            elseif ref <= 20
                println("conditions are quite moist")
            elseif ref <= 30
                println("conditions are quite dry")
            elseif ref <= 40
                println("conditions are rather dry")
            elseif ref <= 50
                println("conditions are dry")
            elseif ref <= 70
                println("conditions are very dry")
            elseif ref <= 7e10
                println("conditions are exceptionally dry")
            end
        end
        return dropmissing(td)
    end

    """
    kge barplot 
    """
    function kgeval(ds::DataFrame)
        ds.name=map(x->replace(x,r"-qoutjl.*" =>"","_" => " "),ds.name)
        ann = map(x->string.(round(x;sigdigits=2)),ds.KGE)
        p1 = Plots.bar(ds.name, ds.KGE, xlabel = "Name", ylabel = "KGE", legend = false, 
            title = splitpath(pwd())|>last, 
            xrotation = 45, 
            #fmt = :png,
            ylims = (extrema(ds.KGE) .+ [-1.0, 0.5]),  
            size = (800, 600), 
            fillcolor = ifelse.(ds.KGE .> 0, "cornflowerblue", "coral2"),
            #annotations = (ds.name,ds.KGE, ann, :top),
            xaxis = "",
            left_margin = 10mm,
            bottom_margin = 15mm, 
            bar_width = 0.6);

        for i in 1:length(ds.name)
            Plots.annotate!(ds.name[i],ds.KGE[i],(ann[i],11,
                :center,:top,:black))
            #println(ann[i]*" added")
        end
        return p1
    end

    """
    nse barplot with all values
    """
    function nsevalraw(ds::DataFrame)
        ds.name=map(x->replace(x,r"-qoutjl.*" =>"","_" => " "),ds.name)
        ann = map(x->string.(round(x;sigdigits=2)),ds.NSE)
        p1 = Plots.bar(ds.name, ds.NSE, 
            xlabel = "Name", ylabel = "NSE", legend = false, 
            title = splitpath(pwd())|>last, 
            xrotation = 45,
            ylims = (extrema(ds.NSE) .+ [-1.0, 0.5]), 
            #fmt = :png, 
            size = (800, 600), 
            fillcolor = ifelse.(ds.NSE .> 0, "cornflowerblue", "coral2"),
            xaxis = "",
            left_margin = 10mm,
            bottom_margin = 15mm, 
            bar_width = 0.6);

        for i in 1:length(ds.name)
            Plots.annotate!(ds.name[i],ds.NSE[i],
            (ann[i],11,
                :center,:top,:black))
            
        end
        return p1
    end


    function kgegrep()
        #files = glob(r"_output.txt|_outputjl") #non-recurse
        for file in filter(file -> endswith(file, "_output.txt"), readdir())
            match = grep(r"KGE", readlines(file))
            if !isempty(match)
                fn = first(split(file, "_qout"))
                for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                    line = strip(line)  # remove leading and trailing whitespace
                    line = join(split(line), " ")  ##remove inner whitespaces
                    printstyled(rpad("$fn:", 45), lpad("$line\n", 10), color = :green)
                end
            end
        end
    end


    function dfkge(dfs::Vector{Any})
        v = []
        for x in dfs
            df = reorder_df(x)
            dropmissing!(df)
            simulated = df[:,1]
            observed  = df[:,2]
            kge_value = kge1(simulated,observed)
            nse_value = nse(simulated,observed)
            nm = split(basename(last(collect(DataFrames.metadata(df)))[2]),"-")[1]
            println(replace("KGE value is $kge_value on $nm", "\\"  => "/"))
            printstyled(replace("NSE value is $nse_value on $nm\n", "\\"  => "/"),color=:green)
            push!(v,Dict(:KGE=>kge_value,:NSE=>nse_value,:name=>nm))
            v = DataFrame(v)
        end
        return(v)
    end

    """
    getnames(dfs::Vector{Any})
    get metadata names from dfs
    """
    function getnames(dfs::Vector{Any})
        try 
            nms=[]
            for df in dfs
                x=collect(DataFrames.metadata(df))[1][2]|>basename
                println(x)
                push!(nms,x)
            end
            return nms
        catch
            @error "no metadata in $df !"
            return
        end
    end

    """
    automatically sums if only datecolumn is available
    """
    function baryrsum(df::DataFrame)
        v = map(
            (x->occursin(r"date", x) & !occursin(r"year", x)),
            (names(df))
            )

        
        if any(v)
            df = yrsum(df) 
        end
        
        s = Symbol.(filter(x -> !(occursin(r"year|date", x)), names(df)))
        ti = try
            z=DataFrames.metadata(df)|>only|>last|>basename
            basename(pwd())*" $z"
        catch
            @warn "No basename in metadata!"
            ti = "Series of "*basename(pwd())
        end
        @df df groupedbar(df.year,cols(s), 
        legend = :outertopright,
        xticks = df.year,
        size = (1200,800) .* 0.7,
        xrotation = 45,
        xlabel = "", ylabel = "[mm]", 
        title = ti)
    end


   
end #endof funcs



if length(ARGS) == 0
    println("need path for controlfile <infile> <optional: moving xmls with mvwasim2 to current dir>...")
    exit()
elseif !isfile(ARGS[1])
    println("file not readable!")
    exit()
else
    infile = ARGS[1]
end

if length(ARGS) == 2
    println("performing mvwasim2 from folder wasimvzo64_10.06.05...")
    mvwasim2()
end


try
    ca = grep("inpath_hydro",readlines(infile))[end]
    printstyled(ca*"\n",color=:green)        
catch
    @warn "inpath_hydro not found in $infile !"
end

dm = pwd()|>splitpath|>last

begin 
    ofl = "route.txt"
    routeg(infile, ofl)
    #sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
        #sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
    sfn = grep("inpath",readlines(ofl))[1]
    sfn = replace.(split(sfn,"//")[2],r"\s+" => " ")
    sfn = split(sfn," ")[1]|>strip
    
    sfpt ="/mnt/d/Wasim/Tanalys/DEM/Input_V2/meteo/"
    #sfpt ="app/meteo/"
    specfile=joinpath(sfpt,sfn)
    obs = readdf(specfile)
    df = CSV.read(ofl,DataFrame,header=false,
        skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"-" => "_"),df[:,3])
    df.name=map(x->replace(x,r"_>.*" => ""),df.name)
    sort!(df, :sim)
    sim = glob(r"qges")|>first|>readdf  #waread(r"qges")
    #map(x->rm(x),glob("qoutjl"))     #rm all 
    @info "
    taking prefix of sim and colnumber -> more robust merge with regex of obs!"
    outd = []
    for i in eachrow(df)
        println(i[1],"->",i[3])
        try
            xm =innerjoin(
                
                sim[!,Cols("C"*string(i[1]),end)],

                obs[!,Cols(Regex(i[3],"i"),end)],    
                on=:date)

            onam = i[3]*"-qoutjl"
            wawrite(xm,onam)
            println("$onam saved!")
            DataFrames.metadata!(xm, "filename", onam, style=:note);
            push!(outd,xm)
            println(names(xm)," on $onam pushed!")
        catch
            onam = i[3]*"-qoutjl"
            @warn "merge is empty for $onam ! ..."
            continue
        end
    end
end


# fns = glob(r"qoutjl$")
# if length(fns) == 0
#     println("no qoutjl files found!")
#     exit()
# else
#     println(fns)
# end


if length(outd) == 0
    println("no qoutjl files processed... exiting")
    exit()
else
    onames = getnames(outd)
    println(onames)
end
#xo = readall(r"qoutjl$")
#outd=xo

for (cnt,i) in enumerate(outd)
    outname=onames[cnt]*"-lin.svg"     #"-lin.png"
    Plots.savefig(ftplin(i),outname)
    println("saved: $outname")
end


ds = dfkge(outd)
savefig(kgeval(ds),"kgeplot-$dm.svg")
savefig(nsevalraw(ds),"nse-$dm-raw.svg")

# #this is @rgof but takes long.. moved to wsl call
# printstyled("loading RCall now \n",color=:green)
# @time using RCall
# @rimport hydroGOF
# @rimport data.table as dt

# function process_files()
#     nms = filter(file -> occursin(r"qoutjl$", file), readdir())
#     for fcnt in 1:length(nms)
#         fn = nms[fcnt]
#         xj = convert(DataFrame, dt.fread(fn))
#         xj.date = Date.(string.(xj.YY,"-",xj.MM,"-",xj.DD),"yyyy-mm-dd")
#         select!(xj,Not(1:4))

#         if size(xj)[1] <= 5
#             @error("abort: input file has less than 5 lines")
#             return
#         end
#         Rgof = hydroGOF.gof(sim = xj[:, 1], obs = xj[:, 2])
#         println("GOF of Basin $(names(xj)[1]) (sim) to $(names(xj)[2])...")
#         printstyled(Rgof, bold = true, color = :blue)
#         output_file = nms[fcnt] * "_output.txt"
#         open(output_file, "w") do fy    #f fails!
#             redirect_stdout(fy) do
#                 try
#                     println("GOF of Basin $(names(xj)[2]) (sim) to $(names(xj)[end])...")
#                     println("Performance indices obs|sim:")
#                     println(Rgof)
#                 catch e
#                     @error("$e go to next")
#                 end
#             end
#         end
#     end
#     println("..the end")
# end
# process_files()
# kgegrep()
#kgewrite()
#gofbatch()
#gofbatch_nosvg()
#ggofjl_nosvg()
#ggofjl()

##run r in wsl 18
# sc="/mnt/c/Users/Public/Documents/Python_Scripts/rfile/ggof-mon.R"
# run(`wsl -d Ubuntu-18.04 -e Rscript $sc`)
#cd(raw"D:\Wasim\regio\out\rc200\x22\reloc")

odat = grep("otherdata",readdir())
if  length(odat)>0 
    #println("otherdata found!")
    subpt="/mnt/c/Users/Public/Documents/Python_Scripts/bashscripts/sosub.sh"
    odat = first(odat)
    #run(`wsl.exe -d Ubuntu-18.04 -e bash $subpt $odat`)
    run(`bash $subpt $odat`)
end


sc="/mnt/c/Users/Public/Documents/Python_Scripts/rfile/gofdf.R"
run(`wsl.exe -d Ubuntu-18.04 -e Rscript $sc`) #in ubu wsl


########dthydro.py############################################
begin 
    pypt=dirname(src_path)
    m = "qoutjl"
    run(pipeline(`python $pypt/dthydro.py $m`, "$dm-hydro.txt"))
end
#output to screen
lines = readlines("$dm-hydro.txt");
for line in lines
    println(line)
end
###############waba############################################
wpth=src_path*"/waba-rev.jl"

try 
    include(wpth)
catch e
    @error "waba failed! \n $e"
end

###############tdiff############################################
try
    xdt = tdiff();
    wawrite(xdt,"tdiff-$dm-jl.txt")
    @info "tdiff-$dm-jl.txt saved!"
    Plots.theme(:dao)
    ptd = baryrsum(xdt);
    savefig(ptd,"tdiff-$dm-jl.png")
    @info "tdiff-$dm-jl.png saved!"
catch e
    @error "$e \n"
    @error "tdiff failed! \n"
end


###############cleanup############################################
# include("$pypt/julia/rmeq.jl")
# ncrem=src_path*"/ncremover.jl";
# try 
#     include(ncrem)
# catch e
#     @error "ncrem failed! \n $e"
#     #print current dir + size
#     du()
# end

#printstyled("done!\n",color=:green)
kgegrep()