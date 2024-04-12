#"D:/Wasim/regio/out/rc200/x22/loc5/l5spin/"|>cd
#raw"D:\Wasim\regio\out\rc200\x22\loc5\nspin"|>cd
# if length(ARGS) == 0
# 	println("need args! <file>...")
#     exit()
# end
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

@time using DataFrames,CSV,Dates



#"D:/Wasim/regio/out/rc200/x22/loc5/l5spin/"|>cd
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
        df = CSV.read(x, DataFrame; delim="\t", header=1, missingstring=ms, normalizenames=true, types=Float64)
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
        #ignorerepeated = true,
        #delim="\t",
        #skipto=4,
        types = Float64,
        silencewarnings=false,
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
            # Join the root and file name to get the full path
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

    function glob(x::AbstractString)
        """
        greps from current dir iRegex
        """
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

    function te(x::String)
        df = dfr(x)
        dy = yrsum(df)
        for row in eachrow(dy)
            pot=row[2]
            real=row[3]
            year=row[1]
            ref=row[end]
            ref=round((ref/365)*100; digits=2)
        #printstyled("Tdiff for 100 days in $year is $ref [mm]\t|\t",color=:magenta)
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

    function tdiff()
        npot = glob(r"^so_pot_trans")|>first|>readdf
        nreal =  glob(r"^so_real_trans")|>first|>readdf
        td = innerjoin(npot,nreal,on=:date)
        td = reorder_df(td)
        td.Tdiff = td[!,1] .- td[!,2]
        te(td)
        return dropmissing(td)
    end

end #end of funcs

###########
#qba()
#@edit waread(x)
#mvwasim2(;pt=raw"C:\Users\chs72fw\Documents\EFRE_GIS\Hydrologie\WaSiM-ETH\wasimvzo64_10.07.02",force=true)

#mvwasim2()

#@doc mvwasim2(;pt=raw"C:\Users\chs72fw\Documents\EFRE_GIS\Hydrologie\WaSiM-ETH\wasimvzo64_10.07.02")


#infile = ctl2()

# infile = replace(infile,#r"ctl.*"=>"ctl",
# "control"=>"D:/Wasim/regio/control")

if length(ARGS) == 0
    println("need path for controlfile <infile>...")
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

import Grep: grep
try
    ca = grep("inpath_hydro",readlines(infile))[end]
    printstyled(ca*"\n",color=:green)        
catch
    @warn "inpath_hydro not found in $infile !"
end

dm = pwd()|>splitpath|>last

try
    xdt = tdiff();
    wawrite(xdt,"tdiff-$dm-jl.txt")
    @info "tdiff-$dm-jl.txt saved!"
catch #e
#    @error "$e \n
    @error "tdiff failed! \n"
end

begin 
    ofl = "route.txt"
    routeg(infile, ofl)
    sfn = readlines(ofl)[6]|>split|>first|>k->split(k,"/")|>last
    sfpt ="D:/Wasim/Tanalys/DEM/Input_V2/meteo/"
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
            dm =innerjoin(
                
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
end

fns = glob("qoutjl")

if length(fns) == 0
    println("no qoutjl files found!")
    exit()
else
    println(fns)
end

#gofbatch()
#gofbatch_nosvg()
#ggofjl_nosvg()
#ggofjl()
#kgegrep()
##run r in wsl 18
sc="/mnt/c/Users/Public/Documents/Python_Scripts/rfile/ggof-mon.R"
#run(`wsl -d Ubuntu-18.04 -e awk -h`)
run(`wsl -d Ubuntu-18.04 -e Rscript $sc`)
# using RCall
# R"""
# """

begin 
    #pypt="C:/Users/Public/Documents/Python_Scripts"
    pypt=dirname(src_path)
    dm = pwd()|>splitpath|>last
    m = "qoutjl"
    #run(`mamba run -n ds python $pypt/dthydro.py $m`)
    #output to file
    run(pipeline(`mamba run -n ds python $pypt/dthydro.py $m`, "$dm-hydro.txt"))
end
#output to screen
lines = readlines("$dm-hydro.txt");
for line in lines
    println(line)
end
#waba
#wpth=src_path*"/water-balance.jl"
wpth=src_path*"/waba-rev.jl"

try 
    include(wpth)
catch e
    @error "waba failed! \n $e"
end

include("$pypt/julia/rmeq.jl")
ncrem=src_path*"/ncremover.jl";
try 
    include(ncrem)
catch e
    @error "ncrem failed! \n $e"
    #print current dir + size
    du()
end

#printstyled("done!\n",color=:green)



# using Conda
# @doc Conda.runconda()
# cmds = (`run python $pypt/dthydro.py $m`)
# Conda.runconda(cmds)

#pipeline("tst",Conda.runconda(cmds))

# txt = sprint(io -> Conda.runconda(cmds))
# open("delim_file.txt", "w") do io
#     writedlm(io, txt)
#     #writedlm(io,sprint(Conda.runconda(cmds)))
# end

# function run_command_and_save(cmds::Cmd, output_file::AbstractString)
#     # Run the command and capture the output
#     #output = read(pipe(cmds, "r"))
#     #output = read(pipeline(Conda.runconda(cmds)), String)

#     # Display the output on the terminal
#     println(output)

#     # Write the output to a text file
#     open(output_file, "w") do io
#         writedlm(io, output)
#     end
# end

# # Example usage:
# cmds = (`run python $pypt/dthydro.py $m`)
# output_file = "delim_file.txt"

# run_command_and_save(cmds, output_file)

