#functions
#using DataFrames, CSV, Statistics, Dates, StatsPlots, Distributions, Printf
#using DataFrames, CSV, Statistics, Dates, Rasters, StatsPlots, Distributions
using DelimitedFiles, Grep, Printf, Statistics, Dates
#for plots...
#default(show = true)
#plotlyjs()

function ssup()
    include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/smallfuncs-ubu.jl")
end



#source functions
#include("/mnt/c/Users/Public/Documents/Python_Scripts/julia/func-win.jl")
# for file in readdir(dir_path)
#     if occursin(file_ending, file)

function ls()
    p=pwd()
    f=readdir()
    dirs=filter(x -> (isdir(x)),f)
    files=filter(x -> (isfile(x)),f)
    if length(files)>0 && length(dirs)>0
        nf=length(files)
        println("$p\ndirs: $dirs\n $nf files:\n$files")
    elseif length(dirs)>0
        println("$p\ndirs: $dirs\n")
    elseif (length(files)>0 && length(files)<=12)
        nf=length(files)
        println("$p\n$nf files:\n $files\n")
    elseif (length(files)>0 && length(files)>12)
        nf=length(files)
        println("$p\n$nf files:\n",first(f,6),"...",last(f,6))
    else
        println("$p\n")
    end
end

function wcl(x::AbstractString)
    files = filter(file -> occursin(Regex(x,"i"),file), readdir())
    for file in files
        open(file) do f
            ct=(count(_ -> true, eachline(f)))
            #println(file,ct)
            println("$file:\t $ct")
        end
    end
end
#wcl("evar")
#wcl("Ev")
#r"so_snow_evaporation"|>dfp
#wcl(glob("wq")|>first,true)


#faster and pure julia:
function vgrep(regex, file_ending)
    files = filter(file -> endswith(file, file_ending), readdir())
    # loop over each file
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                # check if the line matches the regex
                if occursin(Regex(regex), line)
                    #m=count(_ -> true, line) #das zählt die linechars
                    println("$file: $counter:\t $line")
                end
            end
        end
    end
end

function vg(snippet::AbstractString, file_ending::AbstractString)
    files = filter(file -> endswith(file, file_ending), readdir())
    # loop over each file
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
end

function vgjl(snippet::AbstractString)
    owd=pwd()
    cd("/mnt/c/Users/Public/Documents/Python_Scripts/julia")
    files = filter(file -> endswith(file, ".jl"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(owd)
end

function vgro(snippet::AbstractString)
    owd=pwd()
    cd("D:/Fernerkundungsdaten/Klassifikation/R-Sessions")
    files = filter(file -> endswith(file, ".R"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(owd)
end

function vgr(snippet::AbstractString)
    owd="D:/Fernerkundungsdaten/Klassifikation/R-Sessions"
    files = filter(file -> endswith(file, ".R"), readdir(owd,join=true))
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled(rpad("$file:",50),color=:light_magenta)
                    #printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled(lpad("$line\n",30),color=:green,bold=true)
                    #underline = true 
                end
            end
        end
    end
end

#"dySeries"|>vgr

function vgpyo(snippet::AbstractString)
    owd=pwd()
    cd("/mnt/c/Users/Public/Documents/Python_Scripts")
    files = filter(file -> endswith(file, ".py"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(owd)
end

function vgpy(snippet::AbstractString)
    owd="/mnt/c/Users/Public/Documents/Python_Scripts"
    files = filter(file -> endswith(file, ".py"), readdir(owd,join=true))
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled(rpad("$file:",50),color=:light_magenta)
                    #printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled(lpad("$line\n",30),color=:green,bold=true)
                    #underline = true 
                end
            end
        end
    end
end

function vgctl(snippet::AbstractString)
    """
    hint. vgctl("set \$TS")
    """
    owd=pwd()
    nwd="D:/Wasim/regio/control/"
    nwd2="D:/temp/saale/control/"
    nwd3="D:/Wasim/Tanalys/DEM/brend_fab/control/"
    nwd4="D:/Wasim/regio/control/"
    nwd5="D:/Wasim/streu/control/"
    cd(nwd)
    #println("greps from *ctl from  \n$nwd and \n$nwd2...")
    println("greps from *ctl from  \n$nwd \n$nwd2 \n$nwd3 \n$nwd4 \n$nwd5...")
    files = filter(file -> endswith(file, ".ctl"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter: $nwd",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(nwd2)
    files = filter(file -> endswith(file, ".ctl"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter: $nwd2",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(nwd3)
    files = filter(file -> endswith(file, ".ctl"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter: $nwd3",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(nwd4)
    files = filter(file -> endswith(file, ".ctl"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter: $nwd4",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(nwd5)
    files = filter(file -> endswith(file, ".ctl"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if contains(line,snippet)
                    printstyled("$counter: $nwd5",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
    cd(owd)
end

function rglob(prefix::AbstractString)
    rootdir="."
    results = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
            if (occursin(Regex(prefix,"i"),filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    return results
end

function rglob(prefix::Regex)
    rootdir="."
    results = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (occursin(prefix,filename))
                push!(results, joinpath(looproot, filename)) 
            end
        end
    end
    return results
end

#go dir up
function cdb()
    dirname(pwd())|>cd
    pwd()|>println
end

#like jdd to vector of strings.
function fdd()
    cwd = pwd()
    dirs = readdir(".")
    s = []
    for dir in dirs
        if isdir(dir)
            push!(s,joinpath(cwd, dir))
            size = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    size += stat(joinpath(root, file)).size
                end
            end
	    @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
#        else
#	    @printf("%-40s\n","$(cwd)");
	end
    end
    return(s)
end

#fdd()


function vgr(regex, file_ending)
    rootdir=pwd()
    println("starting on: $rootdir...\n searching for >> $regex << with file ending >> $file_ending <<\n")
    files = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
            #if (occursin(Regex(prefix,"i"),filename))
            if (endswith(filename, file_ending))
                push!(files, joinpath(looproot, filename)) 
            end
        end
    end
    #files = filter(file -> endswith(file, file_ending), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if occursin(Regex(regex,"i"), line)
                    println("$file: $counter:\t $line")
                end
            end
        end
    end
end

function vjl(regex)
    # greps jl from julia folder
    pt="/mnt/c/Users/Public/Documents/Python_Scripts/julia";
    file_ending=".jl"
    files = filter(file -> endswith(file, file_ending), readdir(pt,join=true))
    for file in files
        open(file) do f
            counter = 0
            for line in eachline(f)
                counter += 1
                if occursin(Regex(regex,"i"), line)
                    println("$file: $counter:\t $line")
                end
            end
        end
    end
end


function cnt()
    return(length(readdir(pwd())))
end

function du()
    cwd = pwd()
    n = length(readdir(cwd))
    osize = 0
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
     end
    end 
    println("$(n) files in directory")
    @printf("%-40s %15.2f MB\n","$(cwd):",osize/1024^2)
end 

function route_from_dir(dir::String)
    dirs = readdir(dir)
    routes::Vector{String} = []
    for directory in dirs
        if isfile("$dir/" * directory)
            push!(routes, "$dir/$directory")
        else
            if ~(directory in routes)
                newread = dir * "/$directory"
                newrs = route_from_dir(newread)
                [push!(routes, r) for r in newrs]
            end
        end
    end
    routes
end


function llf()
    cwd = pwd()
    #n = length(readdir(cwd))
    v = []
    v_in_subdir = []
    for (root, dirs, files) in walkdir(cwd)
        for file in files
            file_path = joinpath(root, file)
            if isfile(file_path) && (!occursin(r"xml|qgk|fzt|ftz|log|ini|wq|txt|yrly|nc|tif|jpg|png|svg",file))
                push!(v,file_path)
            elseif isdir(file_path)
                v_in_subdir = llf()
                println("is dir ",file_path)
            end
        end
        #return(vcat(v,v_in_subdir))
        return(v)
    end
end
        

function vg2(regex::AbstractString, ending::AbstractString)
    cmd = `grep --color=always -C2 -rIHn -E "$regex" --include="*.$ending"`
    run(cmd)
end


function tff2(x::Vector{String})
    for filename in x
        b = Dict{String, Float64}()
        m = Dict{String, Float64}()
        h = Dict{String, Float64}()
        cnte = Dict{String, Int64}()
        ncols = 0

        open(filename) do file
            first_line = true
            

            for (i, line) in enumerate(eachline(file))
                fields = split(line)
                if first_line
                    println("filename: $filename")
                    ncols = length(fields)
                    println("no of fields: $ncols")
                    println("year\t$(join(fields[5:end], "\t"))")
                    first_line = false
                    continue
                end

                if match(r"^\d{4}$", fields[1]) != nothing
                    cnte[match(r"^\d+", fields[1]).match] = get(cnte, match(r"^\d+", fields[1]).match, 0) + 1
                    b[fields[1]] = get(b, fields[1], 0.0) + parse(Float64, replace(fields[5], r"-9999" => "0"))
                    m[fields[1]] = get(m, fields[1], 0.0) + parse(Float64, replace(fields[end-1], r"-9999" => "0"))
                    h[fields[1]] = get(h, fields[1], 0.0) + parse(Float64, replace(fields[end], r"-9999" => "0"))
                end
            end
        end

        if ncols <= 5
            for key in sort(collect(keys(b)))
                println("$key\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
            end
        elseif ncols == 6
            for key in sort(collect(keys(b)))
                println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
            end
        elseif ncols >= 7
            for key in sort(collect(keys(b)))
                println("$key\t$(@sprintf("%.2f", b[key]))\t$(@sprintf("%.2f", m[key]))\t$(@sprintf("%.2f", h[key]))\t| means: $(@sprintf("%.2f", b[key] / cnte[key]))\t$(@sprintf("%.2f", m[key] / cnte[key]))\t$(@sprintf("%.2f", h[key] / cnte[key]))\t| counts: $(cnte[key])")
            end
        end
    end
end

function glob(x::AbstractString)
    """
    greps from current dir iRegex
    """
    filter(file -> occursin(Regex(x,"i"),file), readdir())
end

function glob(x::Regex)
    """
    greps from current dir Regex
    """
    filter(file -> occursin(x,file), readdir())
end

function fdi()
    #    cwd==nothing
    #    dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
        cwd = pwd()
        dirs = readdir(cwd)
        for dir in dirs
            if isdir(dir)
                @printf("%-8s\t|","$dir");
            end
        end
    end
    
    function fdi(cwd::AbstractString)   
        dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
        for dir in dirs
            if isdir(dir)
                @printf("%-8s\t|","$dir");
            end
        end
    end
    
function jdd()
    cwd = pwd()
    dirs = readdir(".")
    for dir in dirs
        if isdir(dir)
            size = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    size += stat(joinpath(root, file)).size
                end
            end
	    @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
        end
    end
end

function dd()
    cwd = pwd()
#    dirs = readdir(".")
    osize = 0
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
     end
    end 
    @printf("%-40s %15.3f GB\n","$(cwd):",osize/1024^3);
end 

function ct(ext::AbstractString)
    cwd = pwd() 
    osize = 0
    fz = 0
    m = []
    for (root, dirs, files) in walkdir(cwd)
     for file in files
     if isfile(file) && occursin(Regex(ext),file)
	 nm=joinpath(root, file)
	 osize = stat(nm).size
	 @printf("%-40s %15.2f MB\n","$(nm):",osize/1024^2);
	 fz += stat(nm).size
	 push!(m,(nm))
     end
    end 
    end 
     n=repeat(" - -",10)
     println(n*" sum of ",ext*n)
     @printf("%-40s %15.2f MB\n","$(cwd):",fz/1024^2);
     println(n,length(m)," matches "*n,"\n")
     return(m)
end 

function ct()
    cwd = pwd() 
    osize = 0
    fz = 0
    m = []
#    sizes = Dict()
    for (root, dirs, files) in walkdir(cwd)
     for file in files
     if isfile(file)
        nm=joinpath(root, file)
        osize = stat(nm).size
        #sizes[file] = stat(nm).size/1024^2
	    @printf("%-40s %15.2f MB\n","$(nm):",osize/1024^2);
	    fz += stat(nm).size
	    push!(m,(nm))
     end     
    end 
    end 
    # sort(df, [order(:a), order(:b, rev = true)]) 
    # sorted_files = sort(collect(keys(sizes)), by=x->sizes[x], rev=false)
    # println(sorted_files)
     n=repeat(" - -",10)
     println(n*" total sum ")
     @printf("%-40s %15.2f MB\n","$(cwd):",fz/1024^2);
     println(n,length(m)," files present "*n,"\n")
     return(m)
end

function print_sorted_sizes(dir)
    folders = [joinpath(dir, f) for f in readdir(dir)]
    sizes = Dict()
    for f in folders
        if isdir(f)
            sizes[f] = get_folder_size(f)
        end
    end
    sorted_folders = sort(collect(keys(sizes)), by=x->sizes[x], rev=false)
    for f in sorted_folders
        if sizes[f] >= 1000000
        printstyled(rpad(f,60, ' '), rpad(sizes[f] ÷ 10^6, 6, ' '), "MB\n",color=:green)
        end
    end
end

function get_folder_size(folder)
    files = readdir(folder)
    size = 0
    for file in files
        path = joinpath(folder, file)
        if isfile(path)
            size += stat(path).size
        elseif isdir(path)
            size += get_folder_size(path)
        end
    end
    return size
end

function fz()
    """
    gets sorted DF by size recursivley
    """
    cwd = pwd() 
    osize = 0
#    fn = 0
    m = []
    for (root, dirs, files) in walkdir(cwd)
     for file in files
        if isfile(file)
           nm=joinpath(root, file)
           osize = stat(nm).size/1024^2
#	       fn += stat(nm).size
	       #push!(m,(nm))
           push!(m,Dict(:name=>file,
           :size=>osize,
#           :total=>(fn/1024^2),
           :fullnaname=>nm))
        end
    end 
end
    df = DataFrame(m)     
    sort!(df, [order(:size,rev = true), order(:name)])
    return(df)
    # fn = fn/1024^2
    # printstyled("$fn MB",color=:blue)
end 

function jdd()
    cwd = pwd()
    dirs = readdir(".")
    for dir in dirs
        if isdir(dir)
            size = 0
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    size += stat(joinpath(root, file)).size
                end
            end
	    @printf("%-40s %15.2f MB\n","$(cwd)\\$dir:",size/1024^2);
        end
    end
end

function dd()
    cwd = pwd()
#    dirs = readdir(".")
    osize = 0
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
     end
    end 
    @printf("%-40s %15.3f GB\n","$(cwd):",osize/1024^3);
end 

function ct(ext::AbstractString)
    cwd = pwd() 
    osize = 0
    fz = 0
    m = []
    for (root, dirs, files) in walkdir(cwd)
     for file in files
     if isfile(file) && occursin(Regex(ext),file)
	 nm=joinpath(root, file)
	 osize = stat(nm).size
	 @printf("%-40s %15.2f MB\n","$(nm):",osize/1024^2);
	 fz += stat(nm).size
	 push!(m,(nm))
     end
    end 
    end 
     n=repeat(" - -",10)
     println(n*" sum of ",ext*n)
     @printf("%-40s %15.2f MB\n","$(cwd):",fz/1024^2);
     println(n,length(m)," matches "*n,"\n")
     return(m)
end 

function ct()
    cwd = pwd() 
    osize = 0
    fz = 0
    m = []
#    sizes = Dict()
    for (root, dirs, files) in walkdir(cwd)
     for file in files
     if isfile(file)
        nm=joinpath(root, file)
        osize = stat(nm).size
        #sizes[file] = stat(nm).size/1024^2
	    @printf("%-40s %15.2f MB\n","$(nm):",osize/1024^2);
	    fz += stat(nm).size
	    push!(m,(nm))
     end     
    end 
    end 
    # sort(df, [order(:a), order(:b, rev = true)]) 
    # sorted_files = sort(collect(keys(sizes)), by=x->sizes[x], rev=false)
    # println(sorted_files)
     n=repeat(" - -",10)
     println(n*" total sum ")
     @printf("%-40s %15.2f MB\n","$(cwd):",fz/1024^2);
     println(n,length(m)," files present "*n,"\n")
     return(m)
end

function print_sorted_sizes(dir)
    folders = [joinpath(dir, f) for f in readdir(dir)]
    sizes = Dict()
    for f in folders
        if isdir(f)
            sizes[f] = get_folder_size(f)
        end
    end
    sorted_folders = sort(collect(keys(sizes)), by=x->sizes[x], rev=false)
    for f in sorted_folders
        if sizes[f] >= 1000000
        printstyled(rpad(f,60, ' '), rpad(sizes[f] ÷ 10^6, 6, ' '), "MB\n",color=:green)
        end
    end
end

function get_folder_size(folder)
    files = readdir(folder)
    size = 0
    for file in files
        path = joinpath(folder, file)
        if isfile(path)
            size += stat(path).size
        elseif isdir(path)
            size += get_folder_size(path)
        end
    end
    return size
end

function fz()
    """
    gets sorted DF by size recursivley
    """
    cwd = pwd() 
    osize = 0
#    fn = 0
    m = []
    for (root, dirs, files) in walkdir(cwd)
     for file in files
        if isfile(file)
           nm=joinpath(root, file)
           osize = stat(nm).size/1024^2
#	       fn += stat(nm).size
	       #push!(m,(nm))
           push!(m,Dict(:name=>file,
           :size=>osize,
#           :total=>(fn/1024^2),
           :fullnaname=>nm))
        end
    end 
end
    df = DataFrame(m)     
    sort!(df, [order(:size,rev = true), order(:name)])
    return(df)
end 

function homg()
    pt="D:/Wasim/Goldbach/revision/"
    cd(pt)
    println("you are here: ",pwd())
end

function hombr()
    pt="D:/Wasim/Tanalys/DEM/brend_fab/out/m4/"
    cd(pt)
    println("you are here: ",pwd())
end

function hometeo()
    cd("D:/Wasim/Tanalys/DEM/Input_V2/meteo/")
    println("you are here: ",pwd())
end

function homreg()
    cd("D:/Wasim/regio/out/");
    println("you are here: ",pwd())
    fdi()
end

function homg()
    cd("D:/Wasim/Goldbach/");
    println("you are here: ",pwd())
    fdi()
end

function home()
    cd("D:/Wasim/");
    println("you are here: ",pwd())
    fdi()
end

function lf()
    """
    list_files_sorted_by_last_change
    formerly lat()
    """
    directory = pwd()
    files = readdir(directory)  
    file_times = Dict{String, Dates.DateTime}()
    for file in files
        file_path = joinpath(directory, file)
        stat_info = stat(file_path)
        file_times[file] = Dates.unix2datetime(stat_info.mtime)
    end   
    file_times = collect(file_times)
    xf = sort(file_times, by = x -> x[2], rev=true)
    sz = []                 #::Vector{Float64}
    for i in 1:length(xf)
       push!(sz,round(stat(xf[i][1]).size/1024^2,digits=4))
    end
    merged_vector = []
    for (item1, item2) in zip(xf,sz)
        merged_item = (item1.first, item1.second,item2)
        push!(merged_vector, merged_item)
    end
    df = DataFrame(merged_vector);
    rename!(df,["latest_file","modified","size_in_mb"]);
    df.modified=map(x -> Dates.format(x, "yyyy-mm-dd HH:MM"),df.modified)
    return(df)
end

##very nice meta programming... 
macro vv(s) vgjl(s);end
#@vv "unter"
macro vpy(s) vgpy(s);end
#@vpy "climate"
macro vr(s) vgr(s);end
#@vr "climate"
macro vct(s) vgctl(s);end
#@vct "das ist"
macro rg(s) rglob(s);end
macro gl(s) glob(s)|>first;end
#fastplot
macro fp(s) dfp(Regex(s));end
macro flog(s) dfl(Regex(s));end
macro ncrm() ncrem="/mnt/c/Users/Public/Documents/Python_Scripts/julia/ncremover.jl";include(ncrem);end
macro nco(s) nconly(s);end

#@ncrm #works

function fsize(m::String)
    """
    names and size in MB via readdir
    """
    files = filter(x -> occursin(Regex(m,"i"), x), readdir())
    fzs = [(file, filesize(file) / 2^20) for file in files]
    tot = round(sum(map(x->x[2],fzs));digits=3)
    #println("Total Size: $tot MB")
    printstyled("Total Size: $tot MB\n",color=:green)
    return(fzs)
end

using TOML
#TOML is short for Tom’s Obvious Minimal Language and is a configuration file format 
#that should be “easy to parse into data structures in a wide variety of languages”
macro toml_str(s::String)
    TOML.parse(s)::Dict{String, <:Any}
end

# toml"""
# ki11 = 6.6
# dr11 = 1
# kb11 = 0.4
# q011 = 0.02
# """

# dx = toml"""ki11 = 6.6 
# dr11 = 1
# kb11 = 0.4
# q011 = 0.02"""|>DataFrame
# #DataFrame(dx)
# findall
# findfirst
# findlast
# findnext
# findprev

# fx=raw"D:\Wasim\regio\control\rcm200_r6.ctl"
# for (i, line) in enumerate(eachline(fx))
#     if findfirst("dr11", line) !== nothing
#         println("Line $i: $line")
#     end
# end

# fx = raw"D:\Wasim\regio\control\rcm200_r6.ctl"
# prevline = ""
# for (i, line) in enumerate(eachline(fx))
#     if findfirst("dr11", line) != nothing
#         println("Line $(i-1): $(prevline)")
#         printstyled("Line $i: $line\n",color=:green)
#         nextline = readline(fx)
#         println("Line $(i+1): $nextline")
#     end
#     prevline = line
# end

function ctlg(dir_path::String, file_ending::String, match::String)
    for file in readdir(dir_path)
        if occursin(file_ending, file)
            fx = joinpath(dir_path, file)
            prevline = ""
            for (i, line) in enumerate(eachline(fx))
                if findfirst(match, line) !== nothing
                    printstyled("File: $file\n",color=:red)
                    println("$(prevline)")
                    printstyled("$line\n",color=:green)
                    nextline = readline(fx)
                    println("$(nextline)")
                end
                prevline = line
            end
        end
    end
end

# ctlg("D:\\Wasim\\regio\\control", "v9.local.ctl", "meteo_")

# ctlg(raw"C:\Users\Public\Documents\Python_Scripts\julia\win", 
#     "jl", "surface")


function merge_vectors(vector1::Vector{Any}, vector2::Vector{Pair{String, DateTime}})
    merged_vector = []
    for (item1, item2) in zip(vector1, vector2)
        merged_item = (item1, item2.first, item2.second)
        push!(merged_vector, merged_item)
    end
    return merged_vector
end

function merge_vectors(vector1::Vector{String}, vector2::Dict{String, DateTime})
    merged_vector = []
    for (item1, item2) in zip(vector1, vector2)
        merged_item = (item1, item2.first, item2.second)
        push!(merged_vector, merged_item)
    end
    return merged_vector
end


function latx()
    """
    list_files_sorted_by_last_change
    """
    directory = pwd()
    files = readdir(directory)  
    file_times = Dict{String, Dates.DateTime}()
    for file in files
        file_path = joinpath(directory, file)
        stat_info = stat(file_path)
        file_times[file] = Dates.unix2datetime(stat_info.mtime)
    end
    #xf = sort(file_times, by = x -> x[2], rev=true)
    
    file_times = collect(file_times)
    
    xf = sort(file_times, by = x -> x[2], rev=true)
    xf = first(xf,11)
    sz = []                 #::Vector{Float64}
    for i in 1:length(xf)
       push!(sz,round(stat(xf[i][1]).size/1024^2,digits=4))
    end
    xf2 = merge_vectors(sz,xf)
    for (size, file, dt) in xf2
        datetime = Dates.format(dt, "yyyy-mm-dd HH:MM")
        printstyled(rpad(file,35),color=:yellow),
        printstyled("$datetime\t" ,color=:green),
        printstyled(rpad(size,7)," MB\n",color=:yellow)    
    end
end

println("you are here: ",pwd(),"\nto edit this script:\n @less (ssup()) or @edit (ssup())")
fdi()

# @vv "surface"
# x=raw"/mnt/c/Users/chs72fw/Documents/Promotionsstudium/Dropbox/brendfab/m3/Soil_Temperature_Stack.nc"
# r = Raster(x)
# r = readras(x)
# plot(r[t=4:5])
# surface(r[t=5],camera = (0,-75),legend=false)
