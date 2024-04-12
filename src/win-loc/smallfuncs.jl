#functions
#using DataFrames, CSV, Statistics, Dates, StatsPlots, Distributions, Printf
#using DataFrames, CSV, Statistics, Dates, Rasters, StatsPlots, Distributions
@time using Grep, Printf, Statistics, Dates, DataFrames, CSV #, PrettyTables
import InteractiveUtils.clipboard
import DelimitedFiles 
printstyled("Grep, Printf, Statistics, Dates, 
    DataFrames, CSV loaded\n",color=:green)
#w csv 2sec #PrettyTables
begin
#    if isnothing(src_path)
        if Sys.isapple()
            platform = "osx"
            const homejl = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
            const mybash = "/Users/apfel/.bash_aliases"
            src_path = "/Users/apfel/Library/Mobile Documents/com~apple~CloudDocs/uni/GitHub/Python-Scripts/julia"
        elseif Sys.iswindows()
            platform = "windows"
            src_path = "C:\\Users\\Public\\Documents\\Python_Scripts\\julia"
            macro wasim() pt="C:\\Users\\chs72fw\\.julia\\dev\\WaSiM\\src\\wa.jl";include(pt);end
            try
                #Base.@showtime using RCall #errors in local scope
                using RCall
                #R""".libPaths(new="C:/Users/chs72fw/AppData/Local/R/win-library/4.2")"""
                #libstr="C:/Users/chs72fw/AppData/Local/R/win-library/4.2"
                #RCall.@R_str(".libPaths($libstr)")
                #printstyled("\nRCall loaded and libpath assigned!\n",color=:green,bold=true,underline=true)
                printstyled("\nR 4.2.2 loaded !\n",color=:green,bold=true,underline=true)
            catch e
                printstyled("RCall errors at:\n $e",color=:red)
            end    
        else
            platform = "unix"
            winpt = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
            #pcld = "~/pCloud Drive/Stuff/Python_Scripts/julia"
            #pcld = "/home/cris/pCloudDrive/Stuff/Python_Scripts/julia"
            pcld = "~/pCloudDrive/Stuff/Python_Scripts/julia"
            src_path = isdir(winpt) ? winpt : pcld
            #if (!isdir(src_path) && isdir(pcld))
            if !isdir(src_path) # for docker images
                src_path = "/app/pyscripts/julia"
            end
            println("sourcepath is $src_path")
            if isdir(winpt)
                macro wasim() pt="/mnt/c/Users/chs72fw/.julia/dev/WaSiM/src/wa.jl";include(pt);end
            end
        end   
#    else
#        println("src_path already set to $src_path")
#    end
end
#@time using DataFrames, CSV, StatsPlots

here = pwd()
#@showtime here = pwd()

println("\nto edit this script: 
    @less (ssup()) or @edit (ssup())\n\t",
"fdi() for dirinfo, ls() for fileinfo")

#printstyled("loc $here\n now initializing...\n",color=:green)

#for plots...
#default(show = true)
#plotlyjs()

function ssup()
    thisfile=src_path*"/smallfuncs.jl"
    #include("C:/Users/Public/Documents/Python_Scripts/julia/smallfuncs.jl")
    include(thisfile)
end

function setup()
    thisfile=src_path*"/func-win.jl"
    include(thisfile)
end

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

"""
counts lines like wc -l in bash
"""
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

#wcl("use")
#wcl("evar")
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
                if Base.contains(line,snippet)
                    printstyled("$counter:\t",color=:light_red) 
                    printstyled("$file:\t",color=:light_magenta,underline = true,blink = false,bold=true) 
                    printstyled("$line\n",color=:green,bold=true) 
                end
            end
        end
    end
end

function vgjl(snippet::AbstractString)
    owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "C:/Users/Public/Documents/Python_Scripts/julia"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
    end
    
    cd(script_dir)
    
    # files = []
    files = filter(file -> endswith(file, ".jl"), readdir())
    fwin = filter(file -> endswith(file, ".jl"), 
        readdir(script_dir*"/win", join=true))
    files = vcat(files,fwin)
    for file in files
        open(file) do f
            counter = 0 # Initialize the counter
            for line in eachline(f)
                counter += 1 # Increment the counter
                if Base.contains(line, snippet)
                    printstyled("$counter:\t", color=:light_red)
                    printstyled("$file:\t", color=:light_magenta, underline=true, blink=false, bold=true)
                    printstyled("$line\n", color=:green, bold=true)
                end
            end
        end
    end
    
    cd(owd)
end

function vgro(snippet::AbstractString)
    owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "D:/Fernerkundungsdaten/Klassifikation/R-Sessions"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions"
    end
    
    cd(script_dir)

    files = filter(file -> endswith(file, ".R"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if Base.contains(line,snippet)
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
    #owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "D:/Fernerkundungsdaten/Klassifikation/R-Sessions"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/d/Fernerkundungsdaten/Klassifikation/R-Sessions"
    end

    files = filter(file -> endswith(file, ".R"), readdir(script_dir,join=true))
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if Base.contains(line,snippet)
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

function vgpyo(snippet::AbstractString)
    owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "C:/Users/Public/Documents/Python_Scripts"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/c/Users/Public/Documents/Python_Scripts"
    end
    
    cd(script_dir)


    files = filter(file -> endswith(file, ".py"), readdir())
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if Base.contains(line,snippet)
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
    #owd = abspath(pwd())
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform
    
    if platform == "windows"
        script_dir = "C:/Users/Public/Documents/Python_Scripts"
    else
        # Assuming you want to use a different path for Linux/WSL, adjust as needed
        script_dir = "/mnt/c/Users/Public/Documents/Python_Scripts"
    end
    
    files = filter(file -> endswith(file, ".py"), readdir(script_dir,join=true))
    for file in files
        open(file) do f
            counter = 0 # Zähler initialisieren
            for line in eachline(f)
                counter += 1 # Zähler erhöhen
                if Base.contains(line,snippet)
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

"""
Usage: vgctl("set \$TS")
"""
function vgctl(snippet::AbstractString)
    owd = pwd()
    platform = Sys.iswindows() ? "windows" : "linux"  # Check the platform

    if platform == "windows"
        nwd = "D:/Wasim/regio/control/"
        nwd2 = "D:/temp/saale/control/"
        nwd3 = "D:/Wasim/Tanalys/DEM/brend_fab/control/"
        nwd4 = "D:/Wasim/Goldbach/revision/control/"
        nwd5 = "D:/Wasim/streu/control/"
    else
        # Modify these paths for your WSL setup
        nwd = "/mnt/d/Wasim/regio/control/"
        nwd2 = "/mnt/d/temp/saale/control/"
        nwd3 = "/mnt/d/Wasim/Tanalys/DEM/brend_fab/control/"
        nwd4 = "/mnt/d/Wasim/Goldbach/revision/control/"
        nwd5 = "/mnt/d/Wasim/streu/control/"
    end

    paths = [nwd, nwd2, nwd3, nwd4, nwd5]

    for path in paths
        cd(path)
        println("Searching in directory: $path")
        files = filter(file -> endswith(file, ".ctl"), readdir())
        for file in files
            open(file) do f
                counter = 0
                for line in eachline(f)
                    counter += 1
                    if Base.contains(line, snippet)
                        printstyled("$counter: $path", color=:light_red)
                        printstyled("$file:\t", color=:light_magenta, underline=true, blink=false, bold=true)
                        printstyled("$line\n", color=:green, bold=true)
                    end
                end
            end
        end
    end

    cd(owd)
end

# function rglob(prefix::AbstractString)
#     rootdir=pwd();
#     results = []
#     if (any(x->isdir(x),readdir()))
#         for (looproot, dirs, filenames) in walkdir(rootdir)
#             for filename in filenames
#                 #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                 if (occursin(Regex(prefix,"i"),filename))
#                     push!(results, joinpath(looproot, filename)) 
#                 end
#             end
#         end
#     else
#         printstyled("no dirs in $rootdir !\n",color=:light_red)
#         for filename in (filter(x->isfile(x),readdir(;join=false)))
#             #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#             if (occursin(Regex(prefix,"i"),filename))
#                 push!(results, filename) 
#             end
#         end
#     end
#     return results
# end

# function rglob(prefix::Regex)
#     rootdir="."
#     results = []
#     for (looproot, dirs, filenames) in walkdir(rootdir)
#         for filename in filenames
#             if (occursin(prefix,filename))
#                 push!(results, joinpath(looproot, filename)) 
#             end
#         end
#     end
#     return results
# end
"""
function rglob(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())
"""
function rglob(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())
    #check it rootdir is a dir
    rootdir = isdir(rootdir) ? rootdir : throw(ArgumentError("rootdir is not a directory!"))    
    results = []
     # checks if all items in rootdir are directories
    if (any(x->isdir(joinpath(rootdir, x)), readdir(rootdir)))
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if pattern isa AbstractString
                    if (occursin(Regex("$(pattern)", "i"), filename))
                        push!(results, joinpath(looproot, filename)) 
                    end
                else
                    if (occursin(pattern, filename))
                        push!(results, joinpath(looproot, filename)) 
                    end
                end
            end
        end
    else
        printstyled("no dirs in $rootdir !\n", color=:light_red)
        for filename in (filter(x->isfile(x), readdir(rootdir; join=true)))
            if pattern isa AbstractString
                if (occursin(Regex("$(pattern)", "i"), filename))
                    push!(results, filename) 
                end
            else
                if (occursin(pattern, filename))
                    push!(results, filename) 
                end
            end
        end
    end
    return results
end

"""
go dir up
"""
function cdu()
    dirname(pwd())|>cd
    pwd()|>println
end


visited_dirs = [pwd()]

"""
go to last dir
"""
function cdb()
    if length(visited_dirs) == 1
        println("No previous directory visited.")
    else
        pop!(visited_dirs)
        dir = visited_dirs[end]
        cd(dir)
        println("Went back to directory: ", dir)
    end
end

function visitdir(dir)
    push!(visited_dirs, dir)
    cd(dir)
    println("Visited directory: ", dir)
end

# function cdb()
#     if Sys.islinux()
#         opwd = `cd (OLDPWD)`
#         run(opwd)
#         pwd()|>println
#     else
#         try
#             #opwd=(`powershell -noprofile Set-Location $env:OldPWD`) 
#             #Cmd(["pwsh","-noprofile", "-Command", "Set-Location", "$env:OldPWD"])
#             # run(opwd)
#             # pwd()|>println
#             prevdir = pwd()
#             run(`cmd /c cd /d ..`)
#             println("Previous directory: ", prevdir)
#             println("Current directory: ", pwd())
#         catch e
#             println(e)
#             @error "pwrsh errored..."  
#         end
#     end
# end

# prev_location = ENV["OLDPWD"]
# run(`cd $prev_location`)
# # List PowerShell's Environmental Variables
# Get-Childitem -Path Env:* | Sort-Object Name

#like jdd to vector of strings.
function fdd(;cwd=pwd())
    dirs = readdir(cwd)
    
    if length(dirs) == 0 
        println("$cwd is empty!")
        return
    end
    
    if filter(x -> (isdir(x)),dirs) == []
        bn = basename(cwd)
        @info "no dirs in $bn !"
        dd()
        return
    end
    
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
	end
    end
    return(s)
end


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

"""
topdown dirsize
du(;cwd=pwd())
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

"""
greps from current dir iRegex
"""
function glob(x::AbstractString)
    filter(file -> occursin(Regex(x,"i"),file), readdir())
end

"""
greps from current dir Regex
"""
function glob(x::Regex)
    filter(file -> occursin(x,file), readdir())
end

"""
fdi(;xm::Regex=r"*")
lists dirs if isdir(dir) & occursin(xm,dir)
"""
function fdi(;cwd::AbstractString=pwd(),xm::Regex=r"")
        dirs = (length(cwd)>1) ? readdir(cwd) : readdir(pwd())
        for dir in dirs 
            if isdir(dir) & occursin(xm,dir)
                @printf("%-8s\t|","$dir")
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
    osize = 0
    for (root, dirs, files) in walkdir(cwd)
     for file in files
         osize += stat(joinpath(root, file)).size
     end
    end 
    @printf("%-40s %15.3f GB\n","$(cwd):",osize/1024^3);
end 

"""
prints sorted sizes of dirs and returns a vector of the dirs
"""
function xct(ext::AbstractString)
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

"""
gets sorted DF by size recursivley
"""
function fz()
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



function xct()
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

"""
gets sorted DF by size recursivley
"""
function fz()
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

function homrc()
    if platform == "unix"
        cd("/mnt/d/Wasim/regio/out/rc200");
    else
        cd("D:/Wasim/regio/out/rc200");
    end
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

function homqm()
    pt="D:/remo/qm/"
    cd(pt)
    println("you are here: ",pwd())
end

"""
list_files_sorted_by_last_change
formerly lat()
"""
function lf()
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
macro glb(s) glob(s);end
macro gl(s) glob(s)|>first;end
macro hd(df) df[1:4,:];end
#fastplot
macro fp(s) dfp(Regex(s));end
macro flog(s) dfl(Regex(s));end
macro ncrm() ncrem=src_path*"/ncremover.jl";include(ncrem);end
macro rasrm() remer=src_path*"/raster_remover.jl";include(remer);end
macro nco(s) first(nconly(s));end
macro dfo(s) first(dfonly(s));end
macro wajs() pt=src_path*"/wajs.jl";include(pt);end
macro jljs() pt=src_path*"/wajs.jl";include(pt);end
macro cmk() pt=src_path*"/cairomakie.jl";include(pt);end
macro rcall() pt=src_path*"/rcall.jl";include(pt);end
macro rgof() pt=src_path*"/RCall_gof.jl";include(pt);end
macro pj() pt=src_path*"/pyjl.jl";include(pt);end
macro pyjl() pt=src_path*"/pyjl.jl";include(pt);end

# using TOML
# #TOML is short for Tom’s Obvious Minimal Language and is a configuration file format 
# #that should be “easy to parse into data structures in a wide variety of languages”
# macro toml_str(s::String)
#     TOML.parse(s)::Dict{String, <:Any}
# end

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

function regand(v::Vector{String},x1::AbstractString,y1::AbstractString)
    needle=join([x1,y1],"+.*");
    z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))] 
    return(z)
end

function regand(v::Vector{String},xv::Tuple{String, String})
    needle=join([xv[1],xv[2]],"+.*");
    z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
    return(z)
end

function regand(v::Vector{String},xv::Vector{Symbol})
    needle=join(xv,"+.*");
    z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
    return(z)
end

function regand(v::Vector{String},xv::Vector{String})
    needle=join(xv,"+.*");
    z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
    return(z)
end


"""
join([x1,y1],"+.*")
r"this+.*that"
"""
function regand(x1::AbstractString,y1::AbstractString)
    needle=join([x1,y1],"+.*");
    z = Regex(needle,"i")
    return(z)
end

function regand(v::Vector{String},xv::Tuple{String, String})
    needle=join([xv[1],xv[2]],"+.*");
    z = v[(broadcast(x->occursin(Regex(needle,"i"),x),v))]
return(z)
end

"""
here you can put any regex to filter the Vector
like regand(getnames(dfs),r"tem")
"""
function regand(v::Vector{Any},xv::Regex)
    z = v[(broadcast(x->occursin(xv,x),v))] 
    return(z)
end

"""
dir_path::String, match::String;file_ending="ctl"
"""
function ctlg(dir_path::String, match::String;file_ending="ctl")
    for file in readdir(dir_path)
        #if ( occursin(file_ending, file) &&
        if ( endswith(file,file_ending) &&
            !occursin(r"tar$|tex$|pl$|sh$|csv|html|xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg", file))
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

"""
list_files_sorted_by_last_change
"""
function latx()
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

ll=latx

"""
grabs methods
asin|>getm  
?asin
@code_llvm readf|>getm|>first 
"""
function getm(s::Any)
    methods(s);
end

function calculate_folder_size(directory)
    size = 0
    count = 0
    for (root, dirs, files) in walkdir(directory)
        for file in files
            size += stat(joinpath(root, file)).size
            count += 1
        end
    end
    return size, count
end

function print_folder_size(directory, size, count)
    size_gb = round(size / 1024^3, digits=3)
    printstyled(rpad("$directory: $size_gb GB", 40), color=:green)
    printstyled(lpad("($count files)", 20), "\n", color=:green)
end

function sdf()

    function calculate_folder_size(directory)
        size = 0
        count = 0
        for (root, dirs, files) in walkdir(directory)
            for file in files
                size += stat(joinpath(root, file)).size
                count += 1
            end
        end
        return size, count
    end

    function print_folder_size(directory, size, count)
        size_gb = round(size / 1024^3, digits=3)
        printstyled(rpad("$directory: $size_gb GB", 40), color=:green)
        printstyled(lpad("($count files)", 20), "\n", color=:green)
    end

    printstyled("folder sizes on Julia...\n", color=:red)

    cwd = pwd()
    dirs = readdir(cwd)

    rs, rcnt = calculate_folder_size(cwd)

    print_folder_size(cwd,rs,rcnt)

    n = repeat(" - -", 10)
    printstyled(n * " subfolders of " * basename(cwd) * n, "\n", color=:yellow)

    for dir in dirs
        if isdir(dir)
            size, count = calculate_folder_size(joinpath(cwd, dir))
            print_folder_size(dir, size, count)
        end
    end
end

"""
Read wasim ts with DelimitedFiles.readdlm, skipto line 3 
no header column
"""
function wread(x::String;skip=3)
    df = DelimitedFiles.readdlm(x, '\t', Float64, '\n';
        header=false,skipstart=skip)
    df = DataFrame(df,:auto)
    for i in 5:size(df,2)
        df[!,i]=replace(df[!,i],-9999.0 => missing)
    end 
    for i in 5:size(df,2)
        replace!(df[!,i],-9999.0 => missing)
    end
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df=df[:,Not(1:4)]
    metadata!(df, "filename", x, style=:note);
end

"""
recursively removes duplicates 
uses SHA
"""
function rmdub(;directory=pwd())
    
    # Create a dictionary to store file hashes as keys and file paths as values
    #hash_dict = Dict{String, String}()
    hash_dict = Dict{Vector{UInt8}, String}()

    # Get a list of all files in the directory
    #files = readdir(directory)
    for (root, dirs, files) in walkdir(directory)
        for file in files
            filepath = joinpath(directory, file)
            if isfile(filepath)
                # Calculate the SHA-256 hash of the file's contents
                #hash = string(sha256(open(filepath, "r")))
                io = open(filepath, "r")
                #filehash = string(sha256(io))
                #filehash = string(hash(io))
                filehash = sha256(io)
                #println(filehash)
                close(io)

                # If the hash is not already in the dictionary, add it
                if !haskey(hash_dict, filehash)
                    hash_dict[filehash] = filepath
                else
                    # If a file with the same hash is found, delete it
                    println("Deleting duplicate file: $filepath")
                    rm(filepath)
                end
            end
        end
    end
end

function op()
    #open(`powershell -noprofile explorer . `,"w",stdout) Set-Location $env:OldPWD
    run(`cmd.exe /c start .`)
end 

function cntcolv(x::String)
    # Get a list of all files in the directory
    #x = "."
    files = filter(file -> (occursin(Regex(x, "i"), file) & 
                        (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file))
                        ), readdir())

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

"""
mkdir("route-bak")
cpinto(glob("so_inf"), "route-bak")
rglob("so_inf")
force=true will first remove an existing dst.
"""
function cpinto(src::Vector{String}, dst::AbstractString;force=false)
    map(x->cp(x,"$dst/$x";force=force),src)
end

function getmoduleorder(file1::AbstractString)
    # Run the awk command and capture the STDOUT
    result = read(`awk '$0 ~ /^[[]/{c=1} c&&c--' $file1`, String)

    # Split the result into lines and strip \r and \t characters
    lines = replace(split(result, "\n"), r"[\r\t]" => "")

    # Extract the names inside square brackets using regular expressions
    pattern = r"\[(.*?)\]"
    names_inside_brackets = [match(pattern, line) === nothing ? "" : match(pattern, line).captures[1] for line in lines]

    # Create a DataFrame with the extracted names inside the square brackets
    df = DataFrame(Line = names_inside_brackets)

    return df
end

function ctg(match::String ; dir=".", file_ending=".ctl")
    for file in readdir(dir)
        if occursin(file_ending, file)
            fx = joinpath(dir, file)
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

function kgewrite(;output_file="kge-table.txt")
    path = pwd()
    files = glob(r"_output.txt|_outputjl") #non-recurse
    
    # Create an empty vector to store the matched lines
    matched_lines = Vector{String}()
    
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"KGE", output)
        if !isempty(match)
            fn = first(split(file, "_qout"))
            for line in sort(match, by = x -> parse(Float64, split(x)[end]); rev = true)
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  # remove inner whitespaces
                push!(matched_lines, "$fn: $line")  # collect the matched line
            end
        end
    end
    
    #output_file = "matched_results.csv"
    writedlm(output_file, matched_lines, '\t')
    println("Matched results saved to $output_file")
end

function rsqgrep()
    path = pwd()
    files = glob(r"_output.txt|_outputjl")
    @printf("Searching for R² values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"R2.*[0-9].[3-9]", output)
        if !isempty(match)
            fn = first(split(file, "_qout"))
            for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                line = strip(line)
                line = join(split(line), " ")
                printstyled(rpad("$fn:", 30), lpad("$line\n", 10), color = :green)
            end
        end
    end
end

#rsqgrep()
function rsq(x::AbstractVector, y::AbstractVector)
    cor(x, y)^2
end

function kgegrep()
    path = pwd()
    files = glob(r"_output.txt|_outputjl") #non-recurse
    #@printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        #match = Grep.grep(r"KGE.*[0-9].[3-9]", output)
        match = Grep.grep(r"KGE", output)
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

function nsegrep()
    path = pwd()
    files = glob(r"_output.txt|_outputjl") #non-recurse
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"mNSE", output)
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

function kgegrepr()
    path = pwd()
    files = rglob(r"_output.txt|_outputjl") #recurse
    @printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file,'\t', String)
        match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
        if !isempty(match)
            fn = first(split(file,"_qout"))
            for line in match
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
            end
        end
    end
end

function nsegrepr()
    path = pwd()
    files = rglob(r"_output.txt|_outputjl") #recurse
    @printf("Searching for NSE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file,'\t', String)
        match = Grep.grep(r"NSE.*[0-9].[3-9]",output)
        if !isempty(match)
            fn = first(split(file,"_qout"))
            for line in match
                line = strip(line) 
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
            end
        end
    end
end

function nconly(x1::AbstractString)
    v::Vector{String} = readdir();
    v = v[broadcast(x->endswith(x,"nc"),v)];
    z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
    return(z)
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
    readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
    return(z)
end

function dfon(x1::AbstractString)
    z = filter(file -> occursin(Regex(x1,"i"),file), 
    readdir()[broadcast(x->!endswith(x,"nc"),readdir())]);
    return(z)
end

# x=nothing
# typeof(x)
# isa(x,Any) #true
#Nothing,

function nconly(rx::Union{AbstractString,Regex})
    z = filter(z -> endswith(z,".nc"),readdir()[map(x->occursin(rx,x),readdir())])
    # v = readdir();
    # z = v[map(x->occursin(rx,x),v)]
    # z = z[map(x->endswith(x,"nc"),z)]
    return(z)
end

function nconly()
    z = filter(z -> endswith(z,".nc"),readdir())
    return(z)
end


"""
tff2
"""
function jlt(x::Vector{String})
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
                    #println("filename: $filename")
                    printstyled("filename: $filename\n",color=:yellow)
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

function tf(x::Union{AbstractString,Regex})
    jlt(dfonly(x))
end

"""
water-balance.jl
"""
function waba()
    #wpth="C:/Users/Public/Documents/Python_Scripts/julia/water-balance.jl"
    wpth=src_path*"/water-balance.jl"
    include(wpth)
    yd=waread("waba-input.wa")|>yrsum
    @info "waba done !"
end

"""
filters internal WaSiM stats of routed discharge files
works recursively
"""
function qgk(;rootdir=".", prefix="qgk")
    files = []
    for (looproot, dirs, filenames) in walkdir(rootdir)
        for filename in filenames
            if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
                push!(files, joinpath(looproot, filename)) 
            end
        end
    end
    
    for z in files
        println(raw"file:	",basename(z),"...")
        m = filter(line -> occursin(r"^[LIN. R]|^[LOG. R]|^CO",line), readlines(open(z)))
        for l in m
            x = replace(l, r"\s+" => "\t")
            x = replace(x, ".\t" => " ")
            println(x)
        end
    end
    return nothing
end

function readbetween(io::IO, start::String, stop::String)
    output = Vector{String}()
    while !eof(io)
        line = readline(io)
        if Base.contains(line, start)
            push!(output, line)
            while !eof(io)
                line = readline(io)
                if Base.contains(line, stop)
                    skip(io, length(line))
                    break
                end
                push!(output, line)
            end
            break
        end
    end
    return output
end

function readbetween(io::IO, start::Regex, stop::Regex)
    output = Vector{String}()
    while !eof(io)
        line = readline(io)
        if occursin(start,line)
            push!(output, line)
            while !eof(io)
                line = readline(io)
                if occursin(stop,line)
                    skip(io, length(line))
                    break
                end
                push!(output, line)
            end
            break
        end
    end
    return output
end

"""
prints out the function definition
"""
function zp(func::Any)
    #pt = joinpath(@__DIR__,"func-win.jl")
    pt = src_path*"/func-win.jl"
    _str = "$(func)"
    readbetween(open(pt),Regex(_str),r"^\s*function")
end

macro wajs() pt="C:\\Users\\Public\\Documents\\Python_Scripts\\julia\\wajs.jl";include(pt);end
macro bash_str(s) open(`bash`,"w",stdout) do io; print(io, s); end;end

#@bash_str "which python" #redirects to wsl bash
# @bash_str "python -V"
# @bash_str "python -c 'import sys; print(sys.version_info)'"
# @bash_str @bash_str "python -c 'import rasterio; print(rasterio.__version__)'"

#@bash_str "python -c 'import rasterio; print(rasterio.__version__)'"
#@cmd_str "python -c 'import rasterio; print(rasterio.__version__)\n\n'"

#pwrs""" which python """
#pwrs""" pwd """

"""
this works in wsl too
run(`cmd.exe /c start .`)
"""
function op()
    #pwrs""" explorer . """
    #open(`powershell -noprofile explorer . `,"w",stdout)
    #open(`cmd.exe /c start . `,"w",stdout)
    #this wrks in wsl, too
    run(`cmd.exe /c start .`)
end 

macro pwp_str(s) open(`powershell`,"w",stdout) do io; print(io, s); end;end
macro cmd_str(s) open(`cmd \c`,"w",stdout) do io; print(io, s); end;end
##nope, bad idea
# run(`cmd \c pwd";" exit`)
function rsqgrep()
    path = pwd()
    files = glob(r"_output.txt|_outputjl")
    @printf("Searching for R² values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        match = Grep.grep(r"R2.*[0-9].[3-9]", output)
        if !isempty(match)
            fn = first(split(file, "_qout"))
            for line in sort(match, by = x -> parse(Float64, split(x)[end]);rev=true)
                line = strip(line)
                line = join(split(line), " ")
                printstyled(rpad("$fn:", 30), lpad("$line\n", 10), color = :green)
            end
        end
    end
end
#rsqgrep()
function rsq(x::AbstractVector, y::AbstractVector)
    cor(x, y)^2
end
function kgegrep()
    path = pwd()
    files = glob(r"_output.txt|_outputjl") #non-recurse
    #@printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file, '\t', String)
        #match = Grep.grep(r"KGE.*[0-9].[3-9]", output)
        match = Grep.grep(r"KGE", output)
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

function grep_KGE(path::AbstractString)
    #@printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), readdir(path))
        #output = read(file, String)
        output = DelimitedFiles.readdlm(file,'\t', String)
        #match = occursin(r"KGE.*[0-9].[3-9]", output)
        match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
        if !isempty(match)
            #@printf("%s: %s\n", file,match)
            #@printf("%s:", first(split(file,"_qout")))
            fn = first(split(file,"_qout"))
            for line in match
                #@printf("\t%s\n", line)
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
            end
        end
    end
end

function kgegrepr()
    path = pwd()
    files = rglob(r"_output.txt|_outputjl") #recurse
    @printf("Searching for KGE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file,'\t', String)
        match = Grep.grep(r"KGE.*[0-9].[3-9]",output)
        if !isempty(match)
            fn = first(split(file,"_qout"))
            for line in match
                line = strip(line)  # remove leading and trailing whitespace
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
            end
        end
    end
end

function nsegrepr()
    path = pwd()
    files = rglob(r"_output.txt|_outputjl") #recurse
    @printf("Searching for NSE values > 0.3 in files matching pattern %s\n", path)
    for file in filter(file -> endswith(file, "_output.txt"), files)
        output = DelimitedFiles.readdlm(file,'\t', String)
        match = Grep.grep(r"NSE.*[0-9].[3-9]",output)
        if !isempty(match)
            fn = first(split(file,"_qout"))
            for line in match
                line = strip(line) 
                line = join(split(line), " ")  ##remove inner whitespaces
                printstyled(rpad("$fn:",30),lpad("$line\n",10),color=:green)
            end
        end
    end
end

"""
names and size in MB via readdir
"""
function pfix()
    #files = filter(x -> occursin(Regex(m,"i"), x), readdir())
    files = readdir()
    fzs = [(file, filesize(file) / 2^20) for file in files]
    tot = round(sum(map(x->x[2],fzs));digits=3)
    #d=Dict(fzs) 
    #sort(collect(d), by = x -> x[2];rev=true)
    sv = sort(collect(fzs), by = x -> x[2];rev=true)
    sv = first(sv,8)
    for i in sv
        printstyled(rpad(i[1],40),
        lpad(round(i[2];digits=4),10)," MB\n",
        color=:yellow)
    end
    #sv = first(Dict(sv),4)
    #printstyled("biggest files :\n $sv \n",color=:yellow)
    lf = length(files)
    printstyled("$lf files ...\n Total Size: $tot MB\n",color=:green)
end
    
"""
varinfo(Main,r"ds")
"""
function vars()
    InteractiveUtils.varinfo(;sortby=:size,minsize=1)
end

function rtreeo()
    dir = pwd()
    dirs = readdir(dir)
    routes::Vector{String} = []
    for directory in dirs
        if isdir("$dir/" * directory)
            push!(routes, "$dir/$directory")
        else
            if ~(directory in routes)
                push!(routes, "$routes/$directory")
                #newread = dir * "/$directory"
                #push!(newread, "$dir/$directory")
                #[push!(routes, r) for r in newrs]
            end
        end
    end
    routes
end

function treeo(;root = pwd(), indent::String = "    ")
    println(root)
    for (root, dirs, files) in walkdir(root)
        # for file in files
        #     println(indent, "├── ", file)
        # end
        for d in dirs
            dp = joinpath(root, d)
            if Base.contains(dp, "\\")
                ln = splitpath(dp)|>length
                dp = replace(dp, root=>"","\\" => "──")
                println(indent,"├──",repeat("────",ln-1), dp)
            end
            #println(indent, "├── ", dp)
        end
    end
end

"""
5 biggest folders.
"""
function bigf()

    printstyled("5 biggest folders...\n", color=:red)

    cwd = pwd()
    dirs = filter(isdir, readdir(; sort=false))
    
    rs, rcnt = [], []

    for i in dirs
        ts, tcnt = calculate_folder_size(i)
        push!(rs, round(ts / 1024^2, digits=3))  # Convert size to MB and round to 3 digits
        push!(rcnt, tcnt)
    end

    # Check if the matrix is empty
    if isempty(dirs)
        println("No directories found.")
        return
    end
    
    # Create an array of directories, sizes in MB, and file counts
    data = [(dir, rs[i], rcnt[i]) for (i, dir) in enumerate(dirs)]

    # Sort the array by size (second element) in descending order
    sorted_data = sort(data, by = x -> x[2], rev = true)

    if isempty(sorted_data)
        println("No directories with size information found.")
        return
    end

    # Extract the largest directories and display them in green
    largest_dirs = sorted_data[1:5]
    
    # Create a formatted string for printing
    dirs_to_display = join([string(dir[1], "\t\t Size: ", dir[2], " MB \t Files: ", dir[3], "\t") for dir in largest_dirs], "\n")
    #dirs_to_display = join([string(dir[1], " (Size: ", dir[2], " MB, Files: ", dir[3], ")") for dir in largest_dirs], "\n")

    printstyled("$dirs_to_display\n--> ", color=:green)
end

"""
import InteractiveUtils.clipboard
wslpath()|>clipboard
cb(wslpath())
"""
function wslpath()
    # Run the `wslpath` command to convert the current directory to a WSL path
    #wsl_cmd = `wsl wslpath -a $(pwd)`
    wsl_cmd = `wsl wslpath -a .`
    wsl_path = readchomp(pipeline(wsl_cmd))
    # Return the WSL path
    return wsl_path
end

"""
import InteractiveUtils.clipboard
wslpath()|>clipboard
cb(wslpath())
"""
function pwc()
    if platform == "windows"
        pt = wslpath()
        println("$pt in clipboard...")
        pt |> clipboard
    else
        pt = pwd()
        pt |> clipboard
        println("$pt in clipoard...")
    end
end

"""
return clipboard(x)
"""
function cb(x)
    return clipboard(x)
end

function oplat()
    files = filter(isfile, readdir())
    sorted_files = sort(files, by = mtime, rev = true)
    if !isempty(sorted_files)
        lat = sorted_files[1]
        println("opening: $lat ...")
        if platform == "osx"
            run(`open $lat`)
        else
            try
                run(`cmd /c start $lat`)    
            catch
                @error "could not open $lat via cmd.exe"
            end
        end
    end
end   

"""
removes latest file
"""
function rmlat()
    files = filter(isfile, readdir(;sort=false))
    sorted_files = sort(files, by = mtime, rev = true)
    if !isempty(sorted_files)
        lat = sorted_files[1]
        println("This deletes the latest created file, i.e: ", lat)
        print("continue? (y/n): ")
        reply = readline(stdin)
        println(reply)
        if lowercase(string.(reply)) == "y"
            rm(lat, force=true)
            println("Deleted: ", lat)
        else
            @error "abort...."
        end
    end
end

function lat()
    files = filter(isfile, readdir(;sort=false))
    sf = sort(files, by = mtime, rev = true)
    lat = sf[1]
    if length(sf) < 6
        sf = join(sf,"\n")
    else
        sf = join(sf[1:6],"\n")
    end
    printstyled("$sf\n--> ",color=:green)
    return lat
end

function tree_helper(root, indent)
    
    entries = filter(isdir, readdir(root))
    
    # #das geht auch, aber too much
    # entries = []
    # for (root, dirs, files) in walkdir(root)
    #     for dir in dirs
    #         #dstr = (joinpath(root, dir))
    #         dstr = dir
    #         push!(entries, dstr)
    #     end
    # end
    
    for (i, entry) in enumerate(entries)
        entry_path = joinpath(root, entry)  
        is_last = i == length(entries)
        entry_path = replace(entry_path, "/"=>"\\") #for win.
        
        entry_path = replace(entry_path, "\\mnt"=>"mnt") #
        # Calculate the depth of the current directory
        #depth = count(x -> x == '/', entry_path) - count(x -> x == '/', root)
        depth = count(x -> x == '\\', entry_path) - count(x -> x == '\\', root)
        
        # Print indentation
        if depth > 0
            print(indent)
            print(is_last ? "└── " : "├── ")
        end
        
        # Print the current entry name
        println(entry)
        
        # Recursively process subdirectories
        tree_helper(entry_path, indent * (is_last ? "    " : "│   "))
    end
end

function tree(;root = pwd())
    tree_helper(root, "")
end

function looks_like_number(str::AbstractString)
    try
        parse(Float64, str)
        return true
    catch
        return false
    end
end

function pyread(x::Union{String,Regex})
    """
    pyreader, reads all as stings, conversion later.
    """
    if x isa Regex
        x = dfonly(x)|>first
    end
    pd = pyimport("pandas")
    df = pd.read_table(x, 
        engine="c",
        verbose=true,
        low_memory=false,
        header=0,
        skipinitialspace=true,
        dtype="str",              #new!
        na_values=[-9999]
        )
    col_names = df.columns  # Get the column names from the Python DataFrame
    col_names = convert(Array, col_names)
    col_arrays = [convert(Array, df[col]) for col in col_names]
    filtered_rows = broadcast(x->looks_like_number(x),col_arrays[1])
    
    df = DataFrame(col_arrays, :auto)
    
    df = df[filtered_rows, :]
    
    rename!(df, col_names)
    if "YY" ∉ names(df)
        println("Column 'YY' not found in the CSV file.")
        return nothing
    end

    function parse_date(row)
        try
            year = parse(Int, row[1])
            month = parse(Int, row[2])
            day = parse(Int, row[3])
            return Date(year, month, day)
        catch
            return missing  # Return missing if parsing fails
        end
    end
    
    df.date = map(parse_date, eachrow(df))
    
    df = select(df, Not(1:4))
    
    for x in names(df)
        if looks_like_number(x)
            newname=replace(x,"$x"=>"C$x", count=1)
            rename!(df,Dict(x=>newname))
        end
    end

    #map(y->typeof(y),eachcol(df))

    # Iterate over column names
    for colname in names(df)
        # Check if the column type is not Date
        if eltype(df[!, colname]) != Date
            df[!, colname] .= tryparse.(Float64, df[!, colname])
        end
    end

    DataFrames.metadata!(df, "filename", x, style=:note)
    return df 
end

function read_df(s::Union{String,Regex})
    """
    reads a DataFrame from a file w dlm and tryparse subsetting
    """
    if s isa Regex
        s = dfonly(s)|>first
    end
    data, colnames = DelimitedFiles.readdlm(s, '\t', String, '\n', header=true)
    df = DataFrame(Matrix{Any}(data), :auto)
    rename!(df,Symbol.(colnames)|>vec)
    df = df[map(x->!isnothing(x),tryparse.(Int,df[!,1])),:]
    for i in 1:size(df, 2)
        df[!, i] = map(x -> parse(Float64, x), df[!, i])
    end
    for i in 5:size(df,2)
        df[!,i]=replace(df[!,i],-9999.0 => missing)
    end 
    for i in 5:size(df,2)
        replace!(df[!,i],-9999.0 => missing)
    end
    
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df = df[:,Not(1:4)]
    metadata!(df, "filename", s, style=:note);
    return df
end

"""
Read wasim ts with DelimitedFiles.readdlm, skipto line 3 
no header column
"""
function wread(x::String;skip=3)
    df = DelimitedFiles.readdlm(x, '\t', Float64, '\n';
        header=false,skipstart=skip)
    df = DataFrame(df,:auto)
    for i in 5:size(df,2)
        df[!,i]=replace(df[!,i],-9999.0 => missing)
    end 
    for i in 5:size(df,2)
        replace!(df[!,i],-9999.0 => missing)
    end
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df=df[:,Not(1:4)]
    metadata!(df, "filename", x, style=:note);
end

function wread(x::Regex;skip=3)
    """
    Read wasim ts with DelimitedFiles.readdlm, skipto line 3 
    no header column
    """
    rgx = glob(x)|>first
    println("loading $rgx ...")
    df = DelimitedFiles.readdlm(rgx, '\t', Float64, '\n';
        header=false,skipstart=skip)
    df = DataFrame(df,:auto)
    for i in 5:size(df,2)
        df[!,i]=replace(df[!,i],-9999.0 => missing)
    end 
    for i in 5:size(df,2)
        replace!(df[!,i],-9999.0 => missing)
    end
    for i in 1:3
        df[!,i]=map(x ->Int(x),df[!,i])
    end
    #and parse dates...
    df.date = Date.(string.(df[!,1],"-",df[!,2],"-",df[!,3]),"yyyy-mm-dd");
    df=df[:,Not(1:4)]
    metadata!(df, "filename", x, style=:note);
end

function wintree()
    #run(`cmd /c tree /f`)
    run(`cmd /c tree`)
end

"""
cd into dir of file
"""
function cdof(x::Union{String,DataFrame})
    if x isa DataFrame
        try
            d = collect(DataFrames.metadata(x))[1][2]
            cd(dirname(d))
        catch
            @error "no basename in present!"
            return nothing
        end
    else
        try
            cd(dirname(x))
        catch
            @error "cd errored!"
            return nothing
        end
    end
    d=pwd()
    println("current dir: $d")
end

#k=raw"C:/Users/Public/Documents/Python_Scripts/julia/smallfuncs.jl"
println("script loc is $(src_path)/smallfuncs.jl")

#homedir()|>cd
#home() ##necessary for ssup to work

# @vv "surface"
# x=raw"C:/Users/chs72fw/Documents/Promotionsstudium/Dropbox/brendfab/m3/Soil_Temperature_Stack.nc"
# r = Raster(x)
# r = readras(x)
# plot(r[t=4:5])
# surface(r[t=5],camera = (0,-75),legend=false)

function pyplot_df(df::DataFrame;log=false)
    x = df.date
    ln = (filter(x -> !occursin(r"date|month|year", x), names(df)))

       
    for col in ln
        y = df[!, Symbol(col)]
        PyPlot.plot(x, y, label=col)
    end

    if log
        PyPlot.yscale("log")
    end

    PyPlot.xlabel("Date")
    PyPlot.ylabel("")
    PyPlot.legend()
    ti = only(values(DataFrames.metadata(df)))
    PyPlot.title(ti)
    PyPlot.grid(true)
end

"""
Returns the total size of all files in the current directory non recursively.
"""
function fsz(rootdir::String=".";rec=false)
    total_size = 0
    files = readdir(rootdir)  # Get a list of files in the current directory
    
    for file in files
        filepath = joinpath(rootdir, file)  # Get the full path of the file
        if isfile(filepath)
            size = stat(filepath).size  # Get the file size
            total_size += size
        end
    end
    
    nr=length(files)

    total_size_mb = total_size / (1024 * 1024)  # Convert size to megabytes
    total_size_mb = round(total_size_mb, digits=2)  # Round to 2 decimal places
    printstyled("Total size of $nr files in $(pwd()): $total_size_mb MB\n", color=:green)
    
    if rec
        dirs = readdir(rootdir)
        for dir in dirs
            if isdir(dir)
                cd(dir)
                fsz()
                cd("..")
            end
        end
    end
    #return total_size_mb
end

function fsize()
    """
    ALL names and size in MB via readdir
    """
    files = filter(x -> isfile(x), readdir())
    fzs = [(file, filesize(file) / 2^20) for file in files]
    tot = round(sum(map(x->x[2],fzs));digits=3)
    fzs = sort(fzs, by = x -> x[2], rev = true)
    odf = rename(DataFrame(fzs),["file","size"])
    DataFrames.metadata!(odf, "Total Size", tot, style=:note)
    printstyled("Total Size: $tot MB\n",color=:green)
    return(odf)
end


"""
Read the text file, preserve line 1 as header column
Instead of using CSV.read, we use CSV.File to create a lazy representation of the file.
This avoids reading the entire file into memory at once, 
which can be more memory-efficient for large datasets.
"""
function waread2(x::String)
    ms = ["-9999", "lin", "log", "--"]
    df = CSV.File(x; delim="\t", header=1, normalizenames=true, missingstring=ms, types=Float64) |> DataFrame
    dropmissing!(df,1)
    dt2 = [Date(Int(row[1]), Int(row[2]), Int(row[3])) for row in eachrow(df)]
    select!(df, Not(1:4))
    df.date = dt2
    metadata!(df, "filename", x, style=:note)
    return df
end

"""
basic tsv reader, takes arguments from CSV.File
df = CSV.File(x;kw...)|>DataFrame|>z->dropmissing(z,1)
"""
function tsread(x::Union{String,Regex};kw...)
    if x isa String
        printstyled("reading $x\n",color=:light_red)
    else x isa Regex
        x = first(dfonly(x))
        printstyled("reading $x\n",color=:light_red)
    end
    #ms = ["-9999","-9999.0","lin", "log", "--"]
    #missingstring=ms,
    #df = CSV.read(x,DataFrame;kw...)
    df = CSV.File(x;kw...)|>DataFrame|>z->dropmissing(z,1)
    DataFrames.metadata!(df, "filename", x, style=:note)
    for x in names(df)
        if startswith(x,"_")
            newname=replace(x,"_"=>"C", count=1)
            rename!(df,Dict(x=>newname))
        end
    end
    return df 
end


"""
Fastest Reader. is also dfr.
Read the text file, preserve line 1 as header column
"""
function dfr(x::String)
    ms = ["-9999","lin","log","--"]
    df = try
     CSV.read(x, DataFrame; 
        delim="\t", 
        header=1, 
        missingstring=ms, 
        #maxwarnings = 1, 
        silencewarnings = true,
        normalizenames=true, 
        types=Float64)
    catch e
        @error("error reading $x\nexiting now\n $e")
        return nothing
    end
    
    df = dropmissing(df, 1)
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

"""
dfr = waread
"""
function dfr(x::Regex)
    """
    Read the text file, preserve line 1 as header column
    """
    x = dfonly(x)|>first
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

function rmeq()
    """
    removes empty TS; 
    use with caution!
    """
    #x = pwd()

    # files = filter(file -> (occursin(Regex(x, "i"), file) & 
    # (!occursin(r"xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file))
    # ), readdir())
    
    files = filter(file -> 
    !occursin(r"tex|pl|sh|csv|html|xml|fzt|ftz|log|ini|wq|yrly|nc|png|svg|txt", file)
    , readdir())
    
    ms = ["-9999", "lin", "log", "--"]
    for inF in files
        if isfile(inF)
            df = CSV.File(inF; delim="\t", header=1,
            silencewarnings=true, 
                normalizenames=false, 
                missingstring=ms, 
                types=Float64) |> DataFrame
            dropmissing!(df,ncol(df))
            if nrow(df)==0
                println(basename(inF)," removed!")
                rm(inF)
            end
        end
    end
end

function pew()
    try
        in = clipboard()
        wpath = replace(in, "\\" => "/")
        println("pt=$wpath")
        clipboard("pt=$wpath")
        return wpath
    catch e
        @error "smth errored $e"
    end
end

function pe()
    try
    inp = clipboard()
    wpath = replace(inp, "\\" => "/", "\"" => "")
    cmd = `wslpath -ua $wpath`
    ot = readchomp(pipeline(cmd))
    clipboard("$ot")
    println("$ot in clipboard!")
    return string(ot)
    catch e
        println("Error: $e")
        println("Failed to translate to wslpath.")
    end
end


"""
function to cb windowspath
"""
function pwin()
    if platform == "windows"
        pt = wslpath()
        println("$pt in clipboard...")
        pt |> clipboard
    elseif platform == "osx"
        pt = pwd()
        pt |> clipboard
        println("$pt in clipboard...")
    else
        wsl_cmd = `wslpath -m .`
        wsl_path = readchomp(pipeline(wsl_cmd))
        # Return the WSL path
        clipboard(wsl_path)
        println("$wsl_path in clipboard...")
        #return wsl_path
    end
end

"""
function to open in notepad++
"""
function npp(fl::String;opener="c:/Program Files (x86)/Notepad++/notepad++.exe")
    try
        if platform == "windows"
            run(`$opener $fl`)
        elseif platform == "unix" || platform == "linux" && src_path == "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
            fl = fl
            run(`"/mnt/c/Program Files (x86)/Notepad++/notepad++.exe" $fl`)
        else
            fx = `wslpath -m $fl` #translate from a WSL path to a Windows path, with '/' 
            opcmd = `wslpath -a $opener` #force result to absolute path format
            wslpt = readchomp(pipeline(opcmd))
            run(`$wslpt $fx`)

        end
    catch        
          @error "could not open $fl via notepad++
              check if notepad++ is installed in $opener"
    end
end

"""
function to open latest file in notepad++
"""
function npplat(;opener="c:/Program Files (x86)/Notepad++/notepad++.exe")
    files = filter(isfile, readdir(;sort=false))
    fl = sort(files, by = mtime, rev = true)[1]
    if endswith(fl,".nc")
        @error "$fl cannot be opened in notepad++"
        return
    end
    
    try
        if platform == "windows"
            run(`$opener $fl`)
        elseif platform == "unix" || platform == "linux" && src_path == "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
            fl = fl
            run(`"/mnt/c/Program Files (x86)/Notepad++/notepad++.exe" $fl`)
        else
            fx = `wslpath -m $fl` #translate from a WSL path to a Windows path, with '/' 
            opcmd = `wslpath -a $opener` #force result to absolute path format
            wslpt = readchomp(pipeline(opcmd))
            run(`$wslpt $fx`)

        end
    catch        
          @error "could not open $fl via notepad++
              check if notepad++ is installed in $opener"
    end
end

"""
gets pkgversion
filter(x-> x.second.name == pkg, Pkg.dependencies()) |> 
         x -> first(x)[2].version
"""
function pkgversion(pkg::String)
    #pak = Pkg.dependencies()
    # filter(x-> x.second.name == pkg, Pkg.dependencies()) |> 
    #     x -> first(x)[2].version
    #info = filter(x-> x.second.name == pkg, Pkg.dependencies())
    #
    try
        pkg = strip(pkg)
        info = filter(x-> occursin(pkg,x.second.name),
            #Regex(x.second.name,"i")), 
            Pkg.dependencies())
        nam = first(info)[2].name
        ver = first(info)[2].version
        #nam = map(x->first(x)[2].name,info)
        #ver = map(x->first(x)[2].version,info)
        println("$nam: $ver")    
    catch
        @error ("not found: $pkg")
        return
    end
    
end


"""
newer version with copy df and switched func positions
"""
function wawrite2(df::DataFrame,file::AbstractString;
        hourint::Int=24,
        misval::Union{Number,Missing}=-9999)
    dout = copy(df)
    if in("year",names(dout))
        @warn "yearcol found!"
        CSV.write(file, dout, 
        transform = (col, val) -> something(val, missing), delim="\t")  
        return
    end
    dout.YY = map(x ->year(x),dout.date)
    dout.MM = map(x ->month(x),dout.date)
    dout.DD = map(x ->day(x),dout.date)
    dout[!, "HH"] .= hourint
    dout = dout[!,Cols([:YY,:MM,:DD,:HH],Not(Cols(r"date")))]
    CSV.write(file, dout, 
    transform = (col, val) -> something(val, misval), delim="\t")  
    nothing
end



"""
writes df to file, no date conversion
"""
function writedf(df::Union{DataFrame,String},file::Union{DataFrame,String})
    if df isa String
        @info "write DataFrame to $df !"
        CSV.write(df, file,
        transform = (col, val) -> something(val, missing),
            delim="\t")
        return
    end
    @info "write DataFrame to $file !"
    CSV.write(file, df, 
    transform = (col, val) -> something(val, missing),
        delim="\t")
    nothing
end

"""
writes describe(df) to file, no date conversion
"""
function writedesc(table,file)
    CSV.write(file, describe(table), transform = (col, val) -> something(val, missing),delim="\t")  
    nothing
end



"""
read df to datetime
df = CSV.read(pt,DataFrame) 
df = df[5:end,:]
rename!(df,1 => :date)
fo = dateformat"yyyy mm dd HH MM"
df.date = [DateTime(d, fo) for d in df.date] 
"""
function dfrdt(x::Union{Regex,String})
    if x isa(Regex)
        try 
            x = first(dfonly(x))
        catch
            @error "no match for $x ! "
        end
    end
    ms = ["-9999","lin","log","--"]
    df = CSV.read(x, DataFrame; 
        delim="\t", header=1, missingstring=ms, 
        maxwarnings = 1, #silencewarnings = true,
        normalizenames=true, types=Float64)
    df = dropmissing(df, 1)
    dt2 = map(row -> Dates.DateTime(Int(row[1]), 
        Int(row[2]), Int(row[3]), 
        Int(row[4]),0,0,0),  #note that: Dates.DateTime,dateformat"yyyy mm dd HH MM SS s"
        eachrow(df))
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


"""
function to open latest file in notepad++
"""
function npplat(;opener="c:/Program Files (x86)/Notepad++/notepad++.exe")
    fl = lat()
    try
        if platform == "windows"
            run(`$opener $fl`)
        elseif platform == "unix" || platform == "linux" && src_path == "/mnt/c/Users/Public/Documents/Python_Scripts/julia"
            fl = fl
            run(`"/mnt/c/Program Files (x86)/Notepad++/notepad++.exe" $fl`)
        else
            fx = `wslpath -m $fl` #translate from a WSL path to a Windows path, with '/' 
            opcmd = `wslpath -a $opener` #force result to absolute path format
            wslpt = readchomp(pipeline(opcmd))
            run(`$wslpt $fx`)

        end
    catch        
          @error "could not open $fl via notepad++
              check if notepad++ is installed in $opener"
    end
end

"""
df to clipboard
"""
function ptb(x::DataFrame)
    strs = sprint(io -> pretty_table(io, x, header=uppercasefirst.(names(x)), backend = Val(:text)))
    clipboard(strs)
end

function count_files(path::String, level::Int64)
    n = 0  # Number of files
    s = 0  # Total size in bytes
    
    for entry in readdir(path)
        if entry == "." || entry == ".."
            continue
        end
        
        subpath = joinpath(path, entry)
        stat_result = stat(subpath)
        
        if isfile(stat_result)
            n += 1
            s += stat_result.size
        elseif isdir(stat_result)
            println("$(repeat(" ", level * 2))[$entry]")
            subn, subs = count_files(subpath, level + 1)
            n += subn
            s += subs
        end
    end
    
    return n, s
end

function jlcnt(path=pwd(), level=0)
    n, s = count_files(path, 0)    
    printstyled("Directory: $path\n",color=:yellow)
    printstyled("Number of files: $n\n",color=:yellow)
    printstyled("Total size: $(round(s / (1024 * 1024), digits=2)) MB\n",color=:yellow)
end

"""
wslpath to windows
"""    
function towin(file_path::Union{String, Nothing})
    if isnothing(file_path)
        lw = split(pwd(), '/')[3]
        rst = join(split(pwd(), '/')[4:end], '/')
        win_path = string(uppercase(lw), ":/", rst)
    else
        lw = split(file_path, '/')[3]
        rst = join(split(file_path, '/')[4:end], '/')
        win_path = string(uppercase(lw), ":/", rst)
    end
    return win_path
end

"""
windows path to wsl
"""
function towsl(file_path::Union{String, Nothing})
    if isnothing(file_path)
        #return error("File path cannot be nothing")
        drive = lowercase(pwd()[1])
        rst = replace(pwd()[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    else
        drive = lowercase(file_path[1])
        rst = replace(file_path[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    end
    return wsl_path
end

#import Pkg
#Pkg.gc(; collect_delay=Second(0))

# #pwd is home
# pcmd = `perl -E '$z=$ENV{PWD}=~ s#mnt\S##r =~s/(\w)/\U$1:/mr =~s[/][]r;say $z'`
#X=`$pwd.Path`
#run(pcmd)
# #readchomp(pipeline(pcmd))
# pwd()

function getdirs(kw...)
    printstyled(filter(x->isdir(x),readdir(kw...)),color=:blue)
end

"""
dfroute(;ofl="route.txt")
reads from routeg(infile, ofl) and returns a DataFrame with the following columns:
    - sim: simulated flow
    - obs: observed flow
    - name: name of the station
"""
function dfroute(;ofl="route.txt")
    df = CSV.read(ofl,DataFrame,header=false,skipto=8,delim="\t",footerskip=1,lazystrings=false)
    rename!(df,1=>"sim",2=>"obs",3=>"name")
    df.name=map(x->replace(x,r"#" => "",r" " => "",r"_>.*" => "",r"-" => "_"),df[:,3])
    sort!(df, :sim)
    return df    
end

"""
non-recursively search for control file in current directory
"""
function ctlx()
    matches::Vector{Any} = []
    for file in readdir(".")
        if endswith(file, ".xml")
            path = joinpath(pwd(), file)
            open(path) do f
                for line in eachline(f)
                    if occursin("compiling symbols in control file ", line)
                        fields = split(line)[8:end]
                        println(join(fields, " "))
                        out = join(fields, " ")
                        push!(matches, out)
                    end
                    if occursin("looking for starting date in ", line)
                        fields = split(line)[9:end]
                        str = replace(join(fields, " "),r"\"/>.*" => " -> ")
                        printstyled("discharge file is: $str\n",color=:light_green)
                    end
                end
            end
        end
    end
    if !isempty(matches)
        fl = first(matches)
        fl = split(fl) |> last
        fl = split(fl, "\"") |> first
        if (!occursin("regio",fl) && occursin("regio",pwd()))
            fl = replace(fl,"control"=>"D:/Wasim/regio/control")
        elseif (!occursin("brend",fl) && occursin("brend",pwd()))
            fl = replace(fl,"control"=>"D:/Wasim/Tanalys/DEM/brend_fab/control")            
        elseif (!occursin("temp",fl) && occursin("saale",pwd()))
            fl = replace(fl,"control"=>"D:/temp/saale/control")
        end
        return string(fl)
    else
        @warn "no control file or xml found!..."
        return ""
    end
end

function qtab(;todf=true)
    dfs = getq();
    # z = ctlx()
    # if isempty(z)
    #     printstyled("no control file found!\n",color=:light_red)
    #     pretty_table(dfs,header=uppercasefirst.(names(dfs));)
    #     return
    # end
    ofl = "route.txt"
    if isempty(ofl)
        printstyled("no $ofl found!\n",color=:light_red)
        pretty_table(dfs,header=uppercasefirst.(names(dfs));)
        return
    end   
    dx = dfroute(;ofl=ofl);
    rename!(dx,"sim"=>"Basin");
    dfs.Basin = parse.(Int64,dfs.Basin)
    kd  = innerjoin(dfs, dx, on=:Basin)
    if todf
        return kd
    else
        pretty_table(kd,header=uppercasefirst.(names(kd));)
    end
end

"""
returns DataFame with qgk Vals recursively
"""
function getq(;prefix::AbstractString="qgko")
    rootdir = pwd()
    results = []
    if any(x -> isdir(x), readdir())
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                if occursin(Regex(prefix, "i"), filename) && !occursin(r"txt|yrly|nc|png|svg", filename)
                    push!(results, joinpath(looproot, filename))
                end
            end
        end
    else
        printstyled("eval on: $rootdir !\n", color=:light_red)
        for filename in filter(x -> isfile(x), readdir(; join=false))
            if occursin(Regex(prefix, "i"), filename) && !occursin(r"txt|yrly|nc|png|svg", filename)
                printstyled("collecting $filename...\n", color=:light_yellow)
                push!(results, filename)
            end
        end
    end

    if isempty(results)
        printstyled("no qgk files found!\n", color=:light_red)
        return 
    end

    for file in results
        x = file
        try
            df = CSV.read(x, DataFrame, delim="\t"; header=true, types=String, silencewarnings=true, skipto=364)
            pattern = r"^[LIN. R]|^[LOG. R]|^CO"
            mask = occursin.(pattern, df[!, 1])
            ddd = df[mask, :]
            new = names(ddd)[5:end]
            insert!(new, 1, "basin")
            insert!(new, 2, "timestep")
            ddd = permutedims(ddd)
            dropmissing!(ddd)
            ddd.basin = new
            select!(ddd, :basin, :)
            df = permutedims(ddd)
            columns = uppercasefirst.(getproperty(df, propertynames(df)[1]))
            rename!(ddd, Symbol.(columns))
            kd = ddd[3:end, :]
            return kd
        catch e
            println(e)
            @warn "skipping $x"
        end
    end

    return 
end

"""

tries to parse all columns to float, except of date.
"""
function dfloat(x::DataFrame)
    df = copy(x)
    for colname in names(df)
        if eltype(df[!, colname]) != Date && eltype(df[!, colname]) == String
            #&& eltype(df[!, colname]) != Float64
            df[!, colname] .= tryparse.(Float64, df[!, colname])
        end
    end
    return df
end    

"""
tries to parse all columns to float, except of date.
"""
function dfloat!(df::DataFrame)
    for colname in names(df)
        if eltype(df[!, colname]) != Date && eltype(df[!, colname]) == String
            df[!, colname] .= tryparse.(Float64, df[!, colname])
        end
    end
end


"""
another reader, faster date conversion
"""
function wadd(x::String)
    ms = ["-9999", "lin", "log", "--"]
    df = CSV.File(x; delim="\t", header=1, normalizenames=true, missingstring=ms) |> DataFrame
    dropmissing!(df,1)
    #df.date = Date.(map(row -> string(row...), eachrow(select(df, 1:3))), DateFormat("yyyymmdd"))
    df.date = Date.(
        map(row -> string(join(row,"-")), 
        eachrow(select(df, 1:3))),
        DateFormat("yyyy-mm-dd"))
    select!(df, Not(1:4))
    metadata!(df, "filename", x, style=:note)
    return df
end


"""
    sf(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())

Search for files in a directory and its subdirectories that match a given pattern.

# Arguments
- `pattern::Union{AbstractString, Regex}`: The pattern to match. Can be a string or a regex.
- `rootdir::AbstractString=pwd()`: The root directory to start the search from. Defaults to the current working directory.

# Returns
- A list of file paths that match the pattern.
"""
function sf(pattern::Union{AbstractString, Regex}, rootdir::AbstractString=pwd())
    # Check if rootdir is a directory
    rootdir = isdir(rootdir) ? rootdir : throw(ArgumentError("rootdir is not a directory!"))    

    # Convert pattern to regex if it's a string
    pattern = pattern isa AbstractString ? Regex("$(pattern)", "i") : pattern

    # Function to check if a filename matches the pattern and add it to results
    function check_and_add(filename, results)
        if occursin(pattern, filename)
            push!(results, filename)
        end
    end

    results = []
    if any(isdir, readdir(rootdir))
        for (looproot, dirs, filenames) in walkdir(rootdir)
            for filename in filenames
                check_and_add(joinpath(looproot, filename), results)
            end
        end
    else
        printstyled("no dirs in $rootdir !\n", color=:light_red)
        for filename in filter(isfile, readdir(rootdir; join=true))
            check_and_add(filename, results)
        end
    end

    return results
end

"""
like Grep.grep("x",df)
"""
function findindf(df::DataFrame, x::Union{AbstractString,Regex})
    filter(row -> any(occursin(x, 
        string(value)) for value in row), 
            eachrow(df))
end


#endof func declaration

println("you are here: ")
printstyled(pwd()*"\n",color=:green)

## if gr fails, try to plot GR.histogram(randn(1000))
# import GR
# function fixgr()
#     GR.histogram(randn(1000))
#     GR.GRPreferences.diagnostics()    
# end
# fixgr()
# function gz(prefix::AbstractString)
#     rootdir="."
#     results = []
#     if (any(x->isdir(x),readdir()))
#         for (looproot, dirs, filenames) in walkdir(rootdir)
#             for filename in filenames
#                 #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#                 if (occursin(Regex(prefix,"i"),filename))
#                     push!(results, joinpath(looproot, filename)) 
#                 end
#             end
#         end
#     else
#         printstyled("no dirs in $rootdir !\n",color=:light_red)
#         for filename in (filter(x->isfile(x),
#             readdir(;join=false)))
#             #if (startswith(filename, prefix)) && (!occursin(r"txt|yrly|nc|png|svg",filename))
#             if (occursin(Regex(prefix,"i"),filename))
#                 push!(results, filename) 
#             end
#         end
#     end
#     return results
# end

# ##in NEW session:
# using PackageCompiler
# pt,nm=("C:/Users/chs72fw/.julia/sysimages",
#     "win_ts_plots.so")
# # joinpath(pt,nm)
# create_sysimage([:StatsPlots,
#     :Plots,:PlotThemes,
#     :DataFrames,
#     :CSV,:Dates],sysimage_path=joinpath(pt,nm))
#julia --startup-file=no --sysimage "C:/Users/chs72fw/.julia/sysimages/win_ts_plots.so" C:\Users\Public\Documents\Python_Scripts\julia\waba-rev.jl
#C:/Users/Public/Documents/Python_Scripts


# using PackageCompiler
# pt,nm=("C:/Users/chs72fw/.julia/sysimages",
#     "win_makie.so")
# create_sysimage([
#     :CairoMakie,
#     :DataFrames,
#     :CSV,:Dates],sysimage_path=joinpath(pt,nm))


# pt,nm=("C:/Users/chs72fw/.julia/sysimages",
#     "windf.so")
# create_sysimage([
#     :Grep,
#     :DataFrames,
#     :CSV,:Dates],sysimage_path=joinpath(pt,nm))

#######################
# sqrt(a) === a^0.5
# log(a^0.5) === 0.5*log(a)
# log(a^2) === 2*log(a)
# m,n = 20,15
# log(m * n) === log(m) + log(n)
# log(m ÷ n) === log(m) - log(n) #f
# log(m / n) === log(m) - log(n)  #but rounded, its true
# round(log(m / n),digits=3) === round(log(m) - log(n),digits=3)  #true

# cd(raw"J:\jras")
# using PackageCompiler
# create_sysimage(;sysimage_path="all.so")
