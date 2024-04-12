#soilscript

function readbetween(io::IO, start::String, stop::String)
    output = Vector{String}()
    while !eof(io)
        line = readline(io)
        if contains(line, start)
            push!(output, line)
            while !eof(io)
                line = readline(io)
                if contains(line, stop)
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

using DataFrames
   

    #out=[]
    # for (i,field) in enumerate(data);
    #     if occursin("Name=",field) ;
    #     d = Dict(numbers[i]=>field);
    #     push!(out,d)
    # end;end
    return out
end


"""
skips first line after [soil_table] i.e. no of soil types
now returns a DataFrame
"""
function rd2(filename::String)
    # Read lines between "soil_table" and "substance_transport"
    data = open(filename) do io
        readbetween(io, "soil_table", "substance_transport")
    end

    # Filter out empty lines
    #data = filter(x -> !isempty(x), data[3:end])
    data = data[3:end]
    numbers = extract_numbers(data)
    # Remove comments and unnecessary characters
    data = broadcast(x -> replace(x,    
            r"^#.*" => "",
            r"^[[].*" => "",
            r"method" => "",
            r"MultipleHorizons" => "",
            r"}" => "",
            r" = " => "=",
            #r"[?*.{=;]" => "",      #problems
            r";" => "" ), data)

            # data = replace.(data, [
            #     r"^#.*|^\[\[.*" => "",
            #     r"method|MultipleHorizons|[}]" => "",
            #     r" = |[?*.{=;]" => "",
            #     r";" => ""
            # ])
            
    
    
#    filter!(x -> !isempty(x), data)
    data = Grep.grep(r"^(?! Evap.*$|^[0-9].*$)",data)
    data = strip.(data)
    filter!(x -> !isempty(x), data)
    # initialize a dictionary to store the data for this line
    out = Dict{Any, Any}()
    # store the number in the dictionary
    #out["number"] = numbers #nope in loop!
    #for (i,field) in enumerate(data);if occursin("Name=",field) ;println([field,numbers[i]]);end;end
    for (i,field) in enumerate(data)
        if occursin("=", field) #&& startswith(field, " ") == true
            #field = data[1]
            out["number"] = numbers[i]
            # split the field into key and value
            key, value = split(field, "=")
            key = strip(key)
            #key = strip(";",key)
            value = strip(value)
            
            # check if the key is "Name"
            if key == "Name"
                # keep the value as a string
            else
            try # convert the value to an array of floats
                # check if the value is a number
                if occursin(r"^-?\d+(\.\d+)?$", value)
                    # convert the value to a float
                    value = parse.(Float64, split(value))
                elseif occursin(r"^-?\d+(\.\d+)?\d+?$", value)
                    # convert the value to a float
                    value = parse.(Float64, split(value))
                elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                    # convert the value to a float (scientific notation)
                    value = parse.(Float64, split(value))
                elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                    # convert the value to a float (scientific notation)
                    value = parse.(Float64, split(value))
                elseif occursin(r"^\d+ \d+", value)
                    # convert the value to an array of integers
                    value = parse.(Int, split(value))
                    #value = parse.(Float64, split(value))
                elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)

                    value = parse.(Float64, split(value))
                                            
                end
            catch
                @warn "could not parse $value"
                continue
            end

            end
            # store the key-value pair in the dictionary
            out[key] = value
        else
            # if the field does not contain "=", skip it
            
            continue
        end
    end
    #DataFrame([merge(parse_line(line), Dict("number" => num)) for (line, num) in zip(data[2:end], numbers) if !isempty(line)])
    #result = DataFrame(out)
    #return result
    return #[numbers, DataFrame(out)]
end

function parse_line(line::AbstractString)
    number, fields = split(strip(line), " {", limit=2)

    number = parse(Int, number[1])
    fields = split(fields, ';')

    data = Dict{String, Any}("number" => number)

    for field in fields
        if occursin("=", field)
            key, value = split(strip(field), "=")

            if key != "Name"
                try
                    value = parse_value(value)
                    data[key] = value
                catch
                    @warn "could not parse $value"
                end
            end
        end
    end

    return data
end

function parse_value(value::AbstractString)
    if occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
        return parse(Float64, value)
    elseif occursin(r"^\d+ \d+", value)
        return parse.(Int, split(value))
    elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)
        return parse.(Float64, split(value))
    else
        return value
    end
end

parsed_lines = [parse_line(line) for line in data if !isempty(line)]

infile = "D:/Wasim/regio/control/rcm200_x22-re-v2.ctl"
df = rd2(infile)


filename = infile
data = open(filename) do io
    readbetween(io, "soil_table", "substance_transport")
end

# Filter out empty lines
data = filter(x -> !isempty(x), data[3:end])

function extract_numbers(lines::Vector{String})
    numbers = []
    for line in lines
        match_result = match(r"(\d+)\s*{method", line)
        if match_result !== nothing
            push!(numbers, parse(Int, match_result.captures[1]))
        end
    end
    return numbers
end
z = extract_numbers(data)


ka = rd2(infile)


data = open(filename) do io
    readbetween(io, "soil_table", "substance_transport")
end
#write to a tempfile
dat = join(data[3:end], " ")
tempfile = "tx.txt"
open(tempfile, "w") do io
    write(io, dat)
end

wa.read_soildata_2("tx.txt")
wa.read_soildata_raw("tx.txt")
wa.read_soildata_4("tx.txt")

data = DelimitedFiles.readdlm("tx.txt",' ', String)
fn="tx.txt"
data = CSV.read(fn,DataFrame,delim=";",header=false,transpose=true)
#dat = replace.(data[3:end],";"=>"")
#now transform

data = open(filename) do io
    readbetween(io, "soil_table", "substance_transport")
end
data = data[3:end]
numbers = extract_numbers(data)
# Remove comments and unnecessary characters
data = broadcast(x -> replace(x,    
        r"^#.*" => "",
        r"^[[].*" => "",
        r"method" => "",
        r"MultipleHorizons" => "",
        r"}" => "",
        r" = " => "=",
        r"[?*.{=;]" => "",      #problems
        r";" => "" ), data)
filter!(x -> !isempty(x), data)
#split(data," }")|>DataFrame
data
#th = Grep.grep(r"^th.*|^[0-9]",data)
th = Grep.grep(r"^(?! Evap.*$|^[0-9].*$)",data)



d = Dict()

"""
reads controlfile
uses Grep.grep to select lines
returns a DataFrame
"""
function ctlook(infile::String)
    data = open(filename) do io
        readbetween(io, "soil_table", "substance_transport")
    end
    data = data[3:end]
    numbers = extract_numbers(data)
    # Remove comments and unnecessary characters
    data = broadcast(x -> replace(x,    
            r"^#.*" => "",
            r"^[[].*" => "",
            r"method" => "",
            r"MultipleHorizons" => "",
            r"}" => "",
            r" = " => "=",
            #r"[?*.{=;]" => "",      #problems
            r";" => "" ), data)
    data = Grep.grep(r"^(?! Evap.*$|^[0-9].*$)",data)
    data = strip.(data)
    filter!(x -> !isempty(x), data)
    sel = Grep.grep(r"^Name", data)
    # Split each string into "Name" and "Value"
    split_name_value = split.(sel, '=', limit=2)
    
    ks = split.(Grep.grep(r"^ksat", data), '=', limit=2)
    ths = split.(Grep.grep(r"^theta_sat", data), '=', limit=2)
    thr = split.(Grep.grep(r"^theta_res", data), '=', limit=2)
    thk = split.(Grep.grep(r"^thickness", data), '=', limit=2)
    parn = split.(Grep.grep(r"^Par_n", data), '=', limit=2)
    #ks = Dict(first(ks)[1] => getindex.(ks, 2))
    # Create a DataFrame
    df = DataFrame( id = numbers,
                    Name = getindex.(split_name_value, 2),
                    ksat = getindex.(ks, 2),
                    theta_sat = getindex.(ths, 2),
                    theta_res = getindex.(thr, 2),
                    Par_n = getindex.(parn, 2),
                    Thick = getindex.(thk, 2))
    return df
end


infile = "D:/Wasim/regio/control/st-v1.ctl"
infile = "D:/Wasim/Tanalys/DEM/brend_fab/control/fab_c8-loc3.ctl"
df = ctlook(infile)
r = findindf(df,"1010")
val = DataFrame(r)|>permutedims
(Dict.(names(df)=>string.(val.x1)))
m=hcat(names(df),DataFrames.stack(r))
#Dict(names(df)=>string.(val.x1))
nd = DataFrame(m,:auto)
v = parse.(Float64,nd[end,2]|>split)
sum(v)
cumsum(v)

#thickness barplot:
using Plots
using StatsPlots
using DataFrames
df = ctlook(infile)
df
id, th = df.id, df.Thick
sums = [sum(parse.(Float64, split(line))) for line in th]
bar(sums .* 0.1, legend=false, xlabel="Soil type",
    #xticks = id, 
    ylabel="Thickness [m]", title="Soil thickness")


ndf= DataFrame(
    id=id,thick=sums .* 0.1)

ssel = @rsubset ndf :thick > 5
@df ssel boxplot(:thick)

filename=infile

m="SoilTillage"
Grep.grep(m, readlines(filename))
lines = Grep.grep(r"d+\s*method|^Albedo|^rsc\s*", readlines(filename))
lines = Grep.grep(r"d+\s*method|^Albedo|^rs", readlines(filename))
filter!(x -> !isempty(x), lines)
lines = replace.(lines,r"{.*" => "")


pwd()
filename=raw"C:\Users\Public\Documents\Python_Scripts\lu.txt"
df = read_landuse_table()


# Extract class names
class_names = Grep.grep(r"^[0-9] ", data)
class_names = replace.(class_names, r"{.*" => "")
class_names = strip.(class_names)

# Create a dictionary to store DataFrames
mns = [last(line) for line in (classnames)]
classnames_words = [match(r"\w+", classname).match for classname in classnames]


class_dataframes = Dict{String, DataFrame}()
mns = [line for line in data if occursin(r"(?i)method", line) == true]
#mns = [occursin(r"\W",x) for x in mns]
filter(x->occursin(r"\W",x),mns)
[split(x,r"\W") for x in mns]

classnames = unique(strip.(replace.(Grep.grep(r"^[0-9] ", data), r"{.*" => "",r"\s+" => " ")))
classnames = replace.(classnames, "f W" => "f_W")

#classnames = replace.(classnames, " " => "_","_" => " ",count=2)
#classnames = replace.(classnames, r"\s+" => " ")
#split.(classnames,limit=1)
#[isodd(length(v)) for v in classnames]
classnames_split = split.(classnames)
df = DataFrame(id = parse.(Int, getindex.(classnames_split, 1)),
               name = getindex.(classnames_split, 2))

classnames = [split.(classname) for classname in classnames]
classnames = [match(r"\w+", classname).match for classname in classnames]

# header = ["Parameter Value"]
# filtered_lines = vcat(header, filtered_lines)
# nms = CSV.File(IOBuffer(join(classnames, '\n')),delim=" ",
#      ignorerepeated=true)|>DataFrame

#dat = join(data, " ")
#csv_file = CSV.File(IOBuffer(join(dat, ' ')), delim='=', quotechar=' ', ignorerepeated=true)
#class_df = DataFrame(csv_file)


# Process each class
for class_name in dat
    # Filter lines for the current class
    lines_for_class = [line for line in class_name if occursin(class_name, line)]
    # Filter lines that do not contain "method"
    filtered_lines = [line for line in lines_for_class if occursin(r"(?i)method", line) == false]

    # Add a header
    header = ["Parameter = Value"]
    filtered_lines = vcat(header, filtered_lines)

    # Create a CSV File from the filtered lines
    csv_file = CSV.File(IOBuffer(join(filtered_lines, '\n')), delim='=', quotechar=' ', ignorerepeated=true)

    # Convert the CSV File to a DataFrame
    class_df = DataFrame(csv_file)

    # Store the DataFrame in the dictionary
    class_dataframes[class_name] = class_df
end

# Access the DataFrames for each class
for (class_name, class_df) in class_dataframes
    println("Class Name: $class_name")
    display(class_df)
    println()
end

pt
io = open(pt, "r")
result = readbetween(io, r"{ method", r"}")
close(io)
result

data = open(filename) do io
    readbetween(io, "landuse_table", 
        "special_output")
end
cn = [line for line in data if occursin(r"(?i)method", line) == true]
class_dataframes = []
io = open(pt, "r")
for x in cn
    classmatch = x[1:10] #first 10 chars
    result = readbetween(io, Regex(classmatch), r"}")
    push!(class_dataframes,result)
end
close(io)
class_dataframes

io = open(pt, "r")
k=readbetween(io, Regex("Shrub"), r"}")
close(io)
k

function read_landuse_data(filename::AbstractString)
    data = open(filename) do io
        readbetween(io, "landuse_table", "special_output")
    end

    class_dataframes = []
    io = open(filename, "r")

    for line in data
        if occursin(r"(?i)method", line)
            classmatch = strip(line[3:15])  # first 10 chars, stripping leading/trailing whitespaces
            println(classmatch)
            result = readbetween(io, Regex(classmatch), r"}$")
            push!(class_dataframes, result)
        end
    end
    close(io)
    return class_dataframes
end

# Usage
filename = "lu.txt"  # replace with your actual filename
pt = "C:\\Users\\Public\\Documents\\Python_Scripts\\lu.txt"
result = read_landuse_data(pt); #works!

filename = infile
result = read_landuse_data(infile); #works!
result[end]
result[1]



"""
reads controlfile
select landuse table
returns a Vector of DataFrames
"""
function read_landuse_data2(filename::AbstractString)
    data = open(filename) do io
        readbetween(io, "landuse_table", "special_output")
    end
    data = broadcast(x -> replace(x,
    r"\\t" => "",        #strip comments
    r"^#.*" => "",        #strip comments
                        r";.*$" => "",         # match everything after semicolon until the end of the line
                        r" = " => "=",      #strip spaces around "="
                        r";" => "" ), data)
                        #r"^[[].*" => "",      #strip module names
                        #r"method" => "",     #startflag
                        #r"}" => "",          #this is needed for end flag
                        
                        #r"[?*.{=;]" => "",      #problems


    class_dataframes = []
    io = open(filename, "r")

    for line in data
        if occursin(r"(?i)method", line)
            classmatch = strip(line[3:15])  # first 10 chars, stripping leading/trailing whitespaces
            println(classmatch)
            filtered_lines = readbetween(io, Regex(classmatch), r"}$")
            filtered_lines = replace.(filtered_lines, r";.*" => "",r"\s+"=>" ")
            filtered_lines = strip.(filtered_lines)
            filtered_lines = filter(x -> !isempty(x), filtered_lines)
            # Add a header
            header = ["Parameter=Value"]
            filtered_lines = vcat(header, filtered_lines)
            # Create a CSV File from the filtered lines
            csv_file = CSV.File(IOBuffer(join(filtered_lines, '\n')), delim='=', quotechar=' ', ignorerepeated=true)
            # Convert the CSV File to a DataFrame
            class_df = DataFrame(csv_file)
            push!(class_dataframes, class_df)
        end
    end
    close(io)
    return class_dataframes
end


pt
dfs = read_landuse_data2(pt);
dfs[1]
dfs[9]
#zogr = Grep.grep(r"^[0-9] |Z0", data)
#@rsubset dfs[9] :1

sel = findindf(dfs[7],"Z")|>DataFrame|>q->select(q,2) 
value_str = findindf(dfs[7],"Z").Value
z0 = [parse.(Float64, substr) for substr in split.(value_str)]
using Dates
cls = split(dfs[7][1,1],r"\W")[2]
month_abbr = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"];
plot(month_abbr,only(z0),title="Roughness length in m",label=cls)
    
#"class $cls"

using Plots
# Assuming dfs is a vector of dataframes
plots = []
month_abbr = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]
plot();
for i in 1:size(dfs,1)
    value_str = findindf(dfs[i], "Z").Value
    #value_str = replace.(value_str, r" " => "")
    value_str = string.(value_str)
    value_str = strip.(value_str)
    vs = split.(value_str)
    z0 = [parse.(Float64, z) for z in vs]
    cls = split(dfs[i][1, 1], r"\W")[2] 
    plot!(month_abbr, only(z0), label=cls)
    #push!(plots, plot(month_abbr, only(z0), label=cls))
end
title!("Roughness length in m")

plot();
for i in 1:size(dfs,1)
    value_str = findindf(dfs[i], "LAI").Value
    #value_str = replace.(value_str, r" " => "")
    value_str = string.(value_str)
    value_str = strip.(value_str)
    vs = split.(value_str)
    z0 = [parse.(Float64, z) for z in vs]
    cls = split(dfs[i][1, 1], r"\W")[2] 
    scatter!(month_abbr, only(z0), label=cls,yaxis=:log10)
    #push!(plots, plot(month_abbr, only(z0), label=cls))
end
title!("LAI in [m2/m2]")


function luscatter(dfs::Vector=dfs,x::String="rs_evaporation")
    month_abbr = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]
    p1 = plot(title=x);
    for i in 1:size(dfs,1)
        value_str = findindf(dfs[i], x).Value
        #value_str = replace.(value_str, r" " => "")
        value_str = string.(value_str)
        #value_str = strip.(value_str)
        vs = split.(value_str)
        z0 = [parse.(Float64, z) for z in vs]
        cls = split(dfs[i][1, 1], r"\W")[2] 
        scatter!(month_abbr, first(z0), label=cls,
            yaxis=:identity) #log10
        #push!(plots, plot(month_abbr, only(z0), label=cls))
    end
    #title!("LAI in [m2/m2]")
    Plots.show(p1)
    return p1    
end

luscatter(dfs,"rs")
luscatter(dfs,"rs_evap")
luscatter(dfs,"LAI")

findindf(dfs[i], "evapo")

ot = []
for i in 1:size(dfs,1)
    value_str = findindf(dfs[i], "rs").Value
    value_str = string.(value_str)
    #value_str = strip.(value_str)
    vs = split.(value_str)
    z0 = [parse.(Float64, z) for z in vs]
    cls = split(dfs[i][1, 1], r"\W")[2]
    push!(ot,DataFrame(cls=cls,rs_evap=first(z0)))
end

a = map(x->findindf(x,"VCF"),dfs)
a = vcat(a...)|>DataFrame
a.nm = vcat(map(x->split(x[1, 1], r"\W")[2],dfs)...)
a
## Vegetation covered fraction
plot(a.Value,legend=a.nm,title="VCF")

value_str = string.(a.Value)
vs = split.(value_str)
z0 = [parse.(Float64, z) for z in vs]
m = hcat(z0...)'
heatmap(m)
xticks!(1:12,month_abbr)
yticks!(1:size(m,1),a.nm)
title!("VCF")


pt #lu.txt
dfs = read_landuse_data2(pt);
size(dfs)
for df in dfs
    @show df[1, 1]
end
#now lai heatmap
a = map(x->findindf(x,"LAI"),dfs)
a = vcat(a...)|>DataFrame
# a.nm = vcat(map(x->split(
#     x[1, 1], r"\W")[2],
#     dfs)...)
a.nm = vcat(map(x->split(
        replace(string(x[1, 1]),r"\s+"=>r"\s"),r"\W")[1],
        dfs)...)

theme(:dao)
begin 
    value_str = string.(a.Value)
    vs = split.(value_str)
    z0 = [parse.(Float64, z) for z in vs]
    m = hcat(z0...)'
    heatmap(m,c=:viridis)
    xticks!(1:12,month_abbr)
    yticks!(1:size(m,1),a.nm)
    title!("LAI Heatmap")
    #add annotation for each field
    for i in 1:size(m,1)
        for j in 1:size(m,2)
            tcolor = m[i,j] > 10 ? :black : :white
            annotate!(j,i,text(string(
                Int(round(m[i, j]; digits=0)) #no floats...
                ),8,tcolor))
        end
    end
    plot!()
end
# str=dfs[end][1,1]
# match(r"\S+", str).match
# words = split(str)
# second_word = length(words) >= 2 ? words[2] : "nomatch"
# vcat(map(x->split(x[1, 1], r"\W")[2],dfs)...)


#now Albedo heatmap
a = map(x->findindf(x,"Albedo"),dfs)
a = vcat(a...)|>DataFrame
a.nm = vcat(map(x->split(
        replace(string(x[1, 1]),
        r"\s+"=>" "),
        r"\W")[1],
        dfs)...)
gr()
begin 
    value_str = string.(a.Value)
    vs = split.(value_str)
    z0 = [parse.(Float64, z) for z in vs]
    m = hcat(z0...)'
    myc = colormap("Grays",10;logscale=true)
    #heatmap(m,c=:grays)
    heatmap(m,c=myc[3:9])
    xticks!(1:12,month_abbr)
    yticks!(1:size(m,1),a.nm)
    title!("Snow Free Albedo")
    #add annotation for each field
    for i in 1:size(m,1)
        for j in 1:size(m,2)
            tcolor = m[i,j] < .15 ? :black : :white
            annotate!(j,i,text(string(
                round(m[i, j]; digits=2)
                ),7,tcolor,rotation=-45))
        end
    end
    plot!()
end


vgctl("modis lai")
vgr("MODIS LAI")
vgr("MODIS")
vgctl("MODIS")
vgctl("modis")


# Plotting
pyplot()
begin 
    heatmap(m, c=:grays, color=:viridis, linecolor=:black, linewidth=0.5)

    # Customize x-axis and y-axis ticks
    xticks!(1:12, month_abbr)
    yticks!(1:size(m, 1), a.nm)

    # Set plot title
    title!("Snow Free Albedo")

    # Add annotation for each field
    for i in 1:size(m, 1)
        for j in 1:size(m, 2)
            tcolor = m[i, j] > 10 ? :black : :white
            annotate!(j, i, text(string(round(m[i, j]; digits=2)), 7, tcolor, rotation=-45), color=:viridis)
        end
    end

    # Display the colorbar with adjusted values
    c = plot!()
    colorbar!(c, label="Adjusted Values", c=:viridis, ticks=0:2:20)
    
end

#colormap(cname="Grays",N=20;logscale=true)


vgjl("LAI ")
fn=raw"D:\Wasim\regio\control\lu_table_modis.txt"
da = read_landuse_data2(fn)
#now lai heatmap
a = map(x->findindf(x,"LAI"),da)
a = vcat(a...)|>DataFrame
a.nm = vcat(map(x->split(
        replace(string(x[1, 1]),r"\s+"=>r"\s"),r"\W")[1],
        dfs)...)

luscatter(da,"LAI")
luscatter(da,"rs")

a = map(k->findindf(k,Regex("lai","i")),da)
a = vcat(a...)|>DataFrame
z = vcat(map(x->split(replace(
            string(x[1, 1]),
            r"\s+"=>" ",
            r"{"=>"",
            r"method"=>"")),
            dfs)...)
a.nm = filter(x->length(x)>2,z)

function luheat(dfs::Vector=dfs,x::String="LAI")
    a = map(k->findindf(k,Regex(x,"i")),dfs)
    a = vcat(a...)|>DataFrame
    z = vcat(map(x->split(replace(
        string(x[1, 1]),
        r"\s+"=>" ",
        r"{"=>"",
        r"method"=>"")),
        dfs)...)
    a.nm = filter(x->length(x)>2,z)
    a = a[:, sort(names(a),rev=true)] #sort columns
    begin 
        value_str = string.(a.Value)
        vs = split.(value_str)
        z0 = [parse.(Float64, z) for z in vs]
        m = hcat(z0...)'
        #myc = colormap("Grays",10;logscale=true)
        myc = colormap("Grays",25;logscale=false)
        #heatmap(m,c=:grays)
        heatmap(m,c=myc[10:end-5])
        xticks!(1:12,month_abbr)
        yticks!(1:size(m,1),a.nm)
        #title!("Snow Free Albedo")
        title!("$x")
        #add annotation for each field
        for i in 1:size(m,1)
            for j in 1:size(m,2)
                tcolor = m[i,j] < .15 ? :black : :white
                annotate!(j,i,text(string(
                    round(m[i, j]; digits=2)
                    ),7,tcolor,rotation=-45))
            end
        end
        return plot!()
    end
end

luheat(da,"LAI")
wa.luheat(da,"rsc")
wa.luheat(da,"vcf")
wa.luheat(da,"dept")
wa.luheat(da,"Albedo")
wa.luheat(da,"rs_interception")
da

k=raw"D:/Wasim/regio/control/vt200-s3-loc.ctl"
dfs = read_landuse_data2(k);
dfs[1]
dfs = filter(x->nrow(x)>3,dfs)
wa.luheat(dfs,"Alb")
wa.luheat(dfs,"lai")

hombr()
cdu()
da=findlog()
hydro(da[1,end])
hydro(da[2,end])
hydro(da[52,end])

da[10,end]|>fread|>ftplin
da[52,end]|>fread|>ftplin
da[end,end]|>fread|>ftplin
da[end,end]|>fread|>ftp
cdof(joinpath(pwd(),da[end,end]))
op()
ctl()
infile = wa.ctl2()
xf = @gl "xml"
grep_with_context("routing_model", infile, 1) #2
##non - observed runoff! 
grep_with_context("[SpinUp", infile, 1)
grep_with_context("[groundwater_f", infile, 1)
hrs,min = wa.extract_duration_from_xml(xf)
message = "run took $hrs hrs and $min min..."
println(message)
clipboard(message)
#run took 2 hrs and 27 min...
qbb()
dfp(r"qbas")
npp(infile)
#D:/Wasim/Tanalys/DEM/brend_fab/

"D:/Wasim/Tanalys/DEM/brend_fab/"|>cd
ezplot()

raw"D:/Wasim/Tanalys/DEM/brend_fab/out/w3/penman/re/"|>cd
qbb()
dfp(r"qbas")
waba2()
dfp(r"gwst")
hydromon(r"qges")
hydro(r"qgko")
hydro(r"qges")