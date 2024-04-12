function read_until_flag(file::IOStream, flag::String)
    line = readuntil(file, flag)
    return line[1:end-length(flag)]
end

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

function read_soildata(filename::String)
    """
    skips first line after [soil_table] i.e. no of soil types
    now returns a DataFrame
    """
    data = open(filename) do io
        a = readbetween(io, "soil_table", "substance_transport")
        return(a)
    end
    data = broadcast(x -> replace(x,    
    r"^#.*" => "",
    r"^[[].*" => "",
    r"method" => "",
    r"MultipleHorizons" => "",
    r"}" => ""), data)
    filter!(s -> !isempty(s), data)
    lines = data[2:end]
    result = []
    # iterate over each line
    for line in lines
        # skip empty lines
        if isempty(line)
            continue
        end
        
        # split the line into number and fields
        number, fields = split(line, " {", limit=2)
        
        # convert the number to an integer
        number = parse(Int, number)
        
        # split the fields into individual fields
        fields = split(fields, ';')
        
        # initialize a dictionary to store the data for this line
        data = Dict{String, Any}()
        
        # store the number in the dictionary
        data["number"] = number
        
        # iterate over each field
        for field in fields
            # check if the field contains the " = " substring
            if occursin(" = ", field)
                # split the field into key and value
                key, value = split(field, " = ")
                
                # check if the key is "Name"
                if key == "Name"
                    # keep the value as a string
                else
                    # check if the value is a number
                    if occursin(r"^-?\d+(\.\d+)?$", value)
                        # convert the value to a float
                        value = parse(Float64, value)
                    elseif occursin(r"^-?\d+(\.\d+)?(e-?\d+)?$", value)
                        # convert the value to a float (scientific notation)
                        value = parse(Float64, value)
                    elseif occursin(r"^\d+ \d+", value)
                        # convert the value to an array of integers
                        value = parse.(Int, split(value))
                    elseif occursin(r"^\d+\.\d+ \d+\.\d+", value)
                        # convert the value to an array of floats
                        value = parse.(Float64, split(value))
                    end
                end
                
                # store the key-value pair in the dictionary
                data[key] = value
            end
        end
        
        # append the dictionary to the result array
        push!(result, data)
    end

    return DataFrame(result)
end


fn = "D:/Wasim/Tanalys/DEM/Input_V2/soil_table.txt"

out[!,5:end]
names(out)

##extract certain columns
ks = []
for value in out.ksat
    push!(ks,parse.(Float64, split(value)))
end

#vars = DataFrames.combine(df, names(df) .=> var) ByRow
hcat(out.number,out.Name,ks)


dk = DataFrame(ks,:auto)
rename!(dk,out.Name,makeunique=true)
qplot(dk)

#rowsums:
DataFrames.combine(dk, names(dk) .=> sum)

rowmeans(dk)

function rowsums(df::DataFrame)
    return DataFrames.combine(df, names(df) .=> sum)
end

function rowmeans(df::DataFrame)
    return DataFrames.combine(df, names(df) .=> mean)
end

rs = rowsums(dk)
vector = vec(Matrix(rs))
Plots.boxplot(vector,legend=false,xticks=false)
Plots.violin!(vector,legend=false,xticks=false)
Plots.bar(vector)
Plots.histogram(vector)

#str = [ @sprintf("%02i", x) for x in (month.(df.date)) ];
ln = propertynames(dk)[20:30]
@df dk StatsPlots.boxplot(cols(ln),fillalpha=0.75, 
linewidth=0.25,legend=false,xticks=false)


@df df StatsPlots.violin(str,cols(ln),linewidth=0.01,legend=false);
@df df StatsPlots.dotplot!(str,cols(ln),fillalpha=0.5,marker=(:black,stroke(1)),legend=false)

function varsort(df::DataFrame)
    rs = DataFrames.combine(df, names(df) .=> var)
    rs = vec(Matrix(rs))
    return sort!(rs,rev=true)
end

k = varsort(dk)
#sort(k,o=reverse)

hcat(out.number,out.Name,k) #NO!

df = out
nam= :ksat
rs = DataFrames.combine(df, 4 .=> var)

rs = vec(Matrix(rs))
return sort!(rs,rev=true)

function varsort(df::DataFrame,nam::Symbol)
    rs = DataFrames.combine(df, names(nam) .=> var)
    rs = vec(Matrix(rs))
    return sort!(rs,rev=true)
end

function unbindvariables()
    for name in names(Main)
        if !isconst(Main, name)
            Main.eval(:($name = nothing))
        end
    end
end
#This simple function may be used to unbind the objects of all non-constant globals in Main
unbindvariables()

fdi()
dd()


fn = "D:/Wasim/regio/control/rcm200_x9.ctl"
out = wa.read_soildata(fn)
##extract certain columns
ks = []
for value in out.ksat
    push!(ks,parse.(Float64, split(value)))
end

out.ksat = ks

#combine(gd, :b, AsTable([:b, :c]) => ByRow(extrema) => [:min, :max])
df = out
vardf = DataFrames.combine(df, :ksat => ByRow(var) )

df = hcat(df,vardf)
##Soildat with higest variance
sort!(df,:ksat_var,rev=true)
wa.writedf("D:/Wasim/regio/control/soil_table_var.csv",df)
#using Base: contains

fn = "D:/Wasim/regio/control/soil_table_var.csv"
xd = CSV.File(fn) |> DataFrame
# formatting: alt + shift + f

#vars = DataFrames.combine(df, names(df) .=> var) ByRow
#hcat(out.number,out.Name,ks)
first(df,5)
dbt = CSV.File("D:/Wasim/main/db_table.csv") |> DataFrame
rename!(dbt,1=>"number")
dropmissing!(dbt)
first(dbt,10)
#bk.join<- merge(x = bk200, y=out, by='TKLE_NR')
xj = innerjoin(dbt,df,on=:number)
outerjoin(dbt,df,on=:number)

# using Shapefile
# shapefile_name="D:/Bodendaten/buek200_2020/bk200_sfmerge.shp"
# shp = Shapefile.Handle(shapefile_name)
# shp.shapes[1:end]
#shp|>plot
using GeoDataFrames
gdf = GeoDataFrames.read(shapefile_name)
rename!(gdf,3=>"number")
gdf.number = 

map(x->parse.(Int64,x),gdf.number)
map(x->round(x;digits=0),gdf.number)


wa.vgctl("transient zone for rain-snow")
ctlg("transient zone for rain-snow")
wa.getm(ctlg)
m=("transient zone for rain-snow")
m=("temperature limit for rain")
m="snow"
ctlg("D:/Wasim/Testdata/Control/","ctl",m)

fn="D:/Wasim/regio/rcm200/v6/rcm.aq1"
wa.agheat(fn)