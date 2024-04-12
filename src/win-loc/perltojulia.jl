using Statistics

# some funcs by Tom Christiansen
# tchrist@perl.com

function by_number(a, b)
    if a < b
        return -1
    elseif a > b
        return 1
    else
        return 0
    end
end

##tchrist@perl
n = 0
sum = 0
values = []
seen = Dict()

# while (<>) {
    # while (/($number_rx)/g) {
        # $n++;
        # my $num = 0 + $1;  # 0+ is so numbers in alternate form count as same
        # $sum += $num;
        # push @values, $num;
        # $seen{$num}++;
    # } 
# } 

#
number_rx = r"(?:(?:[+-]?)(?=[0123456789.])
                    (?:(?:(?:[0123456789]+)(?:(?:[.])
                    (?:[0123456789]*))?)|
                    (?:(?:[.])
                    (?:(?:[0123456789]+))))
                    (?:(?:[Ee])(?:(?:[+-]?)(?:[0123456789]+))|))"

#(path, name, suffix) = fileparts(ARGS[1])
(path, name, suffix) = splitpath(infile)


println("alias pa")
println("script is: $(path)$(name)")
x_bar = 0
x_sd = 0
y_bar = 0
y_sd = 0
i = 0
numerator = 0
r = 0
f1_data = []
f2_data = []
row = []
row2 = []
global k, l
if length(ARGS) < 4
    println("you don't have enough parameters!\nusage: pers <file1> <file2> <column of f1> <column of f2>!\nhint: last column is -1 | first column is 0 ... try pnam or cnam ...\n")
else
    println("command:\t -> ")
    println("pers $(ARGS)")
    #f1_data = readdlm(ARGS[1], Any)[1:end, ARGS[3]]
    #f1_data = CSV.read(ARGS, DataFrame)[1:end, 5]          #eher so...
    f1_data = CSV.read(ARGS[1], DataFrame)[1:end, ARGS[3]]          #eher so...
    k = f1_data[1]
    for row in f1_data
        if occursin(r"^[0-9]", row[1])
            push!(values, row[1])
            seen[row[1]] = get(seen, row[1], 0) + 1
        end
    end

    n = length(values)
    sum = sum(values)
    mean = mean(values)
    stdev = std(values)
    max_seen_count = maximum(values(seen))
    modes = filter(x -> seen[x] == max_seen_count, keys(seen))
    mode = length(modes) == 1 ? modes[1] : "($(join(modes, ", ")))"
    mode *= " @ $max_seen_count"
    sorted_values = sort(values, by_number)
    mid = Int(length(values) / 2)
    median = length(values) % 2 == 0 ? (sorted_values[mid-1] + sorted_values[mid]) / 2 : sorted_values[mid]

    println("stats of $(ARGS[1]), field: $k (targets)")
    println("n is $(n)")
    println("min is $(minimum(values))")
    println("max is $(maximum(values))")
    println("mode is $(mode), median is $(median), mean is $(mean), stdev is $(stdev)")
end
