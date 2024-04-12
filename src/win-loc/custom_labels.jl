#custom_labels.jl
using Plots
# Create some test data
x = 1:10
y = rand(10)

# Define a function to generate custom labels
function custom_label(i,y)
    return "Point $i: $y"
end

# Plot the data
Plots.scatter(x, y);

# Annotate every data point
for i in 1:length(x)
    annotate!(x[i], y[i], 
    Plots.text(custom_label(i,y[i]), :bottom, 8, :blue))
end

# Show the plot
plot!()


wa.baryrsum(td)
function custom_label(y)
    f = round(y, digits=2)
    #return "Tdiff[mm]\n $f"
    return "$f [mm]"
end
#u = unique(td.year)
# left_margin = 10mm,
# bottom_margin = 10mm,
tds = yrsum(td)
mx = maximum(tds[!,2])
for i in 1:length(tds.year)
    #annotate!(tds.year[i], tds.Tdiff[i]+maximum(tds.Tdiff[i]), 
    annotate!(tds.year[i], mx, 
    Plots.text(custom_label(tds.Tdiff[i]), :top, 8, :black, rotation=45))
end
plot!()