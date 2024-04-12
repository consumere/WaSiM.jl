#ClimateExplorer check

#this is : sed -i 's/0.3000000E+34/-9999/g' $a
#"/mnt/d/ClimateExplorer/pressure"

x="glob_10.6293___49.7787_30.txt"
file = open(x, "r")
# Read in the contents of the file
contents = readlines(file)
# Close the file
close(file)

# Replace the string "0.3000000E+34" with "-9999"
for i in 1:length(contents)
    contents[i] = replace(contents[i], "0.3000000E+34" => "-9999")
end
# Open the file for writing
file = open("out.txt", "w")
# Write the modified contents back to the file
for line in contents
    write(file, line*"\n")
end
# Close the file
close(file)
#rm("out.txt")

file=open("out.txt", "r")
cs = readlines(file)
close(file)