# 
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
