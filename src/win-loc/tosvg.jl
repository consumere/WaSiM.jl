#List Windows Subsystem für Linux-Distributionen
run(`wsl -l`)

cd(dirname(filename))
inf = "abfluss.pgm"
out = "abfluss.svg"

cmd = `wsl -d Ubuntu-18.04 -e mkbitmap -f 32 -t 0.4 "$inf" -o out.bmp`
run(cmd)
cmd = `wsl -d Ubuntu-18.04 -e potrace "$inf" --svg -o "$out"`
run(cmd)


###JPG to PGM to SVG#######
using RCall
x="C:/Users/chs72fw/Downloads/abfluss.jpg"
of="C:/Users/chs72fw/Downloads/abfluss.pgm"
    @rput x
    @rput of
R"
i = magick::image_read(x)
pg = magick::image_convert(i,format = 'pgm')
magick::image_write(pg,path = of,quality = 100)
"
cd("C:/Users/chs72fw/Downloads/")
cmd1 = `wsl -d Ubuntu-18.04 -e mkbitmap -f 32 -t 0.4 "$inf" -o -`
cmd2 = `wsl -d Ubuntu-18.04 -e potrace --svg -o "$out"`
run(pipeline(cmd1, cmd2))
op()
###########################show image in julia################
#img = load(out) ##er
@rput out
R"
i = magick::image_read(out)
plot(i)
"
using Images
filename = "C:/Users/chs72fw/Downloads/abfluss.pgm"
img = load(filename)


#################wsl julia################
cmd1 = (`wsl -d Ubuntu-18.04 -e ls -l`)
cmd2 = (`wsl -d Ubuntu-18.04 -e head -n 5`)
run(pipeline(cmd1, cmd2)) #works.
run(`wsl -d Ubuntu-18.04 -e ls -l| head -n 5`) #errors.
run(`wsl -d Ubuntu-18.04 -e head -n 4 <(ls -l)`) #errors.

run(`wsl -d Ubuntu-18.04 -e awk -h`)
#run(`wsl -d Ubuntu-18.04 -e source /home/cris/.bash_aliases; pyll`)

#alias pearson='eval $bspt/pearson.pl'
per="/mnt/c/Users/Public/Documents/Python_Scripts/bashscripts/pearson.pl"
run(`wsl -d Ubuntu-18.04 -e $per`) #works.
pwd()
raw"D:\Wasim\regio\out\rc200\x17\full2"|>cd

f1="Bad_Kissingen_Golfplatz-qoutjl"
cmd1 = (`wsl -d Ubuntu-18.04 -e $per $f1 $f1 5 4`)
run(cmd1)

cmd1 = (`wsl -e $per $f1 $f1 5 4`)
run(cmd1)


#########ps in julia###############
cdb()
run(`powershell.exe -NoLogo -Command fdm`) 
run(`powershell.exe -NonInteractive -NoLogo -Command fdm`) 
sdf()



getnames(outd)
typeof(outd)
zp(getnames)
a=wa.getnames(outd)[1]
b=wa.getnames(outd)[2]
pwd()
cmd1 = (`wsl -e $per $a $b 5 5`)
cmd1 = (`wsl -e $per $a $b 4 4`)
cmd1 = (`wsl -e $per $a $b 5 4`)
run(cmd1)
cmd1 = (`wsl -e $per $a $a 5 4`)
run(cmd1)



fl=raw"D:\Wasim\regio\rcm200\v10\kleinhans_pars.txt"
df = CSV.read(fl,DataFrame;delim=" ")


using PyCall
# import the necessary Python modules
statsmodels = pyimport("statsmodels.api")
pylab = pyimport("pylab")
# generate some random data
data = randn(100)
# create the QQ plot
statsmodels.qqplot(data, line="45")



#sys:1: UserWarning: Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.
#--> matplotlib.use("TkAgg")
matplotlib.use("TkAgg")

df=waread(r"Wol")
#@pyimport matplotlib.pyplot as plt

x, y = df[:,1], df[:,2]
ax = plt.subplot(111)
plt.scatter(x, y)
ax.set_xlabel(names(df)[1])
ax.set_ylabel(names(df)[2])
statsmodels.qqline(ax, "r", x, y)
plt.show()


using PyCall

# define a Julia function that calls the Python function waread3
function waread3_julia(x, flag=true)
    # use py"""...""" syntax to access the Python function
    py"""
    import pandas as pd

    def waread3(x, flag=True):
        if flag:
            df = pd.read_csv(x, delim_whitespace=True, header=0,
                             na_values=-9999, verbose=True,engine='c')
            if 'YY' not in df.columns:
                print("Column 'YY' not found in the CSV file.")
                return None
            if df.iloc[:, 0].dtype == 'object':
                print('First column contains strings, subsetting to Int...')
                df = df[~df.iloc[:, 0].str.contains("[A-z]|-", na=False)]
            source_col_loc = df.columns.get_loc('YY')        
            df['date'] = df.iloc[:, source_col_loc:source_col_loc +
                                 3].apply(lambda x: "-".join(x.astype(str)), axis=1)
            df = df.iloc[:, 4:]
            df['date'] = pd.to_datetime(df['date'])
            df.set_index('date', inplace=True)
            df.iloc[:,0:-2] = df.iloc[:,0:-2].apply(lambda x: x.astype(float), axis=1)
            df.filename = x
            print(df.filename,"done!")
            return df
        else:
            print('Date parse failed, try reading without date transformation...')
            df = pd.read_csv(x, delim_whitespace=True, comment="Y", skip_blank_lines=True).dropna()
            df.filename = x
            print(df.filename,"done!")
            return df
    """
    # call the Python function with the arguments and return the result
    return py"waread3($x, $flag)"
end

# test the Julia function with some example arguments
s=@gl "Bad"
df = waread3_julia(s)
djl = wa.pydf(df)

py"pf = $waread3(s, True)"

# import the necessary Python modules
pd = pyimport("pandas")
ddf = pd.read_csv(s, delim_whitespace=true, header=0,na_values=-9999, verbose=true)
#ddf.columns.get_loc("Y")
#ddf["date"] = ddf.iloc[:, 0:3].apply(lambda x: "-".join(x.astype(str)), axis=1)

ddf = wa.pydf(ddf)
dt = Date.(map(x -> join(x, "-"), eachrow(ddf[:, 1:3])))
typeof(dt)
nd = hcat(ddf[!, 5:end], dt)
rename!(nd,ncol(nd) =>"date")
dfp(nd)


nd = pd.read_table(s,
                    delim_whitespace=true,
                    skiprows=0,
                    na_values=-9999,
                    skip_blank_lines=true,
                    #encoding='cp1252',
                    #verbose=true,
                    #low_memory=false,
                    #infer_datetime_format=true,
                    #format="%YY %mm %dd",
                    parse_dates = ["YY", "MM", "DD"])
nd.filename = s
#format="%YY %mm %dd",
nd.reset_index(inplace=true)
#pd.to_datetime(nd["index"], format="%YY %mm %dd")


typeof(nd)
nd.head()
nd.columns

nd.drop(["YY", "MM", "DD", "HH"], axis=1, inplace=true)
pd.plotting.scatter_matrix(nd, alpha=0.2, figsize=(6, 6), diagonal="kde")
plt.show()
pnd = wa.pydf(nd)


m=@gl "Sch"
df = pyread(m)