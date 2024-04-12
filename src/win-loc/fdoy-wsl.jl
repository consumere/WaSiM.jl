
cd("/mnt/d/remo/qm/prec")
ssup()
@pyjl
plt = pyimport("matplotlib.pyplot")
xr = pyimport("xarray")
nconly()
ds = "pre-add.nc"
if ds isa String
    z = xr.open_dataset(ds)
else
    z = ds
end

k = z.keys()|>collect|>last
@info "key is $k"
idx = z.variables|>collect
if (!any(x->occursin("time",x),z.coords|>collect))
    z = z.rename(Dict("t"=>"time"))
    z = z.set_coords("time")
    @info "set_coords from t to time!"
end
if txy
    lng = "x"
    ltt = "y"
    z = z.assign_coords(t=z.time)
else
    lng = ([i for i in idx if endswith(i, "x") || startswith(i, "lon")] )|>first
    ltt = ([i for i in idx if endswith(i, "y") || startswith(i, "lat")] )|>first
end

if tosum
    grouped_ds = z[k].mean(lng).mean(ltt).groupby("time.dayofyear").sum()
else
    grouped_ds = z[k].mean(lng).mean(ltt).groupby("time.dayofyear").mean()
end        
vr = uppercase(first(k))
nm = basename(ds)
plt.rc("font", family="serif", serif=["cmr10"])
plt.figure()
plt.plot(grouped_ds, label="\$$nm _{$vr}\$")

#filter(x -> occursin(r"name", string(x)), collect(keys(z.attrs)))
#filter(x -> occursin(r"on", string(x)), collect(keys(z.attrs)))
plt.title("$(uppercasefirst(nm))"*ti)
plt.xlim(0, 365)
plt.gca().grid(alpha=0.3)
plt.legend()
plt.show()

@pyjl
nconly("Q")
pyjl.fdoy("QDM_adjust.nc";tosum=true)
nconly("rr-c")|>first|>pyjl.fdoy
@vv "rr-co.nc"
pyjl.fdoy("rr_obs_365.nc";tosum=true)
lat()
a=xr.open_dataset("QDM_adjust.nc")
b=xr.open_dataset("obs.nc")
idx = a.variables|>collect

begin
    k = a.keys()|>collect|>last
    lng = ([i for i in idx if endswith(i, "x") || startswith(i, "lon")] )|>first
    ltt = ([i for i in idx if endswith(i, "y") || startswith(i, "lat")] )|>first
    ag = a[k].mean(lng).mean(ltt).groupby("time.dayofyear").sum()
    bk = b.keys()|>collect|>last
    bg = b[bk].mean(lng).mean(ltt).groupby("time.dayofyear").sum()
    # plt.figure()
    # plt.plot(b[bk], label="\$quantiles\$")
    # plt.show()
end

begin   
    nm="QDM_adjust.nc"
    plt.figure()
    plt.plot(ag, label="\$cor\$")
    plt.plot(bg, label="\$obs\$")
    plt.title(nm)
    plt.xlim(0, 365)
    plt.gca().grid(alpha=0.3)
    plt.legend()
    plt.show()
end
