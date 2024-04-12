if length(ARGS) == 0
	println("need args! <ncfileregex>...")
    exit()
end

using PyCall

function nconly(x1::AbstractString)
    v::Vector{String} = readdir();
    v = v[broadcast(x->endswith(x,"nc"),v)];
    z = v[(broadcast(x->occursin(Regex(x1),x),v))] 
    return(z)
end

function frx(x::AbstractString;maskval=0,lyr=0)

    x = nconly(x)|>last

py"""
ncf = $x;
ncf = ncf if ncf.endswith('.nc') else print('not an NetCDF grid!\ntry rioplot',ncf) & exit(0)
lyr=$lyr
maskval = float($maskval)
import xarray as xr
from matplotlib import pyplot
import numpy as np
from pandas import DataFrame
dx=xr.open_dataset(ncf)
m=list(dx.variables)[-1]
print('cf attrs:',dx[m].attrs)
df=dx[m].quantile(q=[0,.05,.1,.25,.5,.75,.9,.95,1]).to_dataframe()
print('fieldmean:',dx[m].mean().to_numpy())
print(list(dx.variables),'\n',df,'\n','-'*18)
idx = list(dx.indexes._coord_name_id)
lng = ''.join([i for i in idx if (i.startswith('lon') or i.endswith('x'))])
lat = ''.join([i for i in idx if (i.startswith('lat') or i.endswith('y'))])
tim = ''.join([i for i in idx if (i.startswith('t'))])
if idx == ['x', 'y', 't']:
    print('this is a wasim grid... transposing now...','\n','-'*18)
    dx = dx.transpose().to_array().isel(t=lyr)
elif idx == ['x', 'y', tim]:
    dx = dx.rename_dims({tim:'t'}).transpose('x', 'y', 't')
    dx = dx.to_array().isel(t=lyr)
else:
    dx = dx.rename_dims({tim:'t',lng:'x',lat:'y'})
    print('renamed coord_name_ids: ',idx,'to',dx.dims,'!!!')
    dx = dx.isel(t=lyr).to_array()
y=dx.where(dx.values>maskval)
y=y if not np.all(np.isnan(y.data)) else print('all data of selection is NaN\n exiting now from',ncf) & exit(0)
print('fieldmean of selection:',y.mean().to_numpy())
df=y.variable.quantile(q=[0,.05,.1,.25,.5,.75,.9,.95,1])
ar=np.array([0,.05,.1,.25,.5,.75,.9,.95,1])
df=DataFrame(data=df,index=ar).reset_index().rename({'index':'q',0:m},axis=1)
print('Quantiles of selection:\n',df,'\n','-'*18)
pyplot.rcParams['figure.figsize']=(10,9)
y.plot(cmap='turbo')
pyplot.tight_layout()
title=(ncf.rsplit('/')[-1]+' Layer: '+str(lyr))
pyplot.title(title)
pyplot.show()
"""
end


fn = ARGS[1];
frx(fn)