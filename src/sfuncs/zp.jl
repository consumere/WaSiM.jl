#  which python 
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
