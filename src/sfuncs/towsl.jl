# 
function towsl(file_path::Union{String, Nothing})
    if isnothing(file_path)
        #return error("File path cannot be nothing")
        drive = lowercase(pwd()[1])
        rst = replace(pwd()[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    else
        drive = lowercase(file_path[1])
        rst = replace(file_path[3:end], "\\" => "/")
        wsl_path = string("/mnt/", drive, rst)
    end
    return wsl_path
end

#import Pkg
#Pkg.gc(; collect_delay=Second(0))

# #pwd is home
# pcmd = `perl -E '$z=$ENV{PWD}=~ s#mnt\S##r =~s/(\w)/\U$1:/mr =~s[/][]r;say $z'`
#X=`$pwd.Path`
#run(pcmd)
# #readchomp(pipeline(pcmd))
# pwd()

