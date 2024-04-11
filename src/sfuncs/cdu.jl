# 
function cdu()
    dirname(pwd())|>cd
    pwd()|>println
end


visited_dirs = [pwd()]

"""
go to last dir
"""
