# 
function cdb()
#     if Sys.islinux()
#         opwd = `cd (OLDPWD)`
#         run(opwd)
#         pwd()|>println
#     else
#         try
#             #opwd=(`powershell -noprofile Set-Location $env:OldPWD`) 
#             #Cmd(["pwsh","-noprofile", "-Command", "Set-Location", "$env:OldPWD"])
#             # run(opwd)
#             # pwd()|>println
#             prevdir = pwd()
#             run(`cmd /c cd /d ..`)
#             println("Previous directory: ", prevdir)
#             println("Current directory: ", pwd())
#         catch e
#             println(e)
#             @error "pwrsh errored..."  
#         end
#     end
# end

# prev_location = ENV["OLDPWD"]
# run(`cd $prev_location`)
# # List PowerShell's Environmental Variables
# Get-Childitem -Path Env:* | Sort-Object Name

#like jdd to vector of strings.
