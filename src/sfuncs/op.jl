#  explorer . 
function op()
    #pwrs""" explorer . """
    #open(`powershell -noprofile explorer . `,"w",stdout)
    #open(`cmd.exe /c start . `,"w",stdout)
    #this wrks in wsl, too
    run(`cmd.exe /c start .`)
end 

macro pwp_str(s) open(`powershell`,"w",stdout) do io; print(io, s); end;end
macro cmd_str(s) open(`cmd \c`,"w",stdout) do io; print(io, s); end;end
##nope, bad idea
# run(`cmd \c pwd";" exit`)
