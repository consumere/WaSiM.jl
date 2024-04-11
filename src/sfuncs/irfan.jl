# 
function irfan(fl::String)
        opener="C:/Program Files/IrfanView/i_view64.exe"
        fl_escaped = realpath(fl)
        cmd = pipeline(`$opener $fl_escaped`)
        run(cmd)
    end

    