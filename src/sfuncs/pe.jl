# 
function pe()
    try
    inp = clipboard()
    wpath = replace(inp, "\\" => "/", "\"" => "")
    cmd = `wslpath -ua $wpath`
    ot = readchomp(pipeline(cmd))
    clipboard("$ot")
    println("$ot in clipboard!")
    return string(ot)
    catch e
        println("Error: $e")
        println("Failed to translate to wslpath.")
    end
end