# if length(ARGS) == 0
# 	println("need args! <xmlfile>...")
#     exit()
# end


#fn = filter(x -> occursin(r"xml", x), readdir())

rootdir=pwd()
fn = []
for (looproot, dirs, filenames) in walkdir(rootdir)
    for filename in filenames
        if occursin(r"xml", filename) && !occursin(r"txt|yrly|nc|png|svg", filename)
            push!(fn, joinpath(looproot, filename))
        end
    end
end

if length(fn) == 0
	println("no xmlfile found...")
    exit()
end


using Dates
#to cb...
function extract_duration_from_xml(xml_file::AbstractString)
    finished_timestamp = raw""
    start_timestamp = raw""

    # Read the content of the XML file
    content = readlines(xml_file)

    # Find the lines containing "WaSiM finished" and "WaSiM start"
    finished_lines = filter(x -> occursin(r"WaSiM finished", x), content)
    start_lines = filter(x -> occursin(r"WaSiM start", x), content)

    # Extract the timestamps from the lines if found
    if !isempty(finished_lines)
        finished_line = first(finished_lines)
        finished_timestamp = match(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}", finished_line).match
    end

    if !isempty(start_lines)
        start_line = first(start_lines)
        start_timestamp = match(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}", start_line).match
    end

    # Calculate the duration in hours and minutes
    if !isempty(finished_timestamp) && !isempty(start_timestamp)
        finished_time = DateTime(finished_timestamp)
        start_time = DateTime(start_timestamp)
        duration_milliseconds = Dates.value(finished_time - start_time)

        # Convert duration to hours and minutes
        duration_hours = div(duration_milliseconds, 3_600_000)  # 1 hour = 3,600,000 milliseconds
        duration_minutes = div(rem(duration_milliseconds, 3_600_000), 60_000)  # 1 minute = 60,000 milliseconds

        hrs,min = (duration_hours, duration_minutes)
        println("eval on: $xml_file")

#        if tocb
            #message = "# run took $hrs hrs and $min min..."
            message = "run took $hrs hrs and $min min...\n"
            printstyled(message, color=:green)
            #clipboard(message)
            #printstyled("msg in clipboard! \n",color=:green)
            return nothing
        # else
        #     return hrs, min
        # end
                    
    else
        return nothing
    end
end

#extract_duration_from_xml(ARGS[1])
if length(fn) > 1
    extract_duration_from_xml.(fn)
else
    extract_duration_from_xml(fn[1])
end
