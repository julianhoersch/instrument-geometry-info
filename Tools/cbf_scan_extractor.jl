# Extract and output scan information from minicbf or ADSC files

# Assumptions:
# 1. files are named in some form of xxxx_scanno(_)frameno.cbf
# 2. frames are sequential
# 3. The first frame will be scan 1 frame 1
# 4. filenames are the same length
# 5. axis names are drawn from the lists below
# 6. "<axis>_increment" signals the increment

using ArgParse

# Configuration information

const rot_axes = ("chi", "phi", "detector_2theta", "two_theta", "omega",
                  "angle", "start_angle", "kappa")
const trans_axes = ("detector_distance", "dx", "trans", "distance", "detector_x", "detector_y",
                    "detector_z")
const always_axes = ("distance", "two_theta", "detector_2theta")

is_trans_axis(axis) = lowercase(axis) in trans_axes
is_rot_axis(axis) = lowercase(axis) in rot_axes

get_scan_info(frame_dir; stem = Regex(".*?_")) = begin

    frame_catcher, all_names, first_scan = get_scan_frame_fmt(frame_dir, stem = stem)

    # index them all for speed

    all_frames = map(all_names) do an
        m = match(frame_catcher, an)
        if haskey(m, "scan")
            (m["scan"],parse(Int, m["frame"])) => an
        else
            ("01", parse(Int, m["frame"])) => an
        end
    end

    all_frames = Dict(all_frames)
    
    all_scans = unique!(map( x -> x[1], collect( keys( all_frames ) ) ) )
    println("$(length(all_scans)) scans found")
    scan_info = Dict()

    axes = union(rot_axes, trans_axes)

    frame_type = determine_frame_type( joinpath(frame_dir, first(all_frames).second))

    @debug "$frame_type determined"
    
    for s in all_scans
        
        snames = filter( x-> x[1] == s, keys(all_frames))
        frames = [ x[2] for x in snames ]

        @debug "$(length(frames)) frames found for scan $s"

        # Get information for first and last frame

        fname = joinpath( frame_dir, all_frames[(s, 1)])

        vals, scan_ax, scan_incr, et, wl = get_frame_info(frame_type, fname, axes)
        start = vals[scan_ax]

        fname = joinpath( frame_dir, all_frames[(s, length(frames))])
        vals, _, _, _, _ = get_frame_info(frame_type, fname, axes)
        finish = vals[scan_ax]

        # Check increment and range match
        
        if !isapprox(start + scan_incr * (length(frames)-1), finish, atol=1e-6)
            throw(error("Scan range does not match increment: $start to $finish, $(length(frames)-1) steps of $scan_incr"))
        end
       
        details = Dict("frames" => length(frames),
                       "axis" => scan_ax,
                       "incr" => scan_incr,
                       "time" => et,
                       "start" => start,
                       "range" => scan_incr * length(frames),
                       "wavelength" => wl)
        
        scan_info[s] = (vals, details)
    end

    prune_scan_info!(scan_info)
    
    return scan_info, all_frames
    
end

# Deduce scan/frame naming convention
get_scan_frame_fmt(frame_dir; stem = Regex(".*?_")) = begin

    first_scan = 1
    
    all_names = readdir(frame_dir)
    filter!(x-> x[end-3:end] in (".cbf",".img"), all_names)
    if !ismissing(stem)
        filter!(x -> startswith(x, stem), all_names)
    end
    sort!(all_names)

    # Analyse number of digits between stem and extension: if less than
    # 5, no scan present

    test_name = all_names[1]
    stem_len = length(match(stem, test_name).match)
    num_digits = count(r"[0-9]", test_name[stem_len:end-4])

    @debug "Num digits for $test_name" num_digits

    if num_digits > 4

        fr_sc_regex = missing
        
        # Allow first scan to not be scan 1

        for t in 1:9
            rr = stem * Regex("(?<scan>[0-9]*$t)(?<sep>0|_)(?<frame>[0-9]+1)\\.(cbf|img)")
            m = match(rr,all_names[1])
            if !isnothing(m)
                
                s_pos = m.offsets[1]
                s_len = length(m["scan"])
                f_len = length(m["frame"])
                if m["sep"] == "0"
                    f_pos = m.offsets[2]
                    f_len = f_len + 1
                else
                    f_pos = m.offsets[3]
                end

                fr_sc_regex = stem * Regex("(?<scan>[0-9]{$s_len})_?(?<frame>[0-9]{$f_len})")

                first_scan = t

                @debug "Found naming scheme for first scan $first_scan" s_len f_len

                break
            end
        end

    else
        
        fr_sc_regex = stem * r"(?<frame>[0-9]+)"

    end

    if ismissing(fr_sc_regex)
        @error("Cannot find scan/frame naming pattern for $test_name. Try using the -s option")
        exit()
    end
    
    @assert !isnothing(match(fr_sc_regex, all_names[end]))
    
    return fr_sc_regex, all_names, first_scan
end

"""
    
Determine the type of frame file: currently SMV (ADSC) and CBF are
recognised.
"""
determine_frame_type(filename) = begin

    header = read(filename, 512)
    if 0x0c in header return Val(:SMV) end
    if occursin("_array_data", String(header)) return Val(:CBF) end
    return missing
end

"""
Return any values found for provided axes. All axes converted to lowercase. This routine
tries to adapt to all of the crazy stuff stashed in miniCBF headers.
"""
get_frame_info(t::Val{:CBF}, fname, axes) = begin
    
    header = readuntil(fname, "--CIF-BINARY-FORMAT-SECTION--")
    lines = lowercase.(split(header, "\n"))

    ax_vals = map(ax -> (lowercase(ax), (get_header_value(t,lines, ax))), axes)
    filter!( x-> x[2] != nothing, ax_vals)
    ax_vals = convert_units.(ax_vals)

    ax_incr = map(ax -> (lowercase(ax), get_header_value(t,lines, ax*"_increment")), axes)
    filter!( x-> x[2] != nothing, ax_incr)
    ax_incr = convert_units.(ax_incr)
    
    et, _ = get_header_value(t, lines, "exposure_time")
    wl, _ = get_header_value(t, lines, "wavelength")

    scan_ax = findfirst( x -> !isapprox(x[2], 0, atol = 1e-6), ax_incr)

    scan_ax_name = indexin([ax_incr[scan_ax][1]], [x[1] for x in ax_vals])

    if scan_ax_name == [nothing]

        # Maybe we have the convention that Oscillation axis gives the name?

        an = get_header_value(t, lines, "oscillation_axis")
        if an != nothing && any(x -> x[1] == an, ax_vals)
            scan_ax_name = an
            filter!(x -> x[1] != "start_angle", ax_vals)
            @debug ax_vals
        else
            @error "Could not match scanned axis $(ax_incr[scan_ax][1]) with an axis name"
        end
    else
        scan_ax_name = ax_vals[scan_ax_name[]][1]
    end

    # Get rid of duplicate names

    ax_vals = Dict(ax_vals)
    println("$ax_vals")

    if haskey(ax_vals, "distance") && haskey(ax_vals, "detector_distance")
        delete!(ax_vals, "detector_distance")
    end

    # Some Pilatus headers do not mention omega, just "start_angle"
    # and "angle_increment".

    if scan_ax_name != "angle"
        delete!(ax_vals, "angle")
    end
    
    if haskey(ax_vals, "start_angle") && haskey(ax_vals, "angle") || scan_ax_name != "angle"
        delete!(ax_vals, "start_angle")
    end
    
    @debug "Extracted info" ax_vals ax_incr scan_ax_name et wl
    
    return Dict(ax_vals), scan_ax_name, ax_incr[scan_ax][2], et, wl
    
end

# For a single-axis diffractometer currently
get_frame_info(t::Val{:SMV}, fname, _) = begin

    header = String(read(fname, 512))
    lines = lowercase.(split(header, "\n"))

    ax_vals = [("phi", get_header_value(t, lines, "phi"))]
    push!(ax_vals, ("trans", get_header_value(t, lines, "distance")))
    ax_incr = [("phi", get_header_value(t, lines, "osc_range"))]
    push!(ax_vals, ("trans", 0.0))

    et = get_header_value(t, lines, "time")
    wl = get_header_value(t, lines, "wavelength")

    return Dict(ax_vals), "phi", ax_incr[1][2], et, wl 
end

"""
    Get the value following the string given in matcher and units if present
"""
get_header_value(::Val{:CBF}, lines, matcher) = begin

    rr = Regex("# $matcher[ =]+")
    one_line = filter( x-> !isnothing(match(rr, x)), lines)
    if length(one_line) != 1
        @debug "Matched lines:" one_line
        return nothing
    end

    one_line = one_line[]
    #@debug "Extracting from" one_line

    m = match(Regex("$matcher[ =]+(?<val>[A-Za-z0-9+-.]+) +(?<units>[A-Za-z.]+)"), one_line)
    val = strip(m["val"])
    units = strip(m["units"])
    #@debug "To get value" val

    try
        return parse(Float64, val), units
    catch
        return lowercase(val), units
    end
end

"""
    Get the value following the string given in matcher and units if present
"""
get_header_value(::Val{:SMV}, lines, matcher) = begin

    rr = Regex("^$matcher[ =]+")
    one_line = filter( x-> !isnothing(match(rr, x)), lines)
    if length(one_line) != 1
        return nothing
    end

    one_line = one_line[]
    #@debug "Extracting from" one_line

    m = match(Regex("$matcher[ =]+(?<val>[A-Za-z0-9+-.]+)"), one_line)
    val = strip(m["val"])

    #@debug "To get value" val

    return parse(Float64, val)
end

"""
   Detect any non-mm translations and convert
"""
convert_units(ax_val) = begin
        name, (val, units) = ax_val
        if units == "m"
            val = val * 1000
        elseif units == "cm"
            val = val * 10
        end
        name, val
end

"""
Axes in `scan_info` and `axis_dict` should already be lowercase.
"""
rename_axes!(scan_info, axis_dict) = begin

    if length(axis_dict) == 0 return end

    for (s, v) in scan_info

        vals, dets = v
        for (a, x) in vals

            if a in keys(axis_dict)
                vals[axis_dict[a]] = x
                delete!(vals, a)
            end
        end

        if dets["axis"] in keys(axis_dict)
            dets["axis"] = axis_dict[dets["axis"]]
        end
    end
    
end

"""
    Remove reference to any axes that do not change position and are
    essentially zero, but are not in `always_axes`.
"""
prune_scan_info!(scan_info) = begin

    initial_vals, deets = first(scan_info).second
    scan_axis = deets["axis"]
    keep_this = [scan_axis]
    for (name, ini_val) in initial_vals
        for (s,(v,d)) in scan_info
            if v[name] != ini_val
                push!(keep_this, name)
                break
            end
        end
    end

    @debug "Changing axes: " keep_this

    for (name, ini_val) in initial_vals
        if !(name in always_axes) && !(name in keep_this) && isapprox(ini_val, 0, atol=0.001)
            @debug "Axis $name left off as unchanging and zero"
            for s in keys(scan_info)
                delete!(scan_info[s][1], name)
            end
        end
    end
    
    
end

#============= Output routines ========================#

"""
We identify each frame by its sequence number overall (not just within
its own scan.)
"""

output_scan_info(scan_info, all_frames, output_file, new_url; prepend_dir = "", arch = nothing, separate = false) = begin

    sl = create_scan_list(scan_info)
    exp_info = Dict([s=>scan_info[s][2]["time"] for s in keys(scan_info)])

    op = isnothing( output_file ) ? stdout : open(output_file, "w")

    generate_wavelength(op, scan_info)
    generate_scan_settings(op, scan_info)
    generate_scan_info(op, sl)
    generate_step_info(op, sl, exp_info)
    generate_array_info(op, sl)
    generate_ids(op, sl)
    generate_external_ids(op, new_url, all_frames, sl, prepend_dir, comp=arch, separate = separate)
    
end

create_scan_list(scan_info) = begin

    # Create scan list of (scanid, frame_no) where
    # frame_no is just the number of frames

    scans = sort(collect(keys(scan_info)))
    slist = [(s, scan_info[s][2]["frames"]) for s in scans]

    @debug "Scan list" slist
    
    return slist
end

generate_wavelength(op, scan_info; rad_type = "xray") = begin
    wl = first(scan_info).second[2]["wavelength"]
    println(op,"_diffrn_radiation_wavelength.id 1")
    println(op,"_diffrn_radiation_wavelength.value $wl")
    println(op,"_diffrn_radiation.type $rad_type")
end

"""

Information we have obtained from the cbf file

        details = Dict("frames"
                       "axis"
                       "incr"
                       "time"
                       "start"
                       "range")
        

"""
generate_scan_settings(op, scan_info) = begin
    header = """loop_
_diffrn_scan_axis.scan_id
_diffrn_scan_axis.axis_id
_diffrn_scan_axis.displacement_start
_diffrn_scan_axis.displacement_increment
_diffrn_scan_axis.displacement_range
_diffrn_scan_axis.angle_start
_diffrn_scan_axis.angle_increment
_diffrn_scan_axis.angle_range
"""
    println(op, header)
    for s in sort(collect(keys(scan_info)))
        
        axes, dets = scan_info[s]
        
        for (a, val) in axes

            step, range = 0, 0
            if a == dets["axis"]
                step = dets["incr"]
                range = dets["range"]
                val = dets["start"]
            end

            if is_trans_axis(a)
                println(op, "SCAN$s  $a $val $step $range . . .")
            else
                println(op, "SCAN$s  $a . . . $val $step $range")
            end
        end
        println(op, "")

    end

    
end

"""

Array_id is just IMAGE (a single detector module). The
binary_id is incremented with frame, and all of them
are located externally so external_id is also incremented
together with binary_id.
"""
generate_ids(op, scan_list) = begin
    header = """loop_
_array_data.array_id
_array_data.binary_id
_array_data.external_data_id
"""
    println(op, header)
    ctr = 0
    for (s,f) in scan_list
        for i in 1:f
            ctr += 1
            println(op, "   IMAGE $ctr $ctr")
        end
    end
    println(op,"")
end

"""

`fulluri` is the location of a single archive file. The individual files
on the local storage are assumed to be at the same locations relative to
this top-level directory.
"""
generate_external_ids(op, fulluri, all_frames, scan_list, prepend_dir; comp="TBZ",fmt="CBF", separate=false) = begin
    header = """loop_
_array_data_external_data.id
_array_data_external_data.format
_array_data_external_data.uri
"""
    if !isnothing(arch)
        header = header *
"""_array_data_external_data.archive_format
_array_data_external_data.archive_path"""
    end

    # If comp is nothing, then each frame has a separate URI and is treated as a local file.
    # If comp is something, then each frame is relative to a single URI and is output
    # separately.
    println(op, header)
    ctr = 0
    sctr = 0
    for (s,f) in scan_list
        sctr += 1
        outuri = fulluri
        if separate
            outuri = insert_scan(fulluri, sctr)
        end
        for i in 1:f
            ctr += 1
            fname = all_frames[(s,i)]
            print(op, "  $ctr $fmt $outuri")

            # A too-clever-by-half way of optionally live constructing a URL
            
            if !isnothing(arch)
                print(op, "  $comp $prepend_dir")
            end
            
            println(op, "/" * "$fname")
        end
    end
    println(op,"")
end

"""
Fill in the scan information. We number the frames from
the start
"""
generate_scan_info(op, scan_list) = begin
    header = """loop_
_diffrn_scan.id 
_diffrn_scan.frame_id_start
_diffrn_scan.frame_id_end
_diffrn_scan.frames"""
    println(op, header)
    stptr = 1
    for (s,f) in scan_list
        endptr = stptr + f - 1
        println(op, "SCAN$s    frm$stptr  frm$endptr  $f")
        stptr = endptr + 1
    end
    println(op,"")
end

"""
Fill in information about steps
"""
generate_step_info(op, scan_list, scan_times) = begin
    header = """loop_
_diffrn_scan_frame.frame_id
_diffrn_scan_frame.scan_id
_diffrn_scan_frame.frame_number
_diffrn_scan_frame.integration_time"""
    println(op, header)
    ctr = 0
    for (s,n) in scan_list
        for f in 1:n
            ctr += 1
            println(op, "frm$ctr   SCAN$s    $f $(scan_times[s])")
        end
    end
    println(op, "")
end

generate_array_info(op, scan_list) = begin
    header = """loop_
_diffrn_data_frame.id
_diffrn_data_frame.detector_element_id
_diffrn_data_frame.array_id
_diffrn_data_frame.binary_id"""
    ctr = 0
    println(op, header)
    for (s,n) in scan_list
        for f in 1:n
            ctr += 1
            println(op, "frm$ctr  ELEMENT  IMAGE $(ctr)")
        end
    end
end

"""
Replace the final `nn...` sequence in `template` with a zero-padded
integer `scan_no`
"""
insert_scan(template, scan_no) = begin
    m = collect(eachmatch(r"n{2,}", template))[end]
    n = length(m.match)
    start = m.offset
    scan_str = "$scan_no"
    zeros = "0"^(n-length(scan_str))
    return template[1:start-1]*zeros*scan_str*template[start+n:end]
end

determine_archive(location) = begin

    arch = nothing
    
    if location == []
        location = "file://" * frame_dir
        
    else
        location = parsed_args["location"][]
        long_end = location[end-6:end]
        short_end = location[end-2:end]
        
        if short_end == "tgz" || long_end == ".tar.gz"
            arch = "TGZ"
        elseif short_end == "tbz" || long_end == "tar.bz2"
            arch = "TBZ"
        elseif short_end == "txz" || long_end == ".tar.xz"
            arch = "TXZ"
        elseif short_end == "zip"
            arch = "ZIP"
        end

    end

    return arch, location
        
end

parse_cmdline(d) = begin

    s = ArgParseSettings(d)

    @add_arg_table! s begin
        "-l", "--location"
        help = "Final URL of files in output. Final 'nn' in name will be replaced by scan number if -n option is used"
        nargs = 1
        "-s", "--stem"
        help = "Constant portion of frame file name. This can help determine the scan/frame file naming convention"
        nargs = 1
        default = [""]
        "-i", "--include"
        help = "Include directory name as part of frame location info in output"
        nargs = 0
        "-n", "--separate"
        help = "Separate scans are in separate archives. Use 'nn' in the -l option to indicate naming."
        nargs = 0
        "-o", "--output"
        help = "Output file to write to, if missing stdout"
        nargs = 1
        "-a", "--axis"
        nargs = 2
        metavar = ["cbf", "new"]
        action = "append_arg"
        help = "Change axis name from <cbf> to <new> in output file to match goniometer axis definitions. May be used multiple times for multiple axis renaming" 
        "directory"
        help = "Directory containing scan frames in minicbf format"
        required = true
    end

    parse_args(s)

end

if abspath(PROGRAM_FILE) == @__FILE__

    # Process arguments
    
    parsed_args = parse_cmdline("Extract scan information from minicbf files")
    frame_dir = parsed_args["directory"]
    do_output = parsed_args["output"] != nothing
    prepend_dir = parsed_args["include"] ? frame_dir : ""
    file_stem = parsed_args["stem"][] == "" ? Regex(".*?_") : Regex(parsed_args["stem"][])
    @debug parsed_args
    
    # Analyse CBF files
    
    scan_info, all_frames = get_scan_info(frame_dir, stem = file_stem)

    # Rename axes

    @debug "axis" parsed_args["axis"]
    if length(parsed_args["axis"]) > 0
        lower_axes = Dict([lowercase(x[1]) => x[2] for x in parsed_args["axis"]])
        @debug "lower" lower_axes
        rename_axes!(scan_info, lower_axes)
    end
    
    # Output CIF fragment

    arch, location = determine_archive(parsed_args["location"])
    
    out_file = parsed_args["output"] != [] ? parsed_args["output"][] : nothing
    if do_output

        output_scan_info(scan_info, all_frames, out_file, location,
                         prepend_dir = prepend_dir,
                         arch = arch,
                         separate = parsed_args["separate"])
    end
end
