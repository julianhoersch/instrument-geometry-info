# This routine takes a Github issue and produces an imgCIF fragment,
# which is then passed to Github to create a PR
using Downloads, JSON, Markdown, CrystalInfoFramework, Printf

#=====================

Routines for getting the data out of Github

=====================#

get_issue(issue_number) = begin
    info = Downloads.download("https://api.github.com/repos/COMCIFS/instrument-geometry-info/issues/$issue_number")
    JSON.parse(open(info))["body"]
end

parse_issue(mdown_text) = begin
    sections = Markdown.parse(mdown_text)
    result = Dict{String,Any}()
    for i in 1:2:length(sections)
        result[extract_info(sections[i])] = extract_info(sections[i+1])
        #println("$(sections[i].text):$(sections[i+1])")
    end
    return result
end

extract_info(x) = x

extract_info(p::Markdown.Paragraph) = begin
    join(extract_info.(p.content)," ")
end

extract_info(p::Markdown.Italic) = begin
    p.text[]
end

extract_info(p::Markdown.Header) = p.text[]
extract_info(p::Markdown.List) = begin
    list_vals = Dict{String,Bool}()
    single_item = p.items[]
    for one_item in single_item
        checked = match(r"(?<check>\[.+\])(?<name>.+)",extract_info(one_item))
        list_vals[strip(checked["name"])] = checked["check"] != "[ ]"
    end
    return list_vals
end

"""
Parse a string of form axis,sense,axis,sense...
"""
parse_axis_string(x) = begin
    s = strip.(split(x,","))
    if length(s)%2 != 0
        raise(error("Axis string is incorrect: %x"))
    end
    axes = []
    senses = []
    for i in 1:2:length(s)
        push!(axes,s[i])
        a = lowercase(s[i+1])
        if a == "clockwise"
            a = "c"
        elseif a == "anticlockwise"
            a = "a"
        elseif !(a in ["a","c"])
            raise(error("Unrecognised rotation sense $a"))
        end
        push!(senses,a)
    end
    return axes,senses
end

"""
Convert contents where necessary. Might be a good idea to introduce simplified
keys?
"""
post_process(md_result) = begin
    if haskey(md_result,"Goniometer axes")
        md_result["Goniometer axes"] = parse_axis_string(md_result["Goniometer axes"])
    end
    if haskey(md_result,"Two theta axis")
        md_result["Two theta axis"] = md_result["Two theta axis"] == "clockwise" ? "c" : "a"
    else
        md_result["Two theta axis"] = missing
    end
    return md_result
end

#=======================================

Routines for constructing an imgCIF block

=======================================#

"""
    make_gonio_axes(raw_info)

Given a list of gonio axes, create their representation in imgCIF. The list of gonio
axes goes in order from top to bottom, meaning that the first "depends on" the second
and so forth. We assume a two theta axis.
The items we have to fill in are:
1. type -> rotation
2. depends_on -> next in list
3. equipment -> goniometer
4. vector -> almost always 1 0 0 (rotation about principal axis)
5. offset -> always [0 0 0] but needed for loop integrity

Note that our questions assume looking from above whereas imgCIF is looking from
below, so the sense of rotation is reversed.

"""
make_gonio_axes(raw_info) = begin
    axes,senses = raw_info["Goniometer axes"]
    n = length(axes)
    axis_type = fill("rotation",n)
    equip = fill("goniometer",n)
    kappa_axis = split(get(raw_info,"kappa","- - -"),[' ',','])[1]

    # Construct axis dependency chain

    depends_on = Union{String,Nothing}[]
    append!(depends_on,axes[2:end])
    push!(depends_on,nothing)

    # Remember the principal axis direction

    principal = senses[end]

    # Create direction vectors
    
    vector = Vector{Union{Missing,Real}}[]
    for (a,d) in zip(axes,senses)
        rotfac = d == principal ? 1 : -1
        if a == kappa_axis
            kv = get_kappa_info(raw_info)
            kv[1] *= rotfac
            push!(vector,kv)
        else
            push!(vector,[1,0,0]*rotfac)
        end
    end
    offset = fill(Union{Missing,Real}[0,0,0],n)
    return raw_info["Goniometer axes"][1],axis_type,equip,depends_on,vector,offset
end

#
#
#
get_kappa_info(raw_info) = begin

    # extract kappa information

    axes,senses = raw_info["Goniometer axes"]
    kinfo = split(raw_info["kappa"],[' ',','],keepempty=false)
    if length(kinfo) == 2 push!(kinfo,"0") end
    axname = kinfo[1]
    if !(axname in axes)
        throw(error("Kappa axis $axname not listed as a goniometer axis"))
    end
    kapang = parse(Float64,kinfo[2])
    kappoise = parse(Float64,kinfo[3])
    
    # Now calculate direction of axis in X-Z plane assuming
    # rotation of kappa is same as principal axis. There is no
    # Y component as home position is always under
    # beam path.

    up_comp = cosd(kapang)
    across_comp = sind(kapang)
    if kappoise == 0
        #is under incident beam collimator in X-Z plane
        z_comp = -1*across_comp
    elseif kappoise == 180
        z_comp = across_comp
    end
    return [up_comp,0.0,z_comp]
end

"""
        make_detector_axes(raw_info)

    Add information concerning the detector axes. We define our own axis names,
    with the detector distance being inserted when the data file is read. We
    choose det_x to be in the horizontal direction, and det_y to be vertical.
    We need to add:
    1. type -> translation
    2. depends_on -> x,y depend on translation
    3. equipment -> detector
    4. vector -> worked out from user-provided info
    5. offset -> beam centre

    
"""
make_detector_axes(raw_info) = begin

    # Insert assumed axis names and orientations
    
    axis_id = ["two_theta","trans","detx","dety"]
    axis_type = fill("translation",4)
    axis_type[1] = "rotation"
    equip = fill("detector",4)
    depends_on = [nothing,"two_theta","trans","detx"]

    # Read necessary information
    
    principal_sense = raw_info["Goniometer axes"][2][end]
    principal = raw_info["Principal axis orientation"]
    corner = raw_info["Image orientation"]

    # Adjust two theta direction

    rotsense = raw_info["Two theta axis"] == principal_sense ? 1 : -1
    vector = [[rotsense,0,0]]

    # Detector translation always opposite to beam direction
    
    push!(vector,[0,0,-1])

    # Work out det_x and det_y

    x_d,y_d = determine_detx_dety(principal,principal_sense,corner)
    push!(vector,x_d)
    push!(vector,y_d)
    
    # Beam centre and distance is unknown for now

    offset = [[0,0,0],[0,0,missing],[missing,missing,0],[0,0,0]]
    return axis_id,axis_type,equip,depends_on,vector,offset
end

# Determine direction of detx (horizontal) and dety (vertical) in
# imgCIF coordinates.

determine_detx_dety(principal_angle,principal_sense,corner) = begin
    
    # Start with basic value and then flip as necessary

    x_direction = [-1,0,0]        # spindle rotates anticlockwise at 0, top_left origin
    y_direction = [0,1,0]       #
    if corner == "top right"
        x_direction *= -1
    elseif corner == "bottom right"
        x_direction *= -1
        y_direction *= -1
    elseif corner == "bottom left"
        y_direction *= -1
    end

    @debug "Before pa adjustment" x_direction y_direction
    
    # The direction of the principal axis flips by 180 if the sense
    # changes

    pa = principal_sense == "a" ? parse(Int64,principal_angle) : parse(Int64,principal_angle) + 180
    if pa >= 360 pa = pa - 360 end

    @debug "Principal axis direction" pa
    
    if pa == 90
        temp = x_direction
        x_direction = y_direction
        y_direction = -1*temp
    elseif principal_angle == 180
        x_direction *= -1
        y_direction *= -1
    elseif principal_angle == 270
        temp = x_direction
        x_direction = -1*y_direction
        y_direction = temp
    end

    @debug "After pa adjustment" x_direction y_direction
    
    return x_direction,y_direction
end

"""
    split_vector

Distribute the vectors in `vector` over CIF data names,
formatting as integers if they are exact.
"""
split_vector(vector::Vector,basename,imgblock) = begin
    for i in 1:3
        imgblock[basename*"[$i]"] = map(vector) do x
            if ismissing(x[i]) x[i]
            else
                if round(x[i]) == x[i]
                    @sprintf "%d" x[i]
                else
                    @sprintf "%1.5f" x[i]
                end
            end
        end
    end
end

"""
    describe_axes(raw_info,imgblock)

Create the goniometer axes corresponding to the data in `raw_info`,
placing the result in CIF block `imgblock`.
"""
describe_axes(raw_info,imgblock) = begin
       
    gon_axes = make_gonio_axes(raw_info)
    det_axes = make_detector_axes(raw_info)
    for i in 1:length(gon_axes)
        @debug "$(gon_axes[i]) -- $(det_axes[i])"
        append!(gon_axes[i],det_axes[i])
    end
    base = "_axis."
    imgblock[base*"id"] = gon_axes[1]
    imgblock[base*"type"] = gon_axes[2]
    imgblock[base*"equipment"] = gon_axes[3]
    imgblock[base*"depends_on"] = gon_axes[4]
    split_vector(gon_axes[5],base*"vector",imgblock)
    split_vector(gon_axes[end],base*"offset",imgblock)
    create_loop!(imgblock,[base*"id",base*"type",base*"equipment",
                          base*"depends_on",
                          base*"vector[1]",base*"vector[2]",base*"vector[3]",
                          base*"offset[1]",base*"offset[2]",base*"offset[3]"])
end

get_pixel_sizes(raw_info) = begin
    both = parse.(Float64,split(raw_info["Pixel size"],","))
    if length(both) == 1 push!(both,both[1]) end
    return both
end

"""
    describe_detector(raw_info)

Produce the information required for array_structure_list. Here
we assume a rectangular detector with x horizontal, y vertical
"""
describe_array(raw_info,imgblock) = begin
    hor,vert = get_pixel_sizes(raw_info)
    # array structure list axis
    base = "_array_structure_list_axis."
    imgblock[base*"axis_id"] = ["detx","dety"]
    imgblock[base*"axis_set_id"] = ["1","2"]
    imgblock[base*"displacement"] = [hor/2,vert/2]   #half pixel size
    create_loop!(imgblock,[base*"axis_id",base*"axis_set_id",
                      base*"displacement"])
    
    # array structure list
    base = "_array_structure_list."
    imgblock[base*"array_id"] = ["1","1"]
    imgblock[base*"index"] = [1,2]
    imgblock[base*"axis_set_id"] = ["1","2"]
    imgblock[base*"dimension"] = [missing,missing]   #number of elements in each direction
    imgblock[base*"direction"] = ["increasing","increasing"]
    if raw_info["Fast direction"] == "horizontal"
        precedence = [1,2]
    else
        precedence = [2,1]
    end
    imgblock[base*"precedence"] = precedence
    create_loop!(imgblock,[base*"array_id",base*"index",base*"axis_set_id",
                      base*"dimension",base*"direction",base*"precedence"])
end

describe_detector(raw_info,imgblock) = begin
    base = "_diffrn_detector."
    imgblock[base*"id"]=["1"]
    imgblock[base*"number_of_axes"]=[2]
    create_loop!(imgblock,[base*"id",base*"number_of_axes"])
    #
    base = "_diffrn_detector_axis."
    imgblock[base*"axis_id"] = ["detx","dety"]
    imgblock[base*"detector_id"] = ["1","1"]
    create_loop!(imgblock,[base*"axis_id",base*"detector_id"])
end

describe_facility(raw_info,imgblock) = begin
    base = "_diffrn_source."
    if haskey(raw_info,"Beamline name")
        imgblock[base*"beamline"] = [raw_info["Beamline name"]]
        imgblock[base*"facility"] = [raw_info["Facility name"]]
    else
        imgblock[base*"make"] = [raw_info["Name of manufacturer"]*"-"*raw_info["Model"]]
        if haskey(raw_info,"Location")
            imgblock[base*"details"] = ["Located at $(raw_info["Location"])"]
        end
    end
end


make_block_id(raw_info) = begin
    repregex = r"[^A-Za-z0-9_]"
    if haskey(raw_info,"Facility name")
        fname = replace(raw_info["Facility name"],repregex=>"_")
        bline = replace(raw_info["Beamline name"],repregex=>"_")
        block_id = fname*"_"*bline
    else
        mname = replace(raw_info["Name of manufacturer"],repregex=>"_")
        model = replace(raw_info["Model"],repregex=>"_")
        loc = replace(get(raw_info,"Location",""),repregex=>"_")
        block_id = mname*"_"*model
    end
    return block_id
end

const output_order = ("_audit.block_id",
                      "_diffrn_source.beamline",
                      "_diffrn_source.facility",
                      "_axis.id","_axis.type",
                      "_axis.equipment",
                      "_axis.depends_on",
                      "_axis.vector[1]",
                      "_axis.vector[2]",
                      "_axis.vector[3]",
                      "_axis.offset[1]",
                      "_axis.offset[2]",
                      "_axis.offset[3]",
                      )

"""
    create_cif_block(issue_number)

Given a Github issue number, create a block of CIF text for inclusion
in the list of beamline/instrument geometries
"""
create_cif_block(issue_number) = begin
    c = post_process(parse_issue(get_issue(issue_number)))
    create_cif_block(c)
end

create_cif_block(d::Dict{String,Any}) = begin  #from markdown
    block_id = make_block_id(d)
    # for future: check for duplicates here
    my_block = CrystalInfoFramework.Block{Any}()
    my_block["_audit.block_id"] = [block_id]
    describe_axes(d,my_block)
    describe_facility(d,my_block)
    describe_array(d,my_block)
    describe_detector(d,my_block)
    return my_block
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia issue_to_imgcif.jl [-r] <issue_number>")
        println("""
<issue_number> is an issue number in the Github instrument-geometry-info
repository that has been created by the automatic submission tool.
     Use option -r to see the raw information downloaded from Github and
     -rr to see the parsed information.""")
    elseif length(ARGS) == 2
        if ARGS[1] == "-r"
            issue_number = ARGS[2]
            println(issue_number |> get_issue)
        elseif ARGS[1] == "-rr"
            issue_number = ARGS[2]
            println(issue_number |> get_issue |> parse_issue |> post_process)
        else
            throw(error("Incorrect option, only -r recognised"))
        end
    else
        issue_number = ARGS[1]
        final_block = issue_number |> get_issue |> parse_issue |> post_process |> create_cif_block
        show(stdout,MIME("text/cif"),final_block,ordering=output_order)
    end
end
