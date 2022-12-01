import re
import numpy as np
import CifFile as cf



OUTPUT_ORDER = (
    "_audit.block_id",
    "_diffrn_source.beamline",
    "_diffrn_source.facility",
    "_diffrn_source.make",
    "_diffrn_source.details",
    "_axis.id",
    "_axis.type",
    "_axis.equipment",
    "_axis.depends_on",
    "_axis.vector[1]",
    "_axis.vector[2]",
    "_axis.vector[3]",
    "_axis.offset[1]",
    "_axis.offset[2]",
    "_axis.offset[3]",
    )


def convert_user_input_to_imgcif(user_input, cif_block):

    user_input = post_process(user_input)

    # cif_file = cf.CifFile()
    # cif_block = cf.CifBlock()
    # cif_file['imgCIF'] = cif_block
    # block_id = make_block_id(user_input)
    # cif_block["_audit.block_id"] = [block_id]

    describe_facility(user_input, cif_block)
    describe_axes(user_input, cif_block)
    describe_array(user_input, cif_block)
    describe_detector(user_input, cif_block)


    print('citems', cif_block.items())

    # for idx, item in enumerate(OUTPUT_ORDER):
    #     try:
    #         cif_file['imgCIF'].ChangeItemOrder(item, idx)
    #     except ValueError:
    #         continue

    # print('cblovck', cif_block)
    # print(cif_file.WriteOut(wraplength=1000))#, maxoutlength=2048))
    # with open("mycif.cif", 'w') as file:
    #     file.write(cif_file.WriteOut())



def post_process(user_input):
    """
    Convert contents where necessary. Might be a good idea to introduce simplified
    keys?
    """
    if "Goniometer axes" in user_input:
        user_input["Goniometer axes"] = parse_axis_string(user_input["Goniometer axes"])
    if "Two theta axis" in user_input:
        # TODO is the julia code correct??
        # user_input["Two theta axis"] = user_input["Two theta axis"] == "clockwise" ? "c" : "a"
        if user_input["Two theta axis"] == "clockwise":
            user_input["Two theta axis"] = "c"
        elif user_input["Two theta axis"] == "anti-clockwise":
            user_input["Two theta axis"] == "a"
    else:
        user_input["Two theta axis"] = 'none'

    return user_input


def parse_axis_string(x):
    """
    Parse a string of form axis,sense,axis,sense...
    """
    s = [sub.strip() for sub in x.split(',')] #strip.(split(x, ","))
    # print('mysub', s)
    if len(s) % 2 != 0:
        raise Exception("Axis string is incorrect: %x")

    axes = []
    senses = []
    for i in range(0, len(s), 2): #1:2:length(s)
        axes.append(s[i])
        a = s[i+1].lower()
        # print('my a', a)
        if a == "clockwise":
            a = "c"
        elif a == "anticlockwise":
            a = "a"
        elif a not in ["a", "c"]:
            raise Exception(f"Unrecognised rotation sense {a}")

        senses.append(a)


    return axes, senses


def make_block_id(raw_info):
    repregex = r"[^A-Za-z0-9_]"
    if "Facility name" in raw_info:
        facility = re.sub(repregex, '_', raw_info["Facility name"])
        beamline = re.sub(repregex, '_', raw_info["Beamline name"])
        block_id = facility + "_" + beamline
    elif "Name of manufacturer" in raw_info:
        manufacturer = re.sub(repregex, '_', raw_info["Name of manufacturer"])
        model = re.sub(repregex, '_', raw_info["Model"])
        location = re.sub(repregex, '_', raw_info["Location"])
        block_id = manufacturer + "_" + model

    return block_id



def describe_axes(raw_info, imgblock):
    """
    describe_axes(raw_info,imgblock)

    Create the goniometer axes corresponding to the data in `raw_info`,
    placing the result in CIF block `imgblock`.
    """

    if 'gonio_axes_nxsmx' in raw_info:
        gon_axes = make_detector_axes_nxmx(raw_info)
        print('gonaxes2', gon_axes, type(gon_axes))
    else:
        gon_axes = make_gonio_axes(raw_info)
        print('gonaxes', gon_axes, type(gon_axes))
        det_axes = make_detector_axes(raw_info)
        print('detax', det_axes)
        # structure:  [['Phi', 'Chi', 'Omega'], ['rotation', 'rotation', 'rotation']]
        # adds the detector axes
        for i, _ in enumerate(gon_axes):
            gon_axes[i] += det_axes[i]

    print("gonaxes2", gon_axes)
    base = "_axis."
    imgblock[base + "id"] = gon_axes[0]
    imgblock[base + "type"] = gon_axes[1]
    imgblock[base + "equipment"] = gon_axes[2]
    imgblock[base + "depends_on"] = gon_axes[3]
    split_vector(gon_axes[4], base + "vector", imgblock)
    split_vector(gon_axes[-1], base + "offset", imgblock)
    imgblock.CreateLoop([base + "id", base + "type", base + "equipment",
                         base + "depends_on",
                         base + "vector[1]", base + "vector[2]", base + "vector[3]",
                         base + "offset[1]", base + "offset[2]", base + "offset[3]"])
    print("nnt")
    # print(imgblock.items())


    # [[1, 0, 0], [-1, 0, 0], [1, 0, 0], [[-1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]]], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [[0, 0, 0], [0, 0, 0], [None, None, 0], [0, 0, 0]]])
    # [[1, 0, 0], [-1, 0, 0], [1, 0, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [missing, missing, 0], [0, 0, 0]])


def make_gonio_axes(raw_info):

    """
    make_gonio_axes(raw_info)

    Given a list of gonio axes, create their representation in imgCIF. The list
    of gonio axes goes in order from top to bottom, meaning that the first
    "depends on" the second and so forth. We assume a two theta axis.
    The items we have to fill in are:
    1. type -> rotation
    2. depends_on -> next in list
    3. equipment -> goniometer
    4. vector -> almost always 1 0 0 (rotation about principal axis)
    5. offset -> always [0 0 0] but needed for loop integrity

    Note that our questions assume looking from above whereas imgCIF is looking from
    below, so the sense of rotation is reversed.

    """

    axes, senses = raw_info["Goniometer axes"]
    # n = len(axes)
    axis_type = ["rotation" for _ in axes]
    equip = ["goniometer" for _ in axes]

    if raw_info.get("kappa") == None:
        kappa_axis = "- - -"
    else:
        kappa_axis = raw_info["kappa"]
    kappa_axis = re.split(r',| ', kappa_axis)[0]

    if raw_info.get("chi") == None:
        chi_axis = "- - -"
    else:
        chi_axis = raw_info["chi"]
    chi_axis = re.split(r',| ', chi_axis)[0]

    # kappa_axis = split(get(raw_info, "kappa","- - -"),[' ',','])[1]
    # chi_axis = split(get(raw_info,"chi","- -"),[' ',','])[1]

    # Construct axis dependency chain
    depends_on = []
    depends_on += axes[1:]
    depends_on.append('none')

    # Remember the principal axis direction
    principal = senses[-1]

    # Create direction vectors
    vector = []
    offset = []
    for (a, d) in zip(axes, senses):
        print(a, d)
        rotfac = 1 if d == principal else -1
        if a.lower() == kappa_axis.lower():
            kv = get_kappa_info(raw_info)
            kv[1] *= rotfac
            vector.append(kv)
        elif a.lower() == chi_axis.lower():
            vector.append(get_chi_info(raw_info))
        else:
            vector.append([i * rotfac for i in [1, 0, 0]])
        print('vec', vector)
        offset.append([0, 0, 0])

    # offset = fill([0,0,0],n)
    return [raw_info["Goniometer axes"][0], axis_type, equip, depends_on, vector, offset]


def make_detector_axes_nxmx(raw_info):

    axes = []
    axis_type = []
    equip = []
    depends_on = []
    vector = []
    offset = []

    print('raw egg', raw_info['gonio_axes_nxmx'])

    for axis in raw_info['gonio_axes_nxmx']:
        axes.append(axis)
        ax_type = raw_info['gonio_axes_nxmx'][axis]['trafo']
        axis_type.append(ax_type)
        #TODO can i do that?
        if ax_type == 'rotation':
            equip.append('goniometer')
        elif ax_type == 'translation':
            equip.append('detector')
        else:
            equip.append('equipment')
        depends_on.append(raw_info['gonio_axes_nxmx'][axis]['dep'])
        vector.append(list(raw_info['gonio_axes_nxmx'][axis]['vec']))
        offset.append(list(raw_info['offsets_nxmx'][axis]['vec']))

    return [axes, axis_type, equip, depends_on, vector, offset]



def make_detector_axes(raw_info):
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

    # Insert assumed axis names and orientations
    axis_id = ["two_theta", "trans", "detx", "dety"]
    axis_type = ["translation" for _ in axis_id]
    axis_type[1] = "rotation"
    equip = ["detector" for _ in axis_id]
    depends_on = ['none', "two_theta", "trans", "detx"]

    # Read necessary information
    print('prinsens', raw_info["Goniometer axes"])
    principal_sense = raw_info["Goniometer axes"][1][-1]
    principal = raw_info["Principal axis orientation"]
    corner = raw_info["Image orientation"]

    # Adjust two theta direction
    rotsense = 1 if raw_info["Two theta axis"] == principal_sense else -1
    vector = [[rotsense, 0, 0]]

    # Detector translation always opposite to beam direction
    vector.append([0, 0, -1])

    # Work out det_x and det_y
    x_d, y_d = determine_detx_dety(principal, principal_sense, corner)
    vector.append(x_d)
    vector.append(y_d)

    # Beam centre is unknown for now
    # offset = [[0, 0, 0], [0, 0, 0], [np.nan, np.nan, 0], [0, 0, 0]]
    offset = [[0, 0, 0], [0, 0, 0], ['none', 'none', 0], [0, 0, 0]]

    return [axis_id, axis_type, equip, depends_on, vector, offset]


def list_depth(list_of_lists):
    if not isinstance(list_of_lists, list):
        return 0
    return max(map(list_depth, list_of_lists), default=0) + 1


def split_vector(vector, basename, imgblock):
    """
        split_vector

    Distribute the vectors in `vector` over CIF data names,
    formatting as integers if they are exact.
    """

    print('veccoi', vector)

    # TODO is the index correct? julia indexes different than python
    def mapping_func(x):
        # print('mpd',x)
        # print('i', i)
        # print('val', x[i - 1])
        if x[i - 1] == 'none':
            return 'none'#pass #TODO x[i] don't get the julia code
        else:
            if round(x[i - 1]) == x[i - 1]:
                # print('balla')
                return f"{x[i - 1]:d}" # x[i]

            else:
                # print('yalla')
                return f"{x[i - 1]:1.5f}" # x[i] @sprintf "%1.5f" x[i]

    for i in range(1, 4):
        imgblock[basename + f"[{i}]"] = map(mapping_func, vector)

        #  do x
        #     if ismissing(x[i]) x[i]
        #     else
        #         if round(x[i]) == x[i]
        #             @sprintf "%d" x[i]
        #         else
        #             @sprintf "%1.5f" x[i]


def get_kappa_info(raw_info):

    # extract kappa information
    axes, senses = raw_info["Goniometer axes"]
    # TODO keepempty=false is ok with python?
    kinfo = raw_info["kappa"].split([' ',','])

    if len(kinfo) == 2:
         kinfo.append("0")

    axname = kinfo[0]
    if not axname in axes:
        # throw(error("Kappa axis $axname not listed as a goniometer axis"))
        # TODO reaise correct here
        raise Exception("Kappa axis $axname not listed as a goniometer axis")

    kapang = float(kinfo[1])
    kappoise = float(kinfo[2])

    # Now calculate direction of axis in X-Z plane assuming
    # rotation of kappa is same as principal axis. There is no
    # Y component as home position is always under
    # beam path.
    up_comp = np.cos(kapang)
    across_comp = np.sin(kapang)
    if kappoise == 0:
        #is under incident beam collimator in X-Z plane
        z_comp = -1 * across_comp
    elif kappoise == 180:
        z_comp = across_comp

    return [up_comp, 0.0, z_comp]


def get_chi_info(raw_info):

    # Extract chi information
    axes, senses = raw_info["Goniometer axes"]
    # TODO keepempty?
    cinfo = raw_info["chi"].split([' ',',']) #,keepempty=false)
    axname = cinfo[0].lower() #lowercase( first( cinfo ) )
    if not axname in axes.lower():
        raise Exception("Chi axis name $axname not listed as a goniometer axis")
    r = -1 * float(cinfo[2])

    # Now turn this into an axis direction
    # At the provided rotation, the chi axis is parallel to z. It is rotated
    # by -chi_rot about omega to bring it to the start position. The sense
    # of omega rotation is always about 1 0 0 by definition
    # chi_sense = indexin([axname], lowercase.(axes))[]
    #TODO correct?
    axes_lower = axes.lower()
    chi_sense = [axes_lower.index(val) for val in [axname]]
    chi_sense = senses[chi_sense]

    # TODO arrays?
    chi_beam_dir = [0, 0, 1] if chi_sense == "a" else [0, 0, -1]
    chi_rot = [
        [1, 0, 0],
        [0, np.cos(r), -np.sin(r)],
        [0, np.sin(r), np.cos(r)]
    ]

    return np.array(chi_rot) * np.array(chi_beam_dir)



def determine_detx_dety(principal_angle, principal_sense, corner):

    # Determine direction of detx (horizontal) and dety (vertical) in
    # imgCIF coordinates.
    # Start with basic value and then flip as necessary

    x_direction = [-1, 0, 0]        # spindle rotates anticlockwise at 0, top_left origin
    y_direction = [0, 1, 0]       #

    #TODO go for np arrays?
    if corner == "top right":
        x_direction = [i * -1 for i in x_direction]

    elif corner == "bottom right":
        x_direction = [i * -1 for i in x_direction]
        y_direction = [i * -1 for i in y_direction]

    elif corner == "bottom left":
        y_direction = [i * -1 for i in y_direction]

    # The direction of the principal axis flips by 180 if the sense
    # changes
    pa = int(principal_angle) if principal_sense == "a" else int(principal_angle) + 180
    if pa >= 360:
        pa = pa - 360

    if pa == 90:
        temp = x_direction
        x_direction = y_direction
        y_direction = [i * -1 for i in temp] #-1*temp
    elif pa == 180:
        x_direction = [i * -1 for i in x_direction]
        y_direction = [i * -1 for i in y_direction]
    elif pa == 270:
        temp = x_direction
        x_direction = [i * -1 for i in y_direction] #-1*y_direction
        y_direction = temp

    return x_direction,y_direction


def describe_facility(raw_info, imgblock):

    base = "_diffrn_source."
    if "Beamline name" in raw_info:
        imgblock[base + "beamline"] = [raw_info["Beamline name"]]
        imgblock[base + "facility"] = [raw_info["Facility name"]]
    else:
        imgblock[base + "make"] = [raw_info["Name of manufacturer"] + "-" + \
            raw_info["Model"]]
        if "Location" in raw_info:
            imgblock[base + "details"] = [f"Located at {raw_info['Location']}"]

def describe_array(raw_info,imgblock):

    """
    describe_detector(raw_info)

    Produce the information required for array_structure_list. Here
    we assume a rectangular detector with x horizontal, y vertical
    """

    hor, vert = get_pixel_sizes(raw_info)
    # array structure list axis
    base = "_array_structure_list_axis."
    imgblock[base + "axis_id"] = ["detx","dety"]
    imgblock[base + "axis_set_id"] = ["1","2"]
    imgblock[base + "displacement"] = [hor/2,vert/2]   #half pixel size
    imgblock[base + "displacement_increment"] = [hor,vert]
    imgblock.CreateLoop([
        base + "axis_id",
        base + "axis_set_id",
        base + "displacement",
        base + "displacement_increment"])

    # array structure list
    base = "_array_structure_list."
    imgblock[base + "array_id"] = ["1","1"]
    imgblock[base + "index"] = [1,2]
    imgblock[base + "axis_set_id"] = ["1","2"]
    imgblock[base + "dimension"] = ['none', 'none']   #number of elements in each direction
    imgblock[base + "direction"] = ["increasing","increasing"]

    if raw_info["Fast direction"] == "horizontal":
        precedence = [1, 2]
    else:
        precedence = [2, 1]

    imgblock[base + "precedence"] = precedence
    imgblock.CreateLoop([
        base + "array_id",
        base + "index",
        base + "axis_set_id",
        base + "dimension",
        base + "direction",
        base + "precedence"])

def describe_detector(raw_info, imgblock):

    base = "_diffrn_detector."
    imgblock[base + "id"]=["1"]
    imgblock[base + "number_of_axes"]=[2]
    imgblock.CreateLoop([base + "id", base + "number_of_axes"])
    #
    base = "_diffrn_detector_axis."
    imgblock[base + "axis_id"] = ["detx","dety"]
    imgblock[base + "detector_id"] = ["1","1"]
    imgblock.CreateLoop([base + "axis_id", base + "detector_id"])


def get_pixel_sizes(raw_info):
    # both = parse.(Float64,split(raw_info["Pixel size"],","))
    both = [float(pxsize) for pxsize in (raw_info["Pixel size"].split(","))]
    if len(both) == 1:
        both.append(both[0])

    print('both', both)

    return both



# if __name__ == '__main__':
#     main()