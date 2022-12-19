import re
from imgCIF_Creator.information_extractors import cbf_smv

# Configuration information
#TODO always?
ROT_AXES = ("chi", "phi", "detector_2theta", "two_theta", "omega",
            "angle", "start_angle", "kappa")
TRANS_AXES = ("detector_distance", "dx", "trans", "distance")
ALWAYS_AXES = ("distance", "two_theta", "detector_2theta")

# #============= Output routines ========================#


def create_imgCIF(filename):

    extractor = cbf_smv.extractor(filename)

    # extractor.get()...
    # generate()...

    pass


def add_scan_info_to_block(scan_info, all_frames, cif_block, new_url,
                           prepend_dir, directory, file_format):
    """
    We identify each frame by its sequence number overall (not just within
    its own scan.)
    """

    # Output CIF fragment
    arch, new_url = get_archive_type(new_url, directory)
    print('new', new_url, arch)
    print('dir', directory)

    # scan_sorted = {key: val for key, val in sorted(scan_info.items(), \
    #     key = lambda ele: ele[0])}
    # print('sort', scan_sorted)
    # scan_info = {f'{idx + 1}' : scan_sorted[scan] for idx, scan in enumerate(scan_info)}


    scan_list = create_scan_list(scan_info)
    # print('scl', scan_list)
    # print('repl', scan_info)
    exposure_time = {s : scan_info[s][1]["time"] for s in scan_info.keys()}
    # print("expin", exp_info)

    # op = isnothing( output_file ) ? stdout : open(output_file, "w")




    # # _diffrn_source block
    # get_facility_info()
    # # describe _axis block
    # get_axes_info()
    # # describe _array_structure_list_axis and _array_structure_list
    # gat_array_info()
    # # describe _diffrn_detector and _diffrn_detector_axis
    # get_detector_info()
    # # generate _diffrn_wavelength block
    # get_wavelength_info()
    # # generate _diffrn_scan_axis block
    # get_scan_settings_info()
    # # generate _diffrn_scan block
    # get_scan_info()
    # # generate _diffrn_scan_frame block
    # get_step_info()
    # # generate _diffrn_data_frame block
    # get_array_info()
    # # generate _array_data block
    # get_ids_info()
    # # generate _array_data_external_data
    # get_external_ids_info()




    # # _diffrn_source block
    # # needs beamline name, facility or manufacturer, model opt location
    # generate_facility(user_input, cif_block)

    # # describe _axis block
    # # needs goniometer axes
    # generate_axes(user_input, cif_block)

    # # describe _array_structure_list_axis and _array_structure_list
    # # needs pixel size, array dimension, fast direction
    # generate_array(user_input, cif_block)

    # # describe _diffrn_detector and _diffrn_detector_axis
    # # needs detector axes,
    # generate_detector(user_input, cif_block)

    # generate _diffrn_wavelength block
    # needs rad_type, wavelength
    generate_wavelength(cif_block, scan_info)

    # generate _diffrn_scan_axis block
    # needs scan info
    generate_scan_settings(cif_block, scan_info)

    # generate _diffrn_scan block
    # needs scan list
    generate_scan_info(cif_block, scan_list)

    # generate _diffrn_scan_frame block
    # needs scan_list, integration time
    generate_step_info(cif_block, scan_list, exposure_time)

    # generate _diffrn_data_frame block
    # needs scan list
    generate_array_info(cif_block, scan_list)

    # generate _array_data block
    # needs scan list
    generate_ids(cif_block, scan_list)

    # generate _array_data_external_data
    # needs uri/zenod?, all frames, scan list, arch, prepend dir, file format
    generate_external_ids(
        cif_block, new_url, all_frames, scan_list, arch, prepend_dir, file_format)









def create_scan_list(scan_info):

    # Create scan list of (scanid, frame_no) where
    # frame_no is just the number of frames

    scans = sorted(scan_info)
    slist = [(s, scan_info[s][1]["frames"]) for s in scans]

    return slist



def generate_wavelength(cifblock, scan_info):

    # print('scaninf', scan_info)
    rad_type = scan_info[list(scan_info.keys())[0]][1].get("rad_type")
    rad_type = "x-ray" if rad_type is None else rad_type
    wl = scan_info[list(scan_info.keys())[0]][1]["wavelength"]
    cifblock["_diffrn_radiation_wavelength.id"] = 1
    cifblock["_diffrn_radiation_wavelength.value"] = wl
    cifblock["_diffrn_radiation.type"] = rad_type



def generate_scan_settings(cif_block, scan_info):
    """

    # Information we have obtained from the cbf file

    #         details = Dict("frames"
    #                        "axis"
    #                        "incr"
    #                        "time"
    #                        "start"
    #                        "range")


    # """

    base = "_diffrn_scan_axis"
    entries = {
        base + '.scan_id' : [],
        base + '.axis_id': [],
        base + '.displacement_start' : [],
        base + '.displacement_increment' : [],
        base + '.displacement_range' : [],
        base + '.angle_start' : [],
        base + '.angle_increment' : [],
        base + '.angle_range' :[],
    }

    for scan in sorted(scan_info):

        # print('scscan', scan_info[scan])

        axes, dets = scan_info[scan]
        # print("keylen", len(axes.keys()))
        #TODO do this for all axes or only the one that changes?
        for axis, val in axes.items():
            step, range = 0, 0
            if axis == dets["axis"]:
                step = dets["incr"]
                range = dets["range"]
                val = dets["start"]

            settings = [val, step, range]
            empty_settings = ['.', '.', '.']

            if axis in TRANS_AXES:
                settings = settings + empty_settings
                # println(op, "SCAN$scan  $a $val $step $range . . .")
            else:
                settings = empty_settings + settings
                # println(op, "SCAN$scan  $a . . . $val $step $range")

            # print('sets', settings)
            # print()

            entries[base + '.scan_id'].append(f"SCAN{scan}")
            entries[base + '.axis_id'].append(axis)
            entries[base + '.displacement_start'].append(settings[0])
            entries[base + '.displacement_increment'].append(settings[1])
            entries[base + '.displacement_range'].append(settings[2])
            entries[base + '.angle_start'].append(settings[3])
            entries[base + '.angle_increment'].append(settings[4])
            entries[base + '.angle_range'].append(settings[5])

    for name, value in entries.items():
        cif_block[name] = value

    cif_block.CreateLoop(list(entries.keys()))



def generate_ids(cif_block, scan_list):
    """

    Array_id is just IMAGE (a single detector module). The
    binary_id is incremented with frame, and all of them
    are located externally so external_id is also incremented
    together with binary_id.
    """

    base = '_array_data'
    entries = {
        base + ".id" : [],
        base + ".binary_id" : [],
        base + ".external_data_id" : [],
        }

    # println(op, header)
    counter = 0
    for _, frames in scan_list:
        for _ in range(1, frames + 1):
            counter += 1
            entries[base + ".id"].append("IMAGE") #TODO is that always IMAGE?
            entries[base + ".binary_id"].append(counter)
            entries[base + ".external_data_id"].append(f'ext{counter}') # TODO put ext here?
            # println(op, "   IMAGE $ctr $ctr")

    for name, value in entries.items():
        cif_block[name] = value

    cif_block.CreateLoop(list(entries.keys()))



def generate_external_ids(cif_block, fulluri, all_frames, scan_list,
                          arch, prepend_dir, file_format):
    """
    `fulluri` is the location of a single archive file. The individual files
    on the local storage are assumed to be at the same locations relative to
    this top-level directory.
    """

    # print('my arch', arch)
    # print('prepdoi', prepend_dir)

    base = '_array_data_external_data'
    entries = {
        base + ".id" : [],
        base + ".format" : [],
        base + ".uri" : [],
        base + ".path" : [],
        base + ".frame" : [],
        }

    if arch is not None:
        entries[base + ".archive_format"] = []
        entries[base + ".archive_path"] = []

    counter = 0
    for scan, frames in scan_list:
        for frame in range(1, frames + 1):
            counter += 1
            frame_name = f"/{all_frames[(scan, frame)][0]}"
            entries[base + ".id"].append(f"ext{counter}")
            entries[base + ".format"].append(file_format)
            entries[base + ".uri"].append(fulluri) #TODO zenodo uri here?
            entries[base + ".path"].append(all_frames[(scan, frame)][1])
            entries[base + ".frame"].append(all_frames[(scan, frame)][2])
            # println(op, "frm$ctr  ELEMENT  IMAGE $(ctr)")

            # A too-clever-by-half way of optionally live constructing a URL
            if arch is not None:
                # print(op, "  $comp $prepend_dir")
                entries[base + ".archive_format"].append(arch)
                entries[base + ".archive_path"].append(prepend_dir + frame_name)


    for name, value in entries.items():
        cif_block[name] = value

    cif_block.CreateLoop(list(entries.keys()))



def generate_scan_info(cif_block, scan_list):
    """
    Fill in the scan information. We number the frames from
    the start
    """
    base = '_diffrn_scan'
    entries = {
        base + ".id" : [],
        base + ".frame_id_start" : [],
        base + ".frame_id_end" : [],
        base + ".frames" : [],
        }

    # println(op, header)
    start_frame = 1
    for scan, frame in scan_list:
        # print(scan, frame)
        end_frame = start_frame + frame - 1
        # println(op, "SCAN$s    frm$start_frame  frm$end_frame  $f")
        entries[base + ".id"].append(f"SCAN{scan}")
        entries[base + ".frame_id_start"].append(f"frm{start_frame}")
        entries[base + ".frame_id_end"].append(f"frm{end_frame}")
        entries[base + ".frames"].append(frame)

        start_frame = end_frame + 1

    for name, value in entries.items():
        cif_block[name] = value

    cif_block.CreateLoop(list(entries.keys()))



def generate_step_info(cif_block, scan_list, scan_times):
    """
    Fill in information about steps
    """

    base = '_diffrn_scan_frame'
    entries = {
        base + ".frame_id" : [],
        base + ".scan_id" : [],
        base + ".frame_number" : [],
        base + ".integration_time" : [],
        }

    # println(op, header)
    start_frame = 1
    counter = 0
    for scan, frames in scan_list:
        for frame_number in range(1, frames + 1):
            counter += 1
            # println(op, "frm$ctr   SCAN$s    $f $(scan_times[s])")
            entries[base + ".frame_id"].append(f"frm{counter}")
            entries[base + ".scan_id"].append(f"SCAN{scan}")
            entries[base + ".frame_number"].append(frame_number)
            entries[base + ".integration_time"].append(scan_times[scan])

    for name, value in entries.items():
        cif_block[name] = value

    cif_block.CreateLoop(list(entries.keys()))



def generate_array_info(cif_block, scan_list):

    base = '_diffrn_data_frame'
    entries = {
        base + ".id" : [],
        base + ".detector_element_id" : [],
        base + ".array_id" : [],
        base + ".binary_id" : [],
        }

    # println(op, header)
    counter = 0
    for _, frames in scan_list:
        for _ in range(1, frames + 1):
            counter += 1
            entries[base + ".id"].append(f"frm{counter}")
            entries[base + ".detector_element_id"].append(f"ELEMENT")
            #TODO element or 1? not present in hdf5 mapper
            entries[base + ".array_id"].append("IMAGE") # is 1 in hdf5 mapper
            entries[base + ".binary_id"].append(counter)
            # println(op, "frm$ctr  ELEMENT  IMAGE $(ctr)")

    for name, value in entries.items():
        cif_block[name] = value

    cif_block.CreateLoop(list(entries.keys()))




def get_archive_type(new_url, frame_dir):

    archive = None

    if new_url == []:
        frame_dir = frame_dir[1:] if frame_dir.startswith('/') else frame_dir
        new_url = "file://" + frame_dir

    else:
        archives = {"TGZ" : r".*((\.tgz\Z)|(\.tar\.gz\Z|))",
                    "TBZ" : r".*((\.tbz\Z)|(\.tar\.bz2\Z|))",
                    "ZIP" : r".*(\.zip\Z)"}

        for archive, regex in archives.items():
            if re.match(regex, new_url):
                return archive, new_url

    return archive, new_url


def generate_facility(raw_info, imgblock):

    base = "_diffrn_source."
    if "Beamline name" in raw_info or "Facility name" in raw_info:
        imgblock[base + "beamline"] = [raw_info["Beamline name"]]
        imgblock[base + "facility"] = [raw_info["Facility name"]]
    else:
        imgblock[base + "make"] = [raw_info["Name of manufacturer"] + "-" + \
            raw_info["Model"]]
        if "Location" in raw_info:
            imgblock[base + "details"] = [f"Located at {raw_info['Location']}"]

def generate_array(raw_info,imgblock):

    """
    describe_detector(raw_info)

    Produce the information required for array_structure_list. Here
    we assume a rectangular detector with x horizontal, y vertical
    """

    hor, vert = get_pixel_sizes(raw_info)
    # array structure list axis
    base = "_array_structure_list_axis."
    # TODO always detx dety?
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
    imgblock[base + "index"] = [1, 2]
    imgblock[base + "axis_set_id"] = ["1","2"]
    imgblock[base + "dimension"] = raw_info['Array dimension']# ['none', 'none']   #number of elements in each direction
    imgblock[base + "direction"] = ["increasing", "increasing"]

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

def generate_detector(raw_info, imgblock):

    base = "_diffrn_detector."
    imgblock[base + "id"]=["1"]
    imgblock[base + "number_of_axes"] = raw_info.get('detector_axes', [2])
    imgblock.CreateLoop([base + "id", base + "number_of_axes"])
    #
    base = "_diffrn_detector_axis."
    # print('olla', imgblock["_diffrn_detector.number_of_axes"])
    if imgblock["_diffrn_detector.number_of_axes"][0] == 2:
        imgblock[base + "axis_id"] = ["detx", "dety"]
        imgblock[base + "detector_id"] = ["1", "1"]
    elif imgblock["_diffrn_detector.number_of_axes"][0] == 1:
        imgblock[base + "axis_id"] = ["trans"]
        imgblock[base + "detector_id"] = ["1"]

    imgblock.CreateLoop([base + "axis_id", base + "detector_id"])


def generate_axes(raw_info, imgblock):
    """
    describe_axes(raw_info,imgblock)

    Create the goniometer axes corresponding to the data in `raw_info`,
    placing the result in CIF block `imgblock`.
    """

    if 'gonio_axes_nxmx' in raw_info:
        gon_axes = make_axes_nxmx(raw_info)
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
