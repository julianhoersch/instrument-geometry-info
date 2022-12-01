import re

# Configuration information
ROT_AXES = ("chi", "phi", "detector_2theta", "two_theta", "omega",
            "angle", "start_angle", "kappa")
TRANS_AXES = ("detector_distance", "dx", "trans", "distance")
ALWAYS_AXES = ("distance", "two_theta", "detector_2theta")

# #============= Output routines ========================#


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
    print('scl', scan_list)
    print('repl', scan_info)
    exposure_time = {s : scan_info[s][1]["time"] for s in scan_info.keys()}
    # print("expin", exp_info)

    # op = isnothing( output_file ) ? stdout : open(output_file, "w")

    # generate _diffrn_wavelength block
    generate_wavelength(cif_block, scan_info)
    # generate _diffrn_scan_axis block
    generate_scan_settings(cif_block, scan_info)
    # generate _diffrn_scan block
    generate_scan_info(cif_block, scan_list)
    # generate _diffrn_scan_frame block
    generate_step_info(cif_block, scan_list, exposure_time)
    # generate _diffrn_data_frame block
    generate_array_info(cif_block, scan_list)
    # generate _array_data block
    generate_ids(cif_block, scan_list)
    # generate _array_data_external_data
    generate_external_ids(
        cif_block, new_url, all_frames, scan_list, arch, prepend_dir, file_format)



def create_scan_list(scan_info):

    # Create scan list of (scanid, frame_no) where
    # frame_no is just the number of frames

    scans = sorted(scan_info)
    slist = [(s, scan_info[s][1]["frames"]) for s in scans]

    return slist



def generate_wavelength(cifblock, scan_info, rad_type="xray"):

    # print('scaninf', scan_info)
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
            entries[base + ".id"].append("IMAGE")
            entries[base + ".binary_id"].append(counter)
            entries[base + ".external_data_id"].append(counter)
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
            entries[base + ".id"].append(counter)
            entries[base + ".format"].append(file_format)
            entries[base + ".uri"].append(fulluri)
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
            entries[base + ".array_id"].append("IMAGE")
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
