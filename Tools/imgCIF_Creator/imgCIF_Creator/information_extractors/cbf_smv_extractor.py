import os
import click
import re
import numpy as np
import imgCIF_Creator.information_extractors.user_input_extractor as user_input_extractor
import imgCIF_Creator.assembler.imgCIF_assembler as imgCIF_assembler


# Configuration information
# ROT_AXES = ("chi", "phi", "detector_2theta", "two_theta", "omega",
#             "angle", "start_angle", "kappa")
# TRANS_AXES = ("detector_distance", "dx", "trans", "distance")
# ALWAYS_AXES = ("distance", "two_theta", "detector_2theta")


def extract_cbf_scan(cif_block, directory, output=None, include_archive_directory=False,
                     file_stem=r".*?_", axis=[], new_url=''):

    # Extract and output scan information from minicbf or ADSC files

    # Assumptions:
    # 1. files are named in some form of xxxx_scanno(_)frameno.cbf
    # 2. frames are sequential
    # 3. The first frame will be scan 1 frame 1
    # 4. filenames are the same length
    # 5. axis names are drawn from the lists below
    # 6. "<axis>_increment" signals the increment

    # print(directory, output, include, stem, axis)

    # Process arguments
    # directory = parsed_args["directory"]
    # output = parsed_args["output"] != nothing
    # include_archive_directory = parsed_args["include"] ? splitpath(directory)[end] : ""
    # print("prependdir", include_archive_directory)
    # file_stem = parsed_args["stem"][] == "" ? Regex(".*?_") : Regex(parsed_args["stem"][])
    # file_stem = f'{stem}' if stem is not None else r".*?_"

    # Analyse CBF files
    print('dir', directory)
    print('stem', file_stem)
    scan_info, all_frames = get_scan_info(directory, stem=file_stem)
    # get_scan_info(directory, stem=file_stem)

    # Rename axes
    if len(axis) > 0:
        axis_to_rename = {axis[0].lower() : axis[1].lower()}
        # print("lo axes ", axis_to_rename)
        rename_axes(scan_info, axis_to_rename)
        # print("lo axes renamed", scan_info)



    # print("myarch ", arch, "myloc ", new_url)


    # cif_file = cf.CifFile()
    # cif_block = cf.CifBlock()
    # cif_file['imgCIF'] = cif_block
    # block_id = make_block_id(user_input)
    # cif_block["_audit.block_id"] = [block_id]

    # out_file = parsed_args["output"] != [] ? parsed_args["output"][] : nothing

    if include_archive_directory:
        prepend_dir = os.path.split(directory)[-1]
    else:
        prepend_dir = ""

    # if output is not N:
    imgCIF_assembler.add_scan_info_to_block(scan_info, all_frames, cif_block, new_url, prepend_dir,
                     directory)




def get_scan_info(frame_dir, stem=r".*?_"):

    scan_frame_regex, all_names, first_scan = get_scan_frame_fmt(frame_dir, stem=stem)

    # print('scan_frame_regex', scan_frame_regex)
    # print('all_names', all_names)
    # print('first_scan', first_scan)

    # index them all for speed
    pattern = re.compile(scan_frame_regex)
    all_frames = {}
    for name in all_names:
        matched = pattern.match(name)
        if matched.groupdict().get("scan"):
            all_frames[(matched["scan"], int(matched["frame"]))] = (name)
        else:
            all_frames[("01", int(matched["frame"]))] = (name)


    # find the unique scans
    unique_scans = set(map(lambda x : x[0], all_frames.keys()))
    # unique_scans = unique!(map( x -> x[1], collect( keys( all_frames ) ) ) )
    print(f"{len(unique_scans)} scans found")


    scan_info = {}
    axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
    frame_type = determine_frame_type(os.path.join(frame_dir, all_frames[list(all_frames)[0]]))

    print(f"Discovered {frame_type} files")

    for scan in unique_scans:

        # print('unic scan', scan)

        scan_names = list(filter(lambda x : x[0] == scan, all_frames.keys()))
        # print('scnn', scan_names)
        frames = [x[1] for x in scan_names]


        # Get information for first
        file_name = os.path.join(frame_dir, all_frames[(scan, 1)][0])
        # axes_single_frame are the axis settings for the individual frame a in scan
        # b
        axes_single_frame, scan_ax, scan_incr, exposure, wl = \
            get_frame_info(frame_type, file_name, axes)
        # get_frame_info(frame_type, file_name, axes)
        # print("axes_single_frame ", axes_single_frame, "scanax ", scan_ax, "scaninc ", scan_incr, "exposure ", exposure, "wl ", wl)
        start = axes_single_frame[scan_ax]

        # Get information and last frame
        file_name = os.path.join(frame_dir, all_frames[(scan, len(frames))][0])
        axes_single_frame, _, _, _, _ = get_frame_info(frame_type, file_name, axes)
        finish = axes_single_frame[scan_ax]

        # Check increment and range match
        if not np.isclose(start + scan_incr * (len(frames)-1), finish, atol=1e-6):
            raise Exception(
                f"Scan range does not match increment: \
{start} to {finish}, {len(frames)-1} steps of {scan_incr}")

        scan_details = {"frames" : len(frames),
                   "axis" : scan_ax,
                   "incr" : scan_incr,
                   "time" : exposure,
                   "start" : start,
                   "range" : scan_incr * len(frames),
                   "wavelength" : wl}

        scan_info[scan] = (axes_single_frame, scan_details)

    prune_scan_info(scan_info)

    return scan_info, all_frames



def get_scan_frame_fmt(frame_dir, stem=r".*?_"):
    # Deduce scan/frame naming convention

    first_scan = 1
    file_pattern = re.compile(stem)

    # all_names = readdir(frame_dir)
    all_names = []
    for _, _, files in os.walk(frame_dir):
        for filename in files:
            all_names.append(filename)

    # filter out only .cbf and .img files
    all_names = list(filter(lambda f_name : f_name[-4:] in [".cbf", ".img"], all_names))

    # if !ismissing(stem)
    # if given a file stem filter out only the files that start with the stem/file_pattern
    all_names = list(filter(lambda f_name : file_pattern.match(f_name), all_names))
    all_names.sort()

    # print('all names', all_names)


    # Analyse number of digits between stem and extension: if less than
    # 5, no scan present the default stem matches everything until the _
    test_name = all_names[1]
    stem_len = len(file_pattern.match(test_name).group(0))
    num_digits = len(re.sub("[^0-9]", "", test_name[(stem_len-1):-4]))
    # num_digits = count(r"[0-9]", test_name[stem_len:end-4])

    if num_digits >= 5:
        scan_frame_regex = None

        # Allow first scan to not be scan 1
        for scan in [str(i) for i in range(1,10)]: #corresponds numbers '1'-'9'
            # adding two regex in julia means that both expressions are enclosed in
            # (?:....) which matches everything enclosed in the bracket
            # julia uses PCRE flavour, in python the regular expression for a named
            # capturing group is (?P<name>...) instead of (?<name>...) in PCRE
            # the following regex matches three expressions
            regex = r"(?:" + stem + r")" +\
                r"(?:(?P<scan>[0-9]*" + re.escape(scan) + \
                r")(?P<sep>0|_)(?P<frame>[0-9]+1)(?P<ext>\.cbf|img))"
            # rr = stem * r"(?<scan>[0-9]*$scan)(?<sep>0|_)(?<frame>[0-9]+1)\\.(cbf|img)"
            # print('regex', regex)
            # keep in mind that re.match matche only the beginning of string
            match = re.match(regex, all_names[0])
            # m = match(rr,all_names[1])

            if match:
                # s_pos = m.offsets[1] # this expression matches always from the
                # beginning
                scan_len = len(match.group("scan"))
                frame_len = len(match.group("frame"))

                if match.group("sep") == "0":
                    # if the separator is a 0, include it
                    frame_pos = match.start("sep")
                    frame_len += 1
                else:
                    frame_pos = match.start("frame")

                # scan_frame_regex = stem * Regex("(?<scan>[0-9]{$scan_len})_?(?<frame>[0-9]{$frame_len})")
                scan_frame_regex = r"(?:" + stem + r")" +\
                    r"(?:(?P<scan>[0-9]{" + re.escape(str(scan_len)) + r"})" +\
                    r"_?(?P<frame>[0-9]{" + re.escape(str(frame_len)) + r"}))"

                # print('scrfreg', scan_frame_regex)

                first_scan = scan

                # break if match found
                break

    else:
        # scan_frame_regex = stem * r"(?<frame>[0-9]+)"
        scan_frame_regex = r"(?:" + stem + r")" + r"(?:(?P<frame>[0-9]+))"


    if scan_frame_regex is None:
        raise Exception(f"Cannot find scan/frame naming pattern for {test_name}. Try using the -s option")
        # @error("Cannot find scan/frame naming pattern for $test_name. Try using the -s option")
        # exit()

    # print('matched', re.match(scan_frame_regex, all_names[-1]))
    assert re.match(scan_frame_regex, all_names[-1]), "Regular expression for first \
frame is not matching the last frame."

    return scan_frame_regex, all_names, first_scan





def determine_frame_type(filename):
    """
    Determine the type of frame file: currently SMV (ADSC) and CBF are
    recognised.
    """

    with open(filename, 'rb') as file:
        # read first 512 characters/bytes as byte string
        header = file.read(512)
        # print('the file', readin)
        # TODO ensure this! maybe its also 0x0c for the form feed character
        if b'\f' in header:
            return 'SMV'
        elif b'_array_data' in header:
            return 'CBF'


    # header = read(filename, 512)
    # if 0x0c in header return Val(:SMV) end
    # if occursin("_array_data", String(header)) return Val(:CBF) end

    # return None


def get_frame_info(frame_type, fname, axes):

    if frame_type == "CBF":
        return get_frame_info_CBF(fname, axes)
    elif frame_type == "SMV":
        return get_frame_info_SMV(fname)



def get_frame_info_CBF(fname, axes):
    """
    Return any values found for provided axes. All axes converted to lowercase. This routine
    tries to adapt to all of the crazy stuff stashed in miniCBF headers.
    """

    with open(fname, 'rb') as file:
        # lines = file.readlines()
        lines = []

        for line in file:
            # print('my lane', line)
            if b"--CIF-BINARY-FORMAT-SECTION--" in line:
                break
            else:
                lines.append(line.decode().strip('\n').lower())

    # print('lns', lines)

    #     header = readuntil(fname, "--CIF-BINARY-FORMAT-SECTION--")
    #     lines = lowercase.(split(header, "\n"))

    # print('axes', axes)
    # test = [get_header_values(frame_type, lines, ax) for ax in axes]
    # map(lambda ax : get_header_values(frame_type, lines, ax), axes)
    ax_vals = \
        list(map(lambda ax : (ax.lower(), get_CBF_header_values(lines, ax)),
                 axes))

    ax_vals = filter( lambda x : x[1] != None, ax_vals)
    ax_vals = [convert_units(ax_val) for ax_val in ax_vals]
    # print('ax vals', ax_vals)

    ax_incr = \
        list(map(lambda ax : (ax.lower(), get_CBF_header_values(
            lines, ax + "_increment")), axes))
    #TODO move to get header
    ax_incr = filter( lambda x : x[1] != None, ax_incr)
    ax_incr = [convert_units(ax_val) for ax_val in ax_incr]

    exposure, _ = get_CBF_header_values(lines, "exposure_time")
    wl, _ = get_CBF_header_values(lines, "wavelength")

    # scan_ax = findfirst( x -> !isapprox(x[2], 0, atol = 1e-6), ax_incr)


    # print('axinc', ax_incr)

    # find the first element that is not close
    # scan_ax = next(filter(lambda x: np.isclose(x[1][0], 0, atol=1e-6), ax_incr), None)
    scan_ax = next(filter(lambda x: not np.isclose(x[1], 0, atol=1e-6), ax_incr), None)
    # print('scax', scan_ax)

    # matching_scan_ax = [axis for axis in ax_vals]
    matching_scan_ax = list(filter(lambda ax : scan_ax[0] in ax[0], ax_vals))

    # scan_ax_name = indexin([ax_incr[scan_ax][1]], [x[1] for x in ax_vals])
    if matching_scan_ax == []:
        raise Exception(
            f"Could not match scanned axis {scan_ax} with an axis name.")
    else:
        matching_scan_ax = matching_scan_ax[0][0]

    # print('mtchscax', matching_scan_ax)

    # Get rid of duplicate names
    ax_vals = dict(ax_vals)
    # print('axvals', ax_vals)

    if "distance" in  ax_vals and "detector_distance" in ax_vals:
        del ax_vals["detector_distance"]

    # Some Pilatus headers do not mention omega, just "start_angle"
    # and "angle_increment".
    if "start_angle" in ax_vals and "angle" in ax_vals:
        del ax_vals["start_angle"]

        if matching_scan_ax != "angle":   #we have an actual one
            del ax_vals["angle"]

    return ax_vals, matching_scan_ax, scan_ax[1], exposure, wl



def get_frame_info_SMV(filename):
    # For a single-axis diffractometer currently


    with open(filename, 'rb') as file:
        header = file.read(512)
        lines = header.split("\n").lower()

    ax_vals = [("phi", get_SMV_header_values(lines, "phi"))]
    ax_vals.append(("trans", get_SMV_header_values(lines, "distance")))
    ax_incr = [("phi", get_SMV_header_values(lines, "osc_range"))]
    # TODO notify james about typo
    ax_incr.append(("trans", 0.0))

    exposure = get_SMV_header_values(lines, "time")
    wl = get_SMV_header_values(lines, "wavelength")

    return dict(ax_vals), "phi", ax_incr[0][1], exposure, wl



def get_CBF_header_values(lines, matcher):
    """
    Get the value following the string given in matcher and units if present
    """

    # rr = Regex("$matcher[ =]+")
    # print('callin you')
    pattern = re.compile(re.escape(matcher) + r"[ =]+")
    # matcc = [pattern.search(line) for line in lines]
    # print('atcc', matcc)
    # one_line = filter( x-> !isnothing(match(rr, x)), lines)
    matching_line = list(filter(lambda x : pattern.search(x) is not None, lines))
    # print('mtch lng', matching_line)

    if len(matching_line) != 1:
        return None

    matching_line = matching_line[0]
    #@debug "Extracting from" one_line

    val_unit_regex = re.escape(matcher) + \
        r"[ =]+(?P<val>[A-Za-z0-9+-.]+) +(?P<units>[A-Za-z.]+)"
    val_unit = re.search(val_unit_regex, matching_line)
    val = val_unit["val"].strip()
    units = val_unit["units"].strip()

    # print('v', val, 'u', units)


    # m = match(Regex("$matcher[ =]+(?<val>[A-Za-z0-9+-.]+) +(?<units>[A-Za-z.]+)"), one_line)
    # val = strip(m["val"])
    # units = strip(m["units"])

    #@debug "To get value" val

    return float(val), units


def get_SMV_header_values(lines, matcher):
    """
    Get the value following the string given in matcher and units if present
    """

    pattern = re.compile(r"^" + re.escape(matcher) + r"[ =]+")
    # rr = Regex("^$matcher[ =]+")
    matching_line = list(filter(lambda x : pattern.search(x) is not None, lines))
    # one_line = filter( x-> !isnothing(match(rr, x)), lines)

    if len(matching_line) != 1:
        return None

    matching_line = matching_line[0]
    #@debug "Extracting from" one_line


    val_regex = re.escape(matcher) + \
        r"[ =]+(?P<val>[A-Za-z0-9+-.]+)"
    val_unit = re.search(val_regex, matching_line)
    val = val_unit["val"].strip()

    # m = match(Regex("$matcher[ =]+(?<val>[A-Za-z0-9+-.]+)"), one_line)
    # val = strip(m["val"])

    #@debug "To get value" val

    return float(val), None




def convert_units(ax_val):
    """
    Detect any non-mm translations and convert
    """

    #TODO this does not change the unit string? problem? consistency?
    name, (val, units) = ax_val
    if units == "m":
        val = val * 1000
    elif units == "cm":
        val = val * 10

    return name, val


def rename_axes(scan_info, renaming_scheme):
    """
    Axes in `scan_info` and `renaming_scheme` should already be lowercase.
    """

    if len(renaming_scheme) == 0:
        return

    for content in scan_info.values():
        vals, dets = content

        for axis, value in vals.items():
            if axis in renaming_scheme.keys():
                vals[renaming_scheme[axis]] = value
                del vals[axis]

        if dets["axis"] in renaming_scheme.keys():
            dets["axis"] = renaming_scheme[dets["axis"]]


def prune_scan_info(scan_info):
    """
    Remove reference to any axes that do not change position and are
    essentially zero, but are not in `always_axes`.
    """


    # print('this', scan_info[list(scan_info.keys())[0]] )
    # print('this is scan', scan_info )

    initial_vals, deets = scan_info[list(scan_info.keys())[0]]
    # print('iits', initial_vals, 'dets', deets)
    scan_axis = deets["axis"]
    keep_this = [scan_axis]
    for name, ini_val in initial_vals.items():
        for content in scan_info.values():
            if content[0][name] != ini_val:
                keep_this.append(name)
                break
    # print('wann kepp', keep_this)

    for name, ini_val in initial_vals.items():
        if not (name in imgCIF_assembler.ALWAYS_AXES) and not (name in keep_this) \
            and np.isclose(ini_val, 0, atol=0.001):

            for s in scan_info:
                del(scan_info[s][0], name)

    # print('scan inf', scan_info)


    # return scan_info

# end





@click.command()
@click.option(
    "--location",
    "-l",
    default=None,
    type=str,
    help="Change the _array_data_external_data.uri from the default file path to \
the value to the value given here.",
)

@click.option(
    "--stem",
    "-s",
    default=None,
    type=str,
    help="Constant portion of frame file name. This can help determine the \
scan/frame file naming convention.",
)

@click.option(
    "--include_archive_directory",
    "-i",
    default=False,
    type=bool,
    is_flag=True,
    help="Include directory name as part of frame location info in output. This \
only has an effect if the location is set and has an archive_path \
as .tgz, then _array_data_external_data.archive_path is filled and if selected \
the folder name is prepended to that name",
)

@click.option(
    "--output",
    "-o",
    default=None,
    type=str,
    help="Output file name the result is written to, if it is None then it's\
printed to stdout.",
)

# This changes an axis name and should be done according to the header info
@click.option(
    "--axis",
    "-a",
    default=None,
    type=list,
    help="Change axis name from <cbf> to <new> in output, to match the goniometer\
axis definitions. May be used multiple times for multiple axis renaming."
)

@click.argument("directory", type=str)
def main(directory, output, include, stem, axis):

    extract_cbf_scan(directory, output, include, stem, axis)


if __name__ == '__main__':
    main()
