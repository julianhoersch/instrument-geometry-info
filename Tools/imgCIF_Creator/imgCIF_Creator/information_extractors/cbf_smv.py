import os
import re
import sys
import numpy as np
from imgCIF_Creator.output_assembler import imgCIF_assembler
from . import extractor_interface
from . import full_cbf


class extractor(extractor_interface.ExtractorInterface):

    def __init__(self, directory, stem) -> None:
        super().__init__()

        self._unique_scans, self.all_frames = \
            self._get_scans_and_frames(directory, stem=stem)

        # print('u sca', self._unique_scans)
        # print('all fr', self.all_frames)
        # retrieve mini header info
        # only mini header info contains scan details regarding increment etc
        # otherwise it would be duplicated
        self._scan_info_mini_header = \
            self._get_scan_info_mini_header(directory, self._unique_scans, self.all_frames)

        # TODO optional axis renaming
        # print('mini info:', self._scan_info_mini_header)

        # retrieve full header info
        self._scan_info_full_header = \
            self._get_scan_info_full_header(directory, self._unique_scans, self.all_frames)
        # print('full info', self._scan_info_full_header)

        self._first_scan = sorted(self._scan_info_mini_header.keys())[0]
        # print('frst sc', self._first_scan)

        self.first_mini_header = \
            self._scan_info_mini_header[self._first_scan][1]['mini_header']
        # print('full info first:', self.first_mini_header)

        full_header_empty = self._scan_info_full_header[self._first_scan].keys() == []

        if full_header_empty:
            self._data_name = None
            self._full_header_dict = {}
        else:
            self._data_name = \
                list(self._scan_info_full_header[self._first_scan].keys())[0]
            self._full_header_dict = \
                self._scan_info_full_header[self._first_scan][self._data_name]

        # print('full header dict', self._full_header_dict)


    def get_all_frames(self):

        # all frames format?
        return self.all_frames


    def get_uncategorized_info(self):

        doi = self._full_header_dict.get('_database.dataset_doi')
        overload = self._full_header_dict.get('_array_intensities.overload')

        # get doi also from mini header?
        # if doi is None:
        #     self._scan_info_mini_header[self._first_scan][1]['x_pixel_size']

        return {'doi' : doi,
                'overload' : overload,
                }


    def get_source_info(self):

        facility = None
        beamline = None

        # TODO can it appear in the mini header
        # print(self._scan_info_full_header[self._first_scan][self._data_name].keys())
        block_ids = ['_diffrn_source.diffrn_id', '_diffrn.id']
        for b_id in block_ids:
            #print('bid', b_id)
            block_id = self._full_header_dict.get(b_id)
            if type(block_id) == list:
                block_id = block_id[0]
            #print('blockz', block_id)
            if block_id is not None:
                facility_short, beamline = block_id.split('_')

        # check beamline is the same?
        source_string = self._full_header_dict.get('_diffrn_source.type')
        if type(source_string) == list:
                source_string = source_string[0]
        if source_string is not None:
            if 'beamline' in source_string:
                splitter = 'beamline'
            elif 'Beamline' in source_string:
                splitter = 'Beamline'
            facility, beamline_source = source_string.split(splitter)
            if beamline is None:
                beamline = beamline_source

            # print('fac', facility, 'bl', beamline)

        # TODO match facility and manufacturer
        manufacturer = None
        model = None
        location = None

        make_string = self._full_header_dict.get('_diffrn_source.make')
        if type(make_string) == list:
                make_string = make_string[0]
        if make_string is not None:
            manufacturer, model = make_string.split('-')

        # is this meant like general details?
        location_string = self._full_header_dict.get('_diffrn_source.details')
        if type(location_string) == list:
            location_string = location_string[0]

        if location_string is not None:
            location = location_string

        source_info = {'beamline' : beamline,
                       'facility' : facility,
                       'manufacturer' : manufacturer,
                       'model' : model,
                       'location' : location}

        return source_info


    def get_axes_info(self):

        axes = self._full_header_dict.get('_axis.id')
        # print('axe', axes)
        axis_type = self._full_header_dict.get('_axis.type')
        equip = self._full_header_dict.get('_axis.equipment')
        depends_on = self._full_header_dict.get('_axis.depends_on')

        vector = []
        offset = []
        if axes is not None:
            for idx, _ in enumerate(axes):
                sub_vector = []
                sub_offset = []

                # possible with pycifrw but uninituitive
                for i in [1, 2, 3]:
                    entry = self._full_header_dict.get(f'_axis.vector[{i}]')
                    entry = entry if entry is None else entry[idx]
                    sub_vector.append(entry)
                    entry = self._full_header_dict.get(f'_axis.offset[{i}]')
                    entry = entry if entry is None else entry[idx]
                    sub_offset.append(entry)

                vector.append(sub_vector)
                offset.append(sub_offset)

        # print('vec', vector)
        # print('off', offset)
        axes_info = {'axes' : axes,
                     'axis_type' : axis_type,
                     'equip' : equip,
                     'depends_on' : depends_on,
                     'vector' : vector,
                     'offset' : offset}

        return axes_info


    def get_array_info(self):

        #TODO what about the mini header information??
        # get pixel size, array dimension and fast direction
        # array structure list axis
        base = "_array_structure_list_axis."
        # TODO always detx dety?
        axis_id = self._full_header_dict.get(base + "axis_id")
        # TODO is the axis_set_id the same for structure list and list_axis?
        axis_set_id = self._full_header_dict.get(base + "axis_set_id")
        pixel_size = self._full_header_dict.get(base + "displacement_increment")

        if pixel_size is None:
            x_px = self._scan_info_mini_header[self._first_scan][1]['x_pixel_size']
            y_px = self._scan_info_mini_header[self._first_scan][1]['y_pixel_size']
            pixel_size = [x_px, y_px]

        # array structure list
        base = "_array_structure_list."
        array_id = self._full_header_dict.get(base + "array_id")
        array_index = self._full_header_dict.get(base + "index")
        array_dimension = self._full_header_dict.get(base + "dimension")
        # ['none', 'none']   #number of elements in each direction
        array_direction = self._full_header_dict.get(base + "direction")
        array_precedence = self._full_header_dict.get(base + "precedence")

        array_info = {
            'axis_id' : axis_id,
            'axis_set_id': axis_set_id,
            'pixel_size' : pixel_size,
            'array_id' : array_id,
            'array_index' : array_index,
            'array_dimension' : array_dimension,
            'array_direction' : array_direction,
            'array_precedence' : array_precedence,
        }
        return array_info


    def get_detector_info(self):

        #TODO multiple detectors are not supportet (yet?)
        detector_id = \
            self._full_header_dict.get('_diffrn_detector.id')
        number_of_axes = \
            self._full_header_dict.get('_diffrn_detector_axis.number_of_axes')

        axis_id = self._full_header_dict.get('_diffrn_detector_axis.axis_id')
        detector_axis_id = self._full_header_dict.get('_diffrn_detector_axis.detector_id')

        detector_info = {
            'detector_id' : detector_id,
            'number_of_axes' : number_of_axes,
            'axis_id' : axis_id,
            'detector_axis_id' : detector_axis_id
        }
        return detector_info


    def get_wavelength_info(self):

        # TODO not creating a list of wl id's here (yet)
        # TODO could infer rad type from wl...
        rad_type = \
            self._full_header_dict.get('_diffrn_radiation.type')

        #TODO units correct?
        # get from full header
        wavelength = None
        base = '_diffrn_radiation_wavelength.'
        for identifier in ['wavelength', 'value']:
            if wavelength is None:
                wavelength = self._full_header_dict.get(base + identifier)

        # get from mini header
        # TODO can this differ from scan to scan? it better shouldnt...\
        # assert that it is the same wl?
        if wavelength is None:
            wavelength = self._scan_info_mini_header[self._first_scan][1].get('wavelength')

        return {'rad_type' : rad_type,
                'wavelength' : wavelength}


    def get_scan_settings_info(self):

        # is already in a dict format
        return self._scan_info_mini_header


    def _get_scans_and_frames(self, frame_dir, stem=r".*?_"):

        # Extract and output scan information from minicbf or ADSC files

        # Assumptions:
        # 1. files are named in some form of xxxx_scanno(_)frameno.cbf
        # 2. frames are sequential
        # 3. The first frame will be scan 1 frame 1
        # 4. filenames are the same length
        # 5. axis names are drawn from the lists below
        # 6. "<axis>_increment" signals the increment

        scan_frame_regex, all_names, first_scan = self._get_scan_frame_fmt(frame_dir, stem=stem)

        # print('scan_frame_regex', scan_frame_regex)
        # print('all_names', all_names)
        # print('first_scan', first_scan)

        # index them all for speed
        pattern = re.compile(scan_frame_regex)
        all_frames = {}
        for name in all_names:
            matched = pattern.match(name)
            if matched.groupdict().get("scan"):
                all_frames[(matched["scan"], int(matched["frame"]))] = \
                    {'filename' : name}
            else:
                all_frames[("01", int(matched["frame"]))] = {'filename' : name}

        # find the unique scans
        unique_scans = set(map(lambda x : x[0], all_frames.keys()))
        # unique_scans = unique!(map( x -> x[1], collect( keys( all_frames ) ) ) )
        print(f"{len(unique_scans)} scans found")

        return unique_scans, all_frames


    def _get_scan_info_full_header(self, frame_dir, unique_scans, all_frames):

        scan_info = {}
        for scan in unique_scans:
            # scan_names = list(filter(lambda x : x[0] == scan, all_frames.keys()))
            # Get information for first
            file_name = os.path.join(frame_dir, all_frames[(scan, 1)]['filename'])
            full_header = full_cbf.extract_full_cbf_header_information(file_name)

            scan_info[scan] = full_header#[list(full_header.keys())[0]]

        return scan_info


    def _get_scan_info_mini_header(self, frame_dir, unique_scans, all_frames):


        scan_info = {}
        axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
        frame_type = self._determine_frame_type(
            os.path.join(frame_dir, all_frames[list(all_frames)[0]]['filename']))

        print(f"Discovered {frame_type} files")

        for scan in unique_scans:

            # print('unic scan', scan)
            scan_names = list(filter(lambda x : x[0] == scan, all_frames.keys()))
            # print('scnn', scan_names)
            frames = [x[1] for x in scan_names]

            # Get information for first
            # print('frdir', frame_dir)
            # print('allfrs', all_frames[(scan, 1)])
            file_name = os.path.join(frame_dir, all_frames[(scan, 1)]['filename'])
            # print('myfilename', file_name)
            # axes_single_frame are the axis settings for the individual frame a in scan
            # b
            mini_header = self._get_mini_header(file_name, frame_type)
            axes_single_frame, scan_ax, scan_incr, exposure, wl,  = \
                self._get_frame_info(mini_header, frame_type, axes)
            #TODO units?
            x_pixel_size, _, y_pixel_size, _ = self._get_pixel_sizes(mini_header)
            # print('pxsz', x_pixel_size, y_pixel_size)

            # self._get_frame_info(frame_type, file_name, axes)
            # print("axes_single_frame ", axes_single_frame, "scanax ", scan_ax, "scaninc ", scan_incr, "exposure ", exposure, "wl ", wl)
            start = axes_single_frame[scan_ax]

            # Get information and last frame
            file_name = \
                os.path.join(frame_dir, all_frames[(scan, len(frames))]['filename'])
            mini_header = self._get_mini_header(file_name, frame_type)
            # print('mini headi', mini_header)
            axes_single_frame, _, _, _, _ = self._get_frame_info(mini_header, frame_type, axes)
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
                            "wavelength" : wl,
                            "x_pixel_size" : x_pixel_size,
                            "y_pixel_size" : y_pixel_size,
                            "mini_header" : mini_header
                            }

            scan_info[scan] = (axes_single_frame, scan_details)

        self._prune_scan_info(scan_info)

        return scan_info


    def _get_scan_frame_fmt(self, frame_dir, stem=r".*?_"):
        # Deduce scan/frame naming convention

        first_scan = 1
        file_pattern = re.compile(stem)

        # all_names = readdir(frame_dir)
        all_names = []
        for _, _, files in os.walk(frame_dir):
            for filename in files:
                all_names.append(filename)

        # filter out only .cbf and .img files
        # print('all_names', all_names)
        all_names = list(filter(lambda f_name : f_name[-4:] in [".cbf", ".img"], all_names))

        # if !ismissing(stem)
        # if given a file stem filter out only the files that start with the stem/file_pattern
        all_names = list(filter(lambda f_name : file_pattern.match(f_name), all_names))
        all_names.sort()

        # print('all names', all_names)
        # Analyse number of digits between stem and extension: if less than
        # 5, no scan present the default stem matches everything until the _
        # print('all_names', all_names)
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
            print(f"Cannot find a scan/frame naming pattern for {test_name}. \
Try to provide the constant stem of the file name using the -s option.\n")
            sys.exit()

        # print('matched', re.match(scan_frame_regex, all_names[-1]))
        assert re.match(scan_frame_regex, all_names[-1]), "Regular expression for first \
    frame is not matching the last frame."

        return scan_frame_regex, all_names, first_scan


    def _determine_frame_type(self, filename):
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


    def _get_frame_info(self, mini_header, frame_type, axes):

        if frame_type == "CBF":
            return self._get_frame_info_CBF(mini_header, axes)
        elif frame_type == "SMV":
            return self._get_frame_info_SMV(mini_header)


    def _get_mini_header(self, fname, frame_type):

        if frame_type == "CBF":
            return self._get_CBF_header(fname)
        elif frame_type == "SMV":
            return self._get_SMV_header(fname)


    def _get_CBF_header(self, fname):

        # TODO do I need to read this as byte string?
        with open(fname, 'rb') as file:
            cbf_header = []

            # TODO same for smv
            found_header = False
            within_section = False
            for line in file:
                # print('within', within_section)
                # print('my lane', line)
                if b'_array_data.header_contents' in line:
                    found_header = True
                    # print('head', found_header)
                # TODO does that catch all CRLF?
                elif (line == b';\n') or (line == b';\r\n'):
                    within_section = not within_section
                elif b"-BINARY-FORMAT-SECTION-" in line:
                    break

                if found_header and within_section and line.startswith(b'#'):
                    cbf_header.append(line.decode('utf-8').lower().strip('\n').strip('\r'))

        # print('lns', cbf_header)
        return cbf_header


    def _get_frame_info_CBF(self, cbf_header, axes):
        """
        Return any values found for provided axes. All axes converted to lowercase.
        This routine tries to adapt to all of the crazy stuff stashed in miniCBF
        headers.
        """

        ax_vals = \
            list(map(lambda ax : (ax.lower(), self._get_CBF_header_values(cbf_header, ax)),
                    axes))

        ax_vals = list(filter( lambda x : x[1][0] != None, ax_vals))
        ax_vals = [self._convert_units(ax_val) for ax_val in ax_vals]
        # print('ax vals', ax_vals)

        ax_incr = \
            list(map(lambda ax : (ax.lower(), self._get_CBF_header_values(
                cbf_header, ax + "_increment")), axes))

        #TODO move to get header
        # print('axinc 1', ax_incr)
        ax_incr = list(filter( lambda x : x[1][0] != None, ax_incr))
        # ax_incr = [incr for incr in ax_incr if incr[1][0] is not None]
        # print('tp', type(ax_incr))
        # print('axinc 2', ax_incr)
        ax_incr = [self._convert_units(incr) for incr in ax_incr]
        # print('axinc 3', ax_incr)

        # TODO unit conversion?
        exposure, _ = self._get_CBF_header_values(cbf_header, "exposure_time")
        wl, _ = self._get_CBF_header_values(cbf_header, "wavelength")

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


    def _get_frame_info_SMV(self, filename):
        # For a single-axis diffractometer currently

        smv_header = self._get_SMV_header(filename)

        ax_vals = [("phi", self._get_SMV_header_values(smv_header, "phi"))]
        ax_vals.append(("trans", self._get_SMV_header_values(smv_header, "distance")))
        ax_incr = [("phi", self._get_SMV_header_values(smv_header, "osc_range"))]
        ax_incr.append(("trans", 0.0))

        exposure = self._get_SMV_header_values(smv_header, "time")
        wl = self._get_SMV_header_values(smv_header, "wavelength")

        return dict(ax_vals), "phi", ax_incr[0][1], exposure, wl


    def _get_CBF_header_values(slef, lines, matcher):
        """
        Get the value following the string given in matcher and units if present
        """

        pattern = re.compile(re.escape(matcher) + r"[ =]+")

        # one_line = filter( x-> !isnothing(match(rr, x)), lines)
        matching_lines = list(filter(lambda x : pattern.search(x) is not None, lines))
        # print('mtch lng', matching_lines)

        # TODO is it now ensured that there are not more than one lines matching?
        if len(matching_lines) < 1:
            # print('didn;t mathc')
            return None, None

        val_unit_regex = re.escape(matcher) + \
            r"[ =]+(?P<val>[A-Za-z0-9+-.]+) +(?P<units>[A-Za-z.]+)"

        val_unit = [re.search(val_unit_regex, matching_line)
            for matching_line in matching_lines]

        # print('valuns', val_unit)
        val_unit = val_unit[0]
        val = val_unit["val"].strip()
        units = val_unit["units"].strip()

        # print('v', val, 'u', units)
        # m = match(Regex("$matcher[ =]+(?<val>[A-Za-z0-9+-.]+) +(?<units>[A-Za-z.]+)"), one_line)

        return float(val), units


    def _get_pixel_sizes(self, lines):

        matcher = 'pixel_size'

        pattern = re.compile(re.escape(matcher) + r"[ =]+")
        matching_lines = list(filter(lambda x : pattern.search(x) is not None, lines))
        if len(matching_lines) < 1:
            return None, None, None, None

        # print('matched', matching_lines)
        dim = r"([A-Za-z0-9+-.]+) +([A-Za-z.]+)"
        val_unit_regex = re.escape(matcher) + r'[ =]+' + dim + r' [A-Za-z] ' + dim
        # print('reg', val_unit_regex)
        val_unit = [re.search(val_unit_regex, matching_line)
            for matching_line in matching_lines]

        # print('valuns', val_unit)
        val_unit = val_unit[0]
        pixel_sizes = [group.strip() for group in val_unit.groups()]
        # print(pixel_sizes)

        return float(pixel_sizes[0]), pixel_sizes[1], float(pixel_sizes[2]), pixel_sizes[3]


    def _get_SMV_header(self, fname):

        # TODO look at smv files
        with open(fname, 'rb') as file:
            header = file.read(512)
            smv_header = header.split("\n").lower()

        return smv_header


    def _get_SMV_header_values(self, lines, matcher):
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

        val_regex = re.escape(matcher) + \
            r"[ =]+(?P<val>[A-Za-z0-9+-.]+)"
        val_unit = re.search(val_regex, matching_line)
        val = val_unit["val"].strip()

        # m = match(Regex("$matcher[ =]+(?<val>[A-Za-z0-9+-.]+)"), one_line)
        # val = strip(m["val"])
        #@debug "To get value" val

        return float(val), None


    def _convert_units(self, ax_val):
        """
        Detect any non-mm translations and convert
        """

        #TODO this does not change the unit string? problem? consistency?
        name, (val, units) = ax_val
        if units == "m":
            val = val * 1000
        elif units == "cm":
            val = val * 10

        # print('conved', ax_val, 'to', name, val)

        return name, val


    def _rename_axes(self, scan_info, renaming_scheme):
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


    def _prune_scan_info(self, scan_info):
        """
        Remove reference to any axes that do not change position and are
        essentially zero, but are not in `always_axes`.
        """

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
