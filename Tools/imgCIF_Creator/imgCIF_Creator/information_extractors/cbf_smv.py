import os
import re
import sys
import numpy as np
from imgCIF_Creator.output_creator import imgCIF_creator
from . import extractor_interface
from . import full_cbf


class extractor(extractor_interface.ExtractorInterface):
    """See also documentation of the init method.

    Args:
        extractor_interface (class): the interface that must be implemented by
            every extractor class
    """

    def __init__(self, directory, stem) -> None:
        """This extractor allows to extract the scan and setup information from
        cbf and smv files. When an instance of the extractor is initialized
        then the information is attempted to be extracted and stored in class
        attributes. The extractor provides public methods to make this information
        accessible.

        Args:
            directory (str): the directory where the extractor tries to extract
                the information
            stem (str): constant portion of a frame file name to help determine the
                scan/frame file naming convention
        """

        self._unique_scans, self.all_frames = \
            self._get_scans_and_frames(directory, stem=stem)

        # print('u sca', self._unique_scans)
        # print('all fr', self.all_frames)
        # retrieve mini header info
        # only mini header info contains scan details regarding increment etc
        # otherwise it would be duplicated
        self._scan_info_mini_header = \
            self._get_scan_info_mini_header(directory, self._unique_scans, self.all_frames)

        # TODO optional axis renaming?
        # print('mini info:', self._scan_info_mini_header)

        # retrieve full header info
        self._scan_info_full_header = \
            self._get_info_full_header(directory, self._unique_scans, self.all_frames)
        # print('full info', self._scan_info_full_header)

        self._first_scan = sorted(self._scan_info_mini_header.keys())[0]
        # print('frst sc', self._first_scan)

        self.first_mini_header = \
            self._scan_info_mini_header[self._first_scan][1]['mini_header']
        # print('full info first:', self.first_mini_header)

        full_header_is_empty = self._scan_info_full_header[self._first_scan].keys() == []

        if full_header_is_empty:
            self._data_name = None
            self._full_header_dict = {}
        else:
            self._data_name = \
                list(self._scan_info_full_header[self._first_scan].keys())[0]
            self._full_header_dict = \
                self._scan_info_full_header[self._first_scan][self._data_name]

        # print('full header dict', self._full_header_dict)


    def get_all_frames(self):
        """Get a dictionary containing an entry for each frame with corresponding
        file name.
        Format:
        {('scan name', frame number): {'filename': 'ciclohexano3_010001.cbf'}, ...}

        Returns:
            dict: a dictionary containing all frames and their corresponding file
        """

        return self.all_frames


    def get_uncategorized_info(self):
        """Return the information that was found about the doi and the array
        intensities overload.

        Returns:
            dict: a dictionary containing the doi and the array intensities overload
        """

        doi = self._full_header_dict.get('_database.dataset_doi')
        overload = self._full_header_dict.get('_array_intensities.overload')

        # get doi also from mini header?
        return {'doi' : doi,
                'overload' : overload,
                }


    def get_source_info(self):
        """Return the information about the facility and beamline or the instrument,
        model and location. Cif block: _diffrn_source

        Returns:
            dict: a dictionary containing the information about the source
        """

        facility = None
        beamline = None

        # TODO can it appear in the mini header?
        # print(self._scan_info_full_header[self._first_scan][self._data_name].keys())
        block_ids = ['_diffrn_source.diffrn_id', '_diffrn.id']
        for b_id in block_ids:
            block_id = self._full_header_dict.get(b_id)
            if type(block_id) == list:
                block_id = block_id[0]
            if block_id is not None:
                beamline = block_id.split('_')[-1]

        # check beamline is the same?
        source_string = self._full_header_dict.get('_diffrn_source.type')
        if type(source_string) == list:
                source_string = source_string[0]
        if source_string is not None:
            if 'beamline' in source_string:
                splitter = 'beamline'
            elif 'Beamline' in source_string:
                splitter = 'Beamline'
            else:
                splitter = None
            if splitter is not None:
                facility, beamline_source = source_string.split(splitter)
                if beamline is None:
                    beamline = beamline_source

            # print('fac', facility, 'bl', beamline)

        # TODO match facility and manufacturer?
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
        """Return the information about the axes settings. Cif block: _axis

        Returns:
            dict: a dictionary containing the information about the axes settings
        """

        axes = self._full_header_dict.get('_axis.id')
        # print('axes', axes)
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
                # assign the vector and offset to the corresponding axis id
                # from the full header we obtain only the column with all axes
                for i in [1, 2, 3]:
                    entry = self._full_header_dict.get(f'_axis.vector[{i}]')
                    entry = entry[idx] if entry is not None else entry
                    sub_vector.append(entry)
                    entry = self._full_header_dict.get(f'_axis.offset[{i}]')
                    entry = entry[idx] if entry is not None else entry
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
        """Return the information about the array. Cif block: _array_structure_list_axis
        and _array_structure_list

        Returns:
            dict: a dictionary containing the information about the array
        """

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
        """Return the information about the detector. Cif block: _diffrn_detector
        and _diffrn_detector_axis.

        Returns:
            dict: a dictionary containing the information about the detector
        """

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


    def get_radiation_info(self):
        """Return the information about the wavelength an type of radiation.
        Cif block: _diffrn_radiation and _diffrn_radiation_wavelength

        Returns:
           dict: a dictionary containing the information about the radiation
        """

        # TODO not creating a list of wl id's here (yet)
        # TODO could infer rad type from wl...
        rad_type = \
            self._full_header_dict.get('_diffrn_radiation.type')

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
        """Return the information about the scans, this is a dictionary containing
        the starting point settings of the axes and the details of each scan.

        For example for scan '08':
        {'08': ({'chi': -60.991, 'phi': 110.0, 'detector_2theta': -12.4,
        'omega': -18.679, 'distance': 40.0}, {'frames': 12, 'axis': 'omega',
        'incr': 2.0, 'time': 1800.0, 'start': -40.679, 'range': 24.0,
        'wavelength': 0.560834, 'x_pixel_size': 0.172, 'y_pixel_size': 0.172,
        'mini_header': ['# detector: pilatus100k',.... ],
        ...}

        Returns:
            dict: a dictionary containing the information about the scans
        """

        return self._scan_info_mini_header


    def _get_scans_and_frames(self, frame_dir, stem=r".*?_"):
        """Extract scan information from minicbf or ADSC files.

        Assumptions:
        1. files are named in some form of xxxx_scanno(_)frameno.cbf
        2. frames are sequential
        3. The first frame will be scan 1 frame 1
        4. filenames are the same length
        5. "<axis>_increment" signals the increment

        Example all frames:
        {('scan name', frame number): {'filename': 'ciclohexano3_010001.cbf'}, ...}
        Example unique scans:
        {'06', '03', '07', '01', '04', '08', '02', '05'}

        Args:
            frame_dir (str): the directory where the frames are located
            stem (regexp, optional): The constant portion of the scan/frame
                naming convention. Defaults to r".*?_".

        Returns:
            tuple: the unique scan names in a set and the scan name / scan frame
                file mapping
        """

        scan_frame_regex, all_names = self._get_scan_frame_fmt(frame_dir, stem=stem)
        # print('scan_frame_regex', scan_frame_regex)
        # print('all_names', all_names)

        pattern = re.compile(scan_frame_regex)
        all_frames = {}
        # if we can't find a matched scan, we assume that there is only one called "01"
        for name in all_names:
            matched = pattern.match(name)
            if matched.groupdict().get("scan"):
                all_frames[(matched["scan"], int(matched["frame"]))] = \
                    {'filename' : name}
            else:
                all_frames[("01", int(matched["frame"]))] = {'filename' : name}

        # find the unique scans
        unique_scans = set(map(lambda x : x[0], all_frames.keys()))
        print(f"{len(unique_scans)} scan(s) found")

        return unique_scans, all_frames


    def _get_info_full_header(self, frame_dir, unique_scans, all_frames):
        """Return the information that can be found in the complete cbf file.

        Args:
            frame_dir (str): the directory containing the files with frames
            unique_scans (set): a set of unique scans
            all_frames (_type_): a dictionary containing the link between frames
                and filenames.

        Returns:
            dict: a dictionary containing the information from the complete cbf
                file per scan
        """

        scan_info = {}
        for scan in unique_scans:
            # Get information from first frame
            file_name = os.path.join(frame_dir, all_frames[(scan, 1)]['filename'])
            full_header = full_cbf.extract_full_cbf_header_information(file_name)

            scan_info[scan] = full_header

        return scan_info


    def _get_scan_info_mini_header(self, frame_dir, unique_scans, all_frames):
        """Get the scan information from the mini header.

        Args:
            frame_dir (str): the directory containing the files with frames
            unique_scans (set): a set of unique scans
            all_frames (_type_): a dictionary containing the link between frames
                and filenames.

        Raises:
            Exception: scan range does not match increment

        Returns:
            dict: a dictionary containing the scan information
        """

        scan_info = {}
        axes = imgCIF_creator.ROT_AXES + imgCIF_creator.TRANS_AXES
        frame_type = self._determine_frame_type(
            os.path.join(frame_dir, all_frames[list(all_frames)[0]]['filename']))

        print(f"Discovered {frame_type} files")
        print('Retrieving scan information...')

        for scan in unique_scans:

            scan_frame_map = list(filter(lambda x : x[0] == scan, all_frames.keys()))
            frames = [x[1] for x in scan_frame_map]

            # Get information for first frame
            # print('frdir', frame_dir)
            # print('allfrs', all_frames[(scan, 1)])
            file_name = os.path.join(frame_dir, all_frames[(scan, 1)]['filename'])
            mini_header = self._get_mini_header(file_name, frame_type)
            # print('mini head', mini_header)
            start_axes_settings, scan_ax, scan_incr, exposure, wl,  = \
                self._get_frame_info(mini_header, frame_type, axes)
            x_pixel_size, y_pixel_size, = self._get_pixel_sizes(mini_header)
            # print('pxsz', x_pixel_size, y_pixel_size)
            # print("start_axes_settings ", start_axes_settings, "scanax ", scan_ax, "scaninc ", scan_incr, "exposure ", exposure, "wl ", wl)
            start = start_axes_settings[scan_ax]

            # Get information for last frame
            file_name = \
                os.path.join(frame_dir, all_frames[(scan, len(frames))]['filename'])
            mini_header = self._get_mini_header(file_name, frame_type)
            # print('mini headi', mini_header)
            start_axes_settings, _, _, _, _ = self._get_frame_info(mini_header, frame_type, axes)
            finish = start_axes_settings[scan_ax]

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
            # print('start_axes_settings', start_axes_settings)
            # print('scan', scan)
            # print('scandet', scan_details)
            scan_info[scan] = (start_axes_settings, scan_details)

        self._prune_scan_info(scan_info)

        return scan_info


    def _get_scan_frame_fmt(self, frame_dir, stem=r".*?_"):
        """Deduce the scan/frame naming convention.

        Args:
            frame_dir (str): the directory containing the files with frames
            stem (regexp, optional): The constant portion of the scan/frame
                naming convention. Defaults to r".*?_".

        Returns:
            scan_frame_regex (regexp): the regular expression to identiy scans and frames
            all_names (list): a list of all filenames
        """

        file_pattern = re.compile(stem)
        all_names = []
        for _, _, files in os.walk(frame_dir):
            for filename in files:
                all_names.append(filename)

        # filter out only .cbf and .img files
        all_names = list(filter(
            lambda f_name : f_name.endswith(".cbf") or f_name.endswith(".img"),
            all_names))
        # print('all_names', all_names)

        # if given a file stem filter out only the files that start with the stem/file_pattern
        all_names = list(filter(lambda f_name : file_pattern.match(f_name), all_names))
        all_names.sort()
        # print('all names', all_names)

        # Analyse number of digits between stem and extension: if less than
        # 5, no scan present the default stem matches everything until the _
        test_name = all_names[0]
        stem_len = len(file_pattern.match(test_name).group(0))
        num_digits = len(re.sub("[^0-9]", "", test_name[(stem_len-1):-4]))
        # num_digits = count(r"[0-9]", test_name[stem_len:end-4])

        if num_digits >= 5:
            scan_frame_regex = None
            # Allow first scan not to be scan 1
            for scan in [str(i) for i in range(1,10)]:
                regex = r"(?:" + stem + r")" +\
                    r"(?:(?P<scan>[0-9]*" + re.escape(scan) + \
                    r")(?P<sep>0|_)(?P<frame>[0-9]+1)(?P<ext>\.cbf|img))"
                # rr = stem * r"(?<scan>[0-9]*$scan)(?<sep>0|_)(?<frame>[0-9]+1)\\.(cbf|img)"
                # keep in mind that re.match matche only the beginning of string
                match = re.match(regex, all_names[0])
                # m = match(rr,all_names[1])

                if match:
                    scan_len = len(match.group("scan"))
                    frame_len = len(match.group("frame"))

                    # if the separator is a 0, include it
                    if match.group("sep") == "0":
                        frame_len += 1

                    # scan_frame_regex = stem * Regex("(?<scan>[0-9]{$scan_len})_?(?<frame>[0-9]{$frame_len})")
                    scan_frame_regex = r"(?:" + stem + r")" +\
                        r"(?:(?P<scan>[0-9]{" + re.escape(str(scan_len)) + r"})" +\
                        r"_?(?P<frame>[0-9]{" + re.escape(str(frame_len)) + r"}))"

                    # print('scrfreg', scan_frame_regex)
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

        print(f'\nFound scan/frame naming convention!')

        return scan_frame_regex, all_names


    def _determine_frame_type(self, filename):
        """Determine the type of a frame file: currently SMV (ADSC) and CBF are
        recognised.

        Args:
            filename (str): the filename for which the frame type should be
                determined

        Returns:
            str: the fileformat
        """

        with open(filename, 'rb') as file:
            # read first 512 characters/bytes as byte string
            header = file.read(512)
            # TODO ensure this! maybe its also 0x0c for the form feed character
            if b'\f' in header:
                return 'SMV'
            elif b'_array_data' in header:
                return 'CBF'


    def _get_frame_info(self, mini_header, frame_type, axes):
        """Choose the method to extract frame information according to the fileformat

        Args:
            mini_header (list): the lines of the miniheader
            frame_type (str): the frame type
            axes (tuple): a tuple of axis names

        Returns:
            method: the suitable method to extract the frame info
        """

        if frame_type == "CBF":
            return self._get_frame_info_CBF(mini_header, axes)
        elif frame_type == "SMV":
            return self._get_frame_info_SMV(mini_header)


    def _get_mini_header(self, filename, frame_type):
        """Choose the method to extract mini header information according to the
        fileformat

        Args:
            filename (str): the filename for which the frame type should be
                determined
            frame_type (str): the type of the frames

        Returns:
            method: the suitable method to extract the mini header
        """

        if frame_type == "CBF":
            return self._get_CBF_header(filename)
        elif frame_type == "SMV":
            return self._get_SMV_header(filename)


    def _get_CBF_header(self, filename):
        """Return the lines of the mini header of a cbf file.

        Args:
            filename (str): the cbf filename

        Returns:
            list: a list containing the lines of the mini header
        """

        with open(filename, 'br') as file:
            cbf_header = []
            found_header = False
            within_mini_header = False
            for line in file:
                if b'_array_data.header_contents' in line:
                    found_header = True
                # consider that the mini header is enclosed within two lines of ;
                elif (line == b';\n') or (line == b';\r\n') or (line == b';\r'):
                    within_mini_header = not within_mini_header
                elif b"-BINARY-FORMAT-SECTION-" in line:
                    break

                if found_header and within_mini_header and line.startswith(b'#'):
                    cbf_header.append(line.decode('utf-8').lower().strip('\n').strip('\r'))

        # print('lns', cbf_header)
        return cbf_header


    def _get_frame_info_CBF(self, cbf_header, axes):
        """Return any values found for provided axes. All axes converted to lowercase.
        Return also the exposure time and the wavelength and determine the scan
        axis.

        Args:
            cbf_header (list): a list of lines from the cbf mini header
            axes (tuple): the axes for which the values should be retrieved

        Raises:
            Exception: if the scan axis found does not match with an axis name

        Returns:
            ax_vals (dict): the retrieved values for the given axes
            matching_scan_ax (str): the axis indentified as scan axis that also
                matched with the given axes
            scan_ax[1] (float): the scan axis increment
            exposure (float): the exposure time in seconds
            wl (float): the wavelength in Angstrom
        """

        ax_vals = \
            list(map(lambda ax : (ax.lower(), self._get_CBF_header_values(cbf_header, ax)),
                    axes))
        ax_vals = list(filter( lambda x : x[1][0] != None, ax_vals))
        ax_vals = [self._convert_units(ax_val, 'length') for ax_val in ax_vals]
        # print('ax vals', ax_vals)

        ax_incr = list(map(lambda ax : (ax.lower(), self._get_CBF_header_values(
            cbf_header, ax + "_increment")), axes))
        ax_incr = list(filter( lambda x : x[1][0] != None, ax_incr))
        ax_incr = [self._convert_units(incr, 'length') for incr in ax_incr]
        # print('axinc 2', ax_incr)

        _, exposure = self._convert_units(('et', self._get_CBF_header_values(
            cbf_header, "exposure_time")), 'time')
        _, wl = self._convert_units(('wl', self._get_CBF_header_values(
            cbf_header, "wavelength")), 'wavelength')

        # find the first element whose increment changes
        # scan_ax = findfirst( x -> !isapprox(x[2], 0, atol = 1e-6), ax_incr)
        scan_ax = next(filter(lambda x: not np.isclose(x[1], 0, atol=1e-6), ax_incr), None)
        # print('scax', scan_ax)

        matching_scan_ax = list(filter(lambda ax : scan_ax[0] in ax[0], ax_vals))
        if matching_scan_ax == []:
            raise Exception(
                f"Could not match the scan axis found ({scan_ax}) with an axis name.")
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
        # TODO
        # For a single-axis diffractometer currently

        smv_header = self._get_SMV_header(filename)

        ax_vals = [("phi", self._get_SMV_header_values(smv_header, "phi"))]
        ax_vals.append(("trans", self._get_SMV_header_values(smv_header, "distance")))
        ax_incr = [("phi", self._get_SMV_header_values(smv_header, "osc_range"))]
        ax_incr.append(("trans", 0.0))

        exposure = self._get_SMV_header_values(smv_header, "time")
        wl = self._get_SMV_header_values(smv_header, "wavelength")

        return dict(ax_vals), "phi", ax_incr[0][1], exposure, wl


    def _get_CBF_header_values(self, lines, matcher):
        """Get the value following the string given in matcher and units if present.

        Args:
            lines (list): the list of lines from the mini header
            matcher (str): the string that should be matched in the lines

        Returns:
            val (float): the value that has been matched
            units (str): the units of the value that has been matched
        """

        pattern = re.compile(re.escape(matcher) + r"[ =]+")
        matching_lines = list(filter(lambda x : pattern.search(x) is not None, lines))
        # print('mtch lng', matching_lines)

        # TODO is it now ensured that there are not more than one lines matching?
        if len(matching_lines) < 1:
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
        return float(val), units


    def _get_pixel_sizes(self, lines):
        """Return the pixel sized found in the given lines.

        Args:
            lines (list): a list of lines found in the mini header

        Returns:
            x_pixel (float): the x pixel size in mm
            y_pixel (float): the y pixel size in mm
        """

        matcher = 'pixel_size'
        pattern = re.compile(re.escape(matcher) + r"[ =]+")
        matching_lines = list(filter(lambda x : pattern.search(x) is not None, lines))
        if len(matching_lines) < 1:
            return None, None, None, None

        dim = r"([A-Za-z0-9+-.]+) +([A-Za-z.]+)"
        val_unit_regex = re.escape(matcher) + r'[ =]+' + dim + r' [A-Za-z] ' + dim
        val_unit = [re.search(val_unit_regex, matching_line)
            for matching_line in matching_lines]
        # print('valuns', val_unit)
        val_unit = val_unit[0]
        pixel_sizes = [group.strip() for group in val_unit.groups()]

        # as it is in mm now, we round to 5 decimal places since otherwise:
        # >>> 0.000172*1000 = 0.17200000000000001
        _, x_pixel = self._convert_units(('x_size',
            (float(pixel_sizes[0]), pixel_sizes[1])), 'length')
        _, y_pixel = self._convert_units(('y_size',
            (float(pixel_sizes[2]), pixel_sizes[3])), 'length')
        # print(x_pixel, y_pixel)

        return round(x_pixel, 5), round(y_pixel, 5)


    def _get_SMV_header(self, filename):

        # TODO look at smv files
        with open(filename, 'rb') as file:
            header = file.read(512)
            smv_header = header.split("\n").lower()

        return smv_header


    def _get_SMV_header_values(self, lines, matcher):
        """
        Get the value following the string given in matcher and units if present
        """
        #TODO

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


    def _convert_units(self, ax_val, val_type):
        """Convert into units used in the imgCIF format (mm, s, Angstrom)

        Args:
            ax_val (tuple): the name, (value, unit) of the entry to convert
            val_type (str): the type of the value (to distingusish between lenghts
                and wavelengths)

        Returns:
            name (str): the name of the converted value
            val (float): the value after conversion
        """

        conversion_map = {('length', 'm') : 1e3, ('length', 'cm') : 10,
            ('time', 'ms') : 1e-3, ('time', 'us') : 1e-6,
            ('time', 'ns') : 1e-9, ('wavelength', 'nm') : 1e-1}

        name, (val, unit) = ax_val
        if (val_type, unit) in conversion_map.keys():
            val = val * conversion_map[(val_type, unit)]

        return name, val


    def _rename_axes(self, scan_info, renaming_scheme):
        """
        Axes in `scan_info` and `renaming_scheme` should already be lowercase.
        """
        #TODO is not used right now

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

        #TODO check that this does the right thing
        start_axes_settings, details = scan_info[list(scan_info.keys())[0]]
        # print('iits', start_axes_settings, 'dets', details)
        scan_axis = details["axis"]
        keep_this = [scan_axis]
        for name, ini_val in start_axes_settings.items():
            for content in scan_info.values():
                if content[0][name] != ini_val:
                    keep_this.append(name)
                    break
        # print('wann kepp', keep_this)

        for name, ini_val in start_axes_settings.items():
            if not (name in imgCIF_creator.ALWAYS_AXES) and not (name in keep_this) \
                and np.isclose(ini_val, 0, atol=0.001):

                for s in scan_info:
                    print('I pruned:', scan_info[s][0], name)
                    del(scan_info[s][0], name)
