import re
import importlib
import numpy as np
from imgCIF_Creator.command_line_interfaces import parser

# Configuration information
#TODO always?
ROT_AXES = ("chi", "phi", "detector_2theta", "two_theta", "omega",
            "angle", "start_angle", "kappa")
TRANS_AXES = ("detector_distance", "dx", "trans", "distance")
ALWAYS_AXES = ("distance", "two_theta", "detector_2theta")


class imgCIFCreator:
    """See documentation of the __init__ method.
    """

    def __init__(self, filename, filetype, stem) -> None:
        """Initialize the imgCIFCreator and dynamically load the appropriate
        extractor module.

        Args:
            filename (str): The filename or directory where the data is located.
            filetype (str): The filetype (smv, cbf or h5)
            stem (str): constant portion of the filenames to determine the scan
                frame naming convention.
        """

        if filetype in ['smv', 'cbf']:
            extractor_module = importlib.import_module(
                'imgCIF_Creator.information_extractors.cbf_smv')
        elif filetype == 'h5':
            extractor_module = importlib.import_module(
                'imgCIF_Creator.information_extractors.hdf5_NxMx')

        self.extractor = extractor_module.extractor(filename, stem)
        self.cmd_parser = parser.CommandLineParser()
        self.generators = imgCIFEntryGenerators()


    def create_imgCIF(self, cif_block, external_url, prepend_dir, filename, filetype):
        """Add the required information to a cif_block and request the missing
        information from the user.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): A cif block created with
                the pycifrw package to which the information is added.
            external_url (str): An external url of the files e.g. a zeondo url
            prepend_dir (str): If the directory name is included as part of the archive
                path name this is the prepended directory name.
            filename (str): The filename or directory where the data is located.
            filetype (str): The filetype (smv, cbf or h5)
        """

        # _diffrn_source block
        source_info = self.extractor.get_source_info()
        source_info = self.check_source_completeness(source_info)
        self.generators.generate_source(cif_block, source_info)

        uncat_info = self.extractor.get_uncategorized_info()
        uncat_info = self.check_uncategorized_completeness(uncat_info)
        self.generators.generate_uncategorized(cif_block, uncat_info)

        # self.generate _diffrn_wavelength block
        radiation_info = self.extractor.get_radiation_info()
        radiation_info = self.check_radiation_completeness(radiation_info)
        self.generators.generate_radiation(cif_block, radiation_info)

        # describe _axis block
        axes_info = self.extractor.get_axes_info()
        axes_info = self.check_axes_completeness(axes_info)
        self.generators.generate_axes(cif_block, axes_info)

        # # describe _array_structure_list_axis and _array_structure_list
        array_info = self.extractor.get_array_info()
        array_info = self.check_array_completeness(array_info)
        self.generators.generate_array(cif_block, array_info)

        # describe _diffrn_detector and _diffrn_detector_axis
        # this correlates with the detector axes in generate axes!
        detector_info = self.extractor.get_detector_info()
        detector_info = self.check_detector_completeness(detector_info)
        self.generators.generate_detector(cif_block, detector_info)

        scan_setting_info = self.extractor.get_scan_settings_info()
        # TODO this is not checking anything right now
        scan_setting_info = self.check_scan_settings_completeness(scan_setting_info)
        scan_list = self.generate_scan_list(scan_setting_info)

        # self.generate _diffrn_scan_axis block
        self.generators.generate_scan_settings(cif_block, scan_setting_info)

        # self.generate _diffrn_scan block
        self.generators.generate_scan_info(cif_block, scan_list)

        # self.generate _diffrn_scan_frame block
        self.generators.generate_step_info(cif_block, scan_setting_info, scan_list)

        # self.generate _diffrn_data_frame block
        self.generators.generate_data_frame_info(cif_block, scan_list)

        # self.generate _array_data block
        self.generators.generate_ids(cif_block, scan_list)

        # self.generate _array_data_external_data
        archive, external_url = self.get_archive_type(external_url, filename)
        self.generators.generate_external_ids(
            cif_block, external_url, self.extractor.all_frames,
            scan_list, archive, prepend_dir, filetype)


    def check_uncategorized_completeness(self, uncat_info):
        """Check if the uncategorized information is complete and request input
        if not.

        Args:
            uncat_info (dict): Some uncategorized information that is needed.

        Returns:
            dict: the information completed
        """

        if self.param_is_none(uncat_info['doi']):
            uncat_info['doi'] = self.cmd_parser.request_input('doi')

        return self.lists_to_values(uncat_info)


    def check_source_completeness(self, source_info):
        """Check if the source information is complete and request input
        if not.
        Args:
            source_info (dict): information about the source

        Returns:
            dict: the information completed
        """

        layout = ''
        if not any(source_info.values()):
            layout = self.cmd_parser.request_input('layout')

        if layout.lower() == 'beamline' or \
            (source_info.get('beamline') is not None) or \
            (source_info.get('facility') is not None):
            required_info = ['facility', 'beamline']
            print('\nCreating a imgCIF file for a beamline.')

        if layout.lower() == 'laboratory' or \
            (source_info.get('manufacturer') is not None) or \
            (source_info.get('model') is not None) or \
            (source_info.get('location') is not None):
            required_info = ['manufacturer', 'model', 'location']
            print('\nCreating a imgCIF file for a laboratory setup.')


        for info in required_info:
            if source_info.get(info) is None:
                source_info[info] = self.cmd_parser.request_input(info)
            else:
                print(f'\nFound the following information about the {info}: \
{source_info[info]}')
        # print('')

        # base = "_diffrn_source."
        # if "Beamline name" in raw_info or "Facility name" in raw_info:
        #     cif_block[base + "beamline"] = [raw_info["Beamline name"]]
        #     cif_block[base + "facility"] = [raw_info["Facility name"]]
        # else:
        #     cif_block[base + "make"] = [raw_info["Name of manufacturer"] + "-" + \
        #         raw_info["Model"]]
        #     if "Location" in raw_info:
        #         cif_block[base + "details"] = [f"Located at {raw_info['Location']}"]

        return self.lists_to_values(source_info)


    def check_axes_completeness(self, axes_info):
        """Check if the axes information is complete and request input if not.

        Args:
            axes_info (dict):

        Returns:
            dict: the information completed
        """

        # print('axinf', axes_info)
        lengths = [len(values) for values in axes_info.values() if values is not None]
        hast_same_len = all([length == lengths[0] for length in lengths])

        missing_information = False
        for prop, value in axes_info.items():

            if self.param_is_none(value):
                missing_information = True

        if missing_information or not hast_same_len:
            print('\nSome information about the goiometer/detector is missing, please enter\
the missing information.')

            print('\nGoniometer information: \n\
Answer the goniometer questions for all axes in zero position.')
            goniometer_axes = self.cmd_parser.request_input('goniometer_axes')
            goniometer_axes = self.cmd_parser.parse_axis_string(goniometer_axes)

            new_regex_stem = r'(' + r'|'.join(goniometer_axes[0]) + r')'
            # make case insensitve
            new_regex = r'(?i)(' + new_regex_stem + r'((\s|,)(\s)*\d{1,3}){1,2}' + r')\Z'
            self.cmd_parser.validation_regex['kappa_axis'] = new_regex
            kappa_axis = self.cmd_parser.request_input('kappa_axis')

            new_regex = r'(?i)(' + new_regex_stem + r'((\s|,)(\s)*\d{1,3})' + r')\Z'
            self.cmd_parser.validation_regex['chi_axis'] = new_regex
            chi_axis = self.cmd_parser.request_input('chi_axis')

            gon_axes = self.cmd_parser.make_goniometer_axes(
                goniometer_axes, kappa_axis, chi_axis)

            print('\nDetector information: \nAnswer the following questions assuming \
all detector positioning axes are at their home positions. If there is more than \
one detector, or detectors are non-rectangular, please describe in the comments \
(not implemented yet).')
            principal_orientation = self.cmd_parser.request_input('principal_orientation')
            image_orientation = self.cmd_parser.request_input('image_orientation')
            two_theta_axis = self.cmd_parser.request_input('two_theta_axis')

            det_axes = self.cmd_parser.make_detector_axes(
                goniometer_axes, principal_orientation, image_orientation, two_theta_axis)

            for key in  gon_axes.keys():
                gon_axes[key] += det_axes[key]

            # print('gohomoni', gon_axes)
            axes_info = gon_axes

        return axes_info


    def check_array_completeness(self, array_info):
        """Check if the array information is complete and request input
        if not.

        Args:
            array_info (dict): information about the data array

        Returns:
            dict: the completed information
        """

        # TODO check hardcoded stuff
        array_structure_labels = ['axis_id', 'axis_set_id', 'pixel_size']
        if any([self.param_is_none(array_info[label]) for label in array_structure_labels]):

            array_info["axis_id"] = ["array_x", "array_y"] #TODO should that be the
            # detector axes detx, dety
            array_info["axis_set_id"] = [1, 2]

            if self.param_is_none(array_info['pixel_size']):
                array_info['pixel_size'] = self.cmd_parser.request_input('pixel_size')

        array_structure_list_labels = ['array_id', 'array_index', 'array_dimension',
            'array_direction', 'array_precedence']
        if any([self.param_is_none(array_info[label]) for label in array_structure_list_labels]):

            array_info["array_id"] = [1, 1]
            array_info["array_index"] = [1, 2]
            array_info["axis_set_id"] = [1, 2]

            if self.param_is_none(array_info['array_dimension']):
                array_info['array_dimension'] = \
                    self.cmd_parser.request_input('pixel_number').replace(' ', '').split(',')

            array_info["array_direction"] = ["increasing", "increasing"]

            if self.param_is_none(array_info['array_precedence']):
                fast_direction = self.cmd_parser.request_input('fast_direction')
                if fast_direction == "horizontal":
                    array_info['array_precedence'] = [1, 2]
                else:
                    array_info['array_precedence'] = [2, 1]

        return array_info


    def check_detector_completeness(self, detector_info):
        """Check if the detector information is complete and request input
        if not.

        Args:
            detector_info (dict): information about the detector

        Returns:
            dict: the information completed
        """

        # TODO does not include multiple detectors yet
        # print('detrinf', detector_info)

        if self.param_is_none(detector_info["detector_id"]):
            detector_info["detector_id"] = ['det1']

        detector_axes_labels = ['number_of_axes', 'axis_id']

        if any([self.param_is_none(detector_info[label]) for label in detector_axes_labels]):
            det_axes = self.cmd_parser.request_input('detector_axes').split(',')
            det_axes = [axis.strip() for axis in det_axes]
            det_axes = [axis for axis in det_axes if axis!='']

            # TODO testing len is always correct
            detector_info["number_of_axes"] = [len(det_axes)]
            # TODO is that true? in the hdf5 this is trans and in the cbf this is
            # detx and dety
            detector_info["axis_id"] = det_axes

        #TODO it's not ensured that this is always a list...
        detector_info["detector_axis_id"] = \
            detector_info["detector_id"] * len(detector_info["axis_id"])

        return detector_info


    def check_radiation_completeness(self, radiation_info):
        """Check if the wavelength information is complete and request input
        if not.

        Args:
            radiation_info (dict): the information about the radiation

        Returns:
            dict: the informatio completed
        """

        if self.param_is_none(radiation_info['rad_type']):
            radiation_info['rad_type'] = \
                self.cmd_parser.request_input('rad_type')
            if radiation_info['rad_type'] == '':
                radiation_info['rad_type'] = 'x-ray'

        if self.param_is_none(radiation_info['wavelength']):
            radiation_info['wavelength'] = \
                self.cmd_parser.request_input('wavelength')

        return self.lists_to_values(radiation_info)


    def check_scan_settings_completeness(self, scan_settings_info):
        """Check if the scan information is complete and request input
        if not.

        Args:
            scan_settings_info (dict): the scan setting information where the
                key is the scan and the informaton is stored in a tuple containing
                dictionaries with information about the axes and the scan details

        Returns:
            dict: the information completed
        """

        # does not need to ask for wl since this is done in radiation_info
        # same for pixel size, which is done in array_info
        # the other parameters are scan dependent? doesn't make sense to request
        # that from the user?
        # TODO

        return scan_settings_info


    def generate_scan_list(self, scan_setting_info):
        """Generate a list contaning tuples with the scan name and the number of
        frames, e.g. [('01', 325)]

        Args:
            scan_settings_info (dict): the scan setting information where the
                key is the scan and the informaton is stored in a tuple containing
                dictionaries with information about the axes and the scan details

        Returns:
            list: a list consiting of tuples with the scan name and the number of
        frames
        """
        # something like scan_list [('01', 325)]
        # Create scan list of (scanid, frame_no) where
        # frame_no is just the number of frames

        scans = sorted(scan_setting_info)
        slist = [(s, scan_setting_info[s][1]["frames"]) for s in scans]

        return slist


    def param_is_none(self, param):
        """Check if the parameter is None or if the string content is 'none'.
        For lists the result is True if already one entry is None or 'none'.

        Args:
            param (list, int, str): the parameter which shall be checked

        Returns:
            bool: True if the parameters is identified as None else False
        """

        if param is None:
            return True
        elif type(param) == str:
            return param.lower() == 'none'
        elif type(param) == list:
            is_none = []
            for entry in param:
                if entry is None:
                    is_none.append(True)
                elif type(entry) == str:
                    is_none.append(entry.lower() == 'none')
                else:
                    is_none.append(False)

            return any(is_none)

        else:
            return False


    def lists_to_values(self, param_dict):
        """Sometimes values parsed from full cbf files are single lenght lists
        which are converted into the value only with this method.

        Args:
            param_dict (dict): the dictionary containing the information

        Returns:
            dict: the same dict with single lenght lists replaced
        """

        for key, value in param_dict.items():
            if (type(value) == list) and (len(value) == 1):
                param_dict[key] = value[0]
                # print('set', value, 'to', param_dict[key])

        return param_dict


    def get_archive_type(self, external_url, filename):
        """Get the archive type of the provided external url or return the default
        filename/path as external url if no external url was provided.

        Args:
            external_url (str): the external url wich possibly points to an archive
                format
            filename (str): The filename or directory where the data is located.

        Returns:
            tuple: the archive format str and the external url
        """

        archive = None

        if external_url == '':
            filename = filename[1:] if filename.startswith('/') else filename
            external_url = "file://" + filename

        else:
            archives = {"TGZ" : r".*((\.tgz\Z)|(\.tar\.gz\Z|))",
                        "TBZ" : r".*((\.tbz\Z)|(\.tar\.bz2\Z|))",
                        "ZIP" : r".*(\.zip\Z)"}

            for archive, regex in archives.items():
                if re.match(regex, external_url):
                    return archive, external_url

        return archive, external_url


class imgCIFEntryGenerators():

    """See the documentation of the __init__ method
    """

    def __init__(self) -> None:
        """This class provides methods to generate the entries in the imgCIF with
        the extracted information.
        """
        pass

    def generate_uncategorized(self, cif_block, uncat_info):
        """Generate the cif_block entries for the uncategorized information.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            uncat_info (dict): the uncategorized information
        """

        cif_block['_database.dataset_doi'] = uncat_info['doi']
        if not uncat_info.get('overload') is None:
            cif_block['_array_intensities.overload'] = uncat_info['overload']


    def generate_radiation(self, cif_block, radiation_info):
        """Generate the cif_block entries for the radiation information.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            radiation_info (dict): the information about the radiation
        """

        # print('scaninf', scan_info)
        cif_block["_diffrn_radiation.type"] = radiation_info['rad_type']

        base = "_diffrn_radiation_wavelength"
        cif_block[base + ".id"] = [1]
        cif_block[base + ".value"] = \
            [radiation_info['wavelength']]
        cif_block.CreateLoop([base + '.id', base + '.value'])


    def generate_scan_settings(self, cif_block, scan_info):
        """Generate the cif_block entries for the scan settings information.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            scan_info (dict): the information about the scan
        """

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
            #TODO do this for all axes or only the one that changes? >> all
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
                else:
                    settings = empty_settings + settings

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


    def generate_scan_info(self, cif_block, scan_list):
        """Generate the cif_block entries for the scan information.
        Fill in the scan information. We number the frames from
        the start

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            scan_list (list): a list consiting of tuples with the scan name and the number of
                frames
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
            entries[base + ".id"].append(f"SCAN{scan}")
            entries[base + ".frame_id_start"].append(f"frm{start_frame}")
            entries[base + ".frame_id_end"].append(f"frm{end_frame}")
            entries[base + ".frames"].append(frame)

            start_frame = end_frame + 1

        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))


    def generate_ids(self, cif_block, scan_list):
        """Generate the cif_block entries for the ids.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            scan_list (list): a list consiting of tuples with the scan name and the number of
                frames
        """
        # Array_id is just IMAGE (a single detector module). The
        # binary_id is incremented with frame, and all of them
        # are located externally so external_id is also incremented
        # together with binary_id.

        base = '_array_data'
        entries = {
            base + ".id" : [],
            base + ".binary_id" : [],
            base + ".external_data_id" : [],
            }

        counter = 0
        for _, frames in scan_list:
            for _ in range(1, frames + 1):
                counter += 1
                entries[base + ".id"].append("IMAGE01") #TODO is that always IMAGE?
                entries[base + ".binary_id"].append(counter)
                entries[base + ".external_data_id"].append(f'ext{counter}') # TODO put ext here?

        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))


    # new url from get arch typem all frames from get all frames, scan linst same
    # prepend dir from creator call, file format from creator
    def generate_external_ids(self, cif_block, external_url, all_frames, scan_list,
                              arch, prepend_dir, file_format):
        """Generate the cif_block entries for the external ids. `external_url` is
        the location of a single archive file. The individual files
        on the local storage are assumed to be at the same locations relative to
        this top-level directory.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            external_url (str): the location of a single archive file
            all_frames (dict): a dictionary with tuples of (scan, frame) as keys
                and the filenames as values
            scan_list (list): a list consisting of tuples with the scan name and
                the number of frames
            arch (str): the archive format of the external url
            prepend_dir (str): the directory that is prepended to the path
            file_format (str): the file format (cbf, smv, h5)
        """
        # print('my arch', arch)
        # print('prepdoi', prepend_dir)
        # print('allfrms', all_frames)

        base = '_array_data_external_data'
        entries = {
            base + ".id" : [],
            base + ".format" : [],
            base + ".uri" : [],
            }

        if file_format == 'h5':
            entries[base + ".path"] = []
            entries[base + ".frame"] = []

        if arch is not None:
            entries[base + ".archive_format"] = []
            entries[base + ".archive_path"] = []

        counter = 0
        for scan, frames in scan_list:
            # print('sc fr scl', scan, frames, scan_list)
            for frame in range(1, frames + 1):
                # print('allfr, scan, fram', all_frames[(scan, frame)])
                counter += 1
                frame_name = f"/{all_frames[(scan, frame)]['filename']}"

                entries[base + ".id"].append(f"ext{counter}")
                file_format = 'HDF5' if file_format == 'h5' else file_format
                entries[base + ".format"].append(file_format.upper())
                entries[base + ".uri"].append(external_url) #TODO zenodo uri here?

                # TODO path and frame only for hdf5? repsectively containerized images?
                if file_format == 'h5' or file_format == 'HDF5':
                    entries[base + ".path"].append(all_frames[(scan, frame)]['path'])
                    entries[base + ".frame"].append(all_frames[(scan, frame)]['frame'])

                # A too-clever-by-half way of optionally live constructing a URL
                if arch is not None:
                    entries[base + ".archive_format"].append(arch)
                    entries[base + ".archive_path"].append(prepend_dir + frame_name)

        for name, value in entries.items():
            # print('name', name, 'value', value)
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))


    def generate_step_info(self, cif_block, scan_setting_info, scan_list):
        """Generate the cif_block entries for the step information.
        Fill in information about steps

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            scan_settings_info (dict): the scan setting information where the
                key is the scan and the informaton is stored in a tuple containing
                dictionaries with information about the axes and the scan details
            scan_list (list): a list consiting of tuples with the scan name and the number of
                frames
        """

        base = '_diffrn_scan_frame'
        entries = {
            base + ".frame_id" : [],
            base + ".scan_id" : [],
            base + ".frame_number" : [],
            base + ".integration_time" : [],
            }

        # println(op, header)
        counter = 0
        for scan, frames in scan_list:
            for frame_number in range(1, frames + 1):
                counter += 1
                # println(op, "frm$ctr   SCAN$s    $f $(scan_times[s])")
                entries[base + ".frame_id"].append(f"frm{counter}")
                entries[base + ".scan_id"].append(f"SCAN{scan}")
                entries[base + ".frame_number"].append(frame_number)
                entries[base + ".integration_time"]\
                    .append(scan_setting_info[scan][1]['time'])

        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))



    def generate_source(self, cif_block, source_info):
        """Generate the cif_block entries for the source information.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            source_info (dict): information about the source
        """

        regex = r"[^A-Za-z0-9_]"
        base = "_diffrn_source."
        if source_info.get("beamline") is not None or\
             source_info.get("facility") is not None:
            # print('facilit', source_info["facility"])
            cif_block[base + "beamline"] = source_info["beamline"]
            cif_block[base + "facility"] = source_info["facility"]
            audit_id_1 = re.sub(regex, '_', source_info["beamline"])
            audit_id_2 = re.sub(regex, '_', source_info["facility"])
        else:
            cif_block[base + "make"] = source_info["manufacturer"] + "-" + \
                source_info["model"]
            cif_block[base + "details"] = 'Located at' + f"{source_info['location']}"
                # f"Located at {source_info['location']}"
            audit_id_1 = re.sub(regex, '_', source_info["manufacturer"])
            audit_id_2 = re.sub(regex, '_', source_info["model"])

        cif_block['_audit.block_code'] = audit_id_1 + '_' + audit_id_2.strip('_')


    def generate_array(self, cif_block, array_info):
        """Generate the cif_block entries for the array information.

        Produce the information required for array_structure_list. Here
        we assume a rectangular detector with x horizontal, y vertical

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            array_info (dict): information about the source
        """
        hor, vert = [float(element) for element in array_info['pixel_size']]
        # array structure list axis
        base = "_array_structure_list_axis."
        # TODO always detx dety?
        cif_block[base + "axis_id"] = array_info['axis_id']
        cif_block[base + "axis_set_id"] = array_info['axis_set_id']
        cif_block[base + "displacement"] = [hor / 2, vert / 2]   #half pixel size
        cif_block[base + "displacement_increment"] = [hor, vert]
        cif_block.CreateLoop([
            base + "axis_id",
            base + "axis_set_id",
            base + "displacement",
            base + "displacement_increment"])

        # array structure list
        base = "_array_structure_list."
        cif_block[base + "array_id"] = array_info['array_id']
        cif_block[base + "index"] = array_info['array_index']
        cif_block[base + "axis_set_id"] = array_info['axis_set_id']
        cif_block[base + "dimension"] = array_info['array_dimension']# ['none', 'none']   #number of elements in each direction
        cif_block[base + "direction"] = array_info['array_direction']
        cif_block[base + "precedece"] = array_info['array_precedence']

        # print('arryinfo', array_info)

        cif_block.CreateLoop([
            base + "array_id",
            base + "index",
            base + "axis_set_id",
            base + "dimension",
            base + "direction",
            base + "precedece"])


    def generate_data_frame_info(self, cif_block, scan_list):
        """Generate the cif_block entries for the data frame information.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            scan_list (list): a list consiting of tuples with the scan name and
                the number of
        """

        base = '_diffrn_data_frame'
        entries = {
            base + ".id" : [],
            base + ".detector_element_id" : [],
            base + ".array_id" : [],
            base + ".binary_id" : [],
            }

        # TODO how to include the counter into element and image? where is the information
        # this probably needs the get_array_info method then
        # println(op, header) for now hardcoded
        counter = 0
        for _, frames in scan_list:
            for _ in range(1, frames + 1):
                counter += 1
                entries[base + ".id"].append(f"frm{counter}")
                #TODO hardcoded element, image
                entries[base + ".detector_element_id"].append(f"ELEMENT01")
                #TODO element or 1? not present in hdf5 mapper
                entries[base + ".array_id"].append("IMAGE01") # is 1 in hdf5 mapper
                entries[base + ".binary_id"].append(counter)

        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))


    def generate_detector(self, cif_block, detector_info):
        """Generate the cif_block entries for the detector information.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            detector_info (dict): information about the detector
        """

        base = "_diffrn_detector."
        cif_block[base + "id"] = detector_info["detector_id"]
        cif_block[base + "number_of_axes"] = detector_info['number_of_axes']
        cif_block.CreateLoop([base + "id", base + "number_of_axes"])
        #
        base = "_diffrn_detector_axis."
        # print('olla', cif_block["_diffrn_detector.number_of_axes"])
        cif_block[base + "axis_id"] = detector_info['axis_id']
        cif_block[base + "detector_id"] = detector_info['detector_axis_id']
        cif_block.CreateLoop([base + "axis_id", base + "detector_id"])


    def generate_axes(self, cif_block, axes_info):
        """Generate the cif_block entries for the axes information.

        Create the goniometer axes corresponding to the data in `raw_info`,
        placing the result in CIF block `cif_block`.

        Args:
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
            axes_info (dict): information about the axes
        """
        # print("gonaxes2", gon_axes)
        base = "_axis."
        cif_block[base + "id"] = axes_info['axes']
        cif_block[base + "type"] = axes_info['axis_type']
        cif_block[base + "equipment"] = axes_info['equip']
        cif_block[base + "depends_on"] = axes_info['depends_on']
        self.split_vectors(axes_info['vector'], base + "vector", cif_block)
        self.split_vectors(axes_info['offset'], base + "offset", cif_block)
        cif_block.CreateLoop([base + "id", base + "type", base + "equipment",
                            base + "depends_on",
                            base + "vector[1]", base + "vector[2]", base + "vector[3]",
                            base + "offset[1]", base + "offset[2]", base + "offset[3]"])


    def split_vectors(self, vectors, basename, cif_block):
        """
        Distribute the vectors in `vectors` over CIF data names, formatting as
        integers if they are exact.

        Args:
            vectors (list): a list of vectors which are lists of 3 entries itself
            basename (str): the basename of the imgCIF entry
            cif_block (CifFile.CifFile_module.CifBlock): the cif block into which
                the information should be written.
        """
        # print('veccoi', vectors)
        # for idx_1, sub_vectors in enumerate(vectors):
        #     for idx_2, element in enumerate(sub_vectors):
        #         vectors[idx_1][idx_2] = float(element)

        # TODO is the index correct? julia indexes different than python
        def mapping_func(x):
            # i is local from the scope of the loop
            try:
                # the vector indices in cif start at 1
                x[i - 1] = float(x[i - 1])
                # display numbers without trailing zero? 1.0 -> 1
                if round(x[i - 1]) == x[i - 1]:
                    return f"{round(x[i - 1]):.0f}"
                else:
                    # don't format if it has less than 5 decimal places
                    if len(str(x[i - 1]).split('.')[1]) < 5:
                        formatted = str(x[i - 1])
                    else:
                        formatted = f"{x[i - 1]:1.5f}"
                    # it can be that some values are by computaional effects
                    # are very close to zero e.g. -1.2246467991473532e-16 -> map to 0
                    if float(formatted) == 0:
                        formatted = "0"

                    return formatted # x[i] @sprintf "%1.5f" x[i]
            except ValueError:
                # print('exepted', x[i-1])
                return x[i - 1]

        for i in range(1, 4):
            # print('vcodor', vectors)
            # mppd = map(mapping_func, vectors)
            # print('mppd', i, list(mppd))
            cif_block[basename + f"[{i}]"] = map(mapping_func, vectors)

