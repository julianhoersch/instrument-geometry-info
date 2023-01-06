import re
import importlib
import os
import sys
# from imgCIF_Creator.information_extractors import cbf_smv
from imgCIF_Creator.command_line_interfaces import parser
from imgCIF_Creator.information_extractors import user_input

# Configuration information
#TODO always?
ROT_AXES = ("chi", "phi", "detector_2theta", "two_theta", "omega",
            "angle", "start_angle", "kappa")
TRANS_AXES = ("detector_distance", "dx", "trans", "distance")
ALWAYS_AXES = ("distance", "two_theta", "detector_2theta")

# #============= Output routines ========================#

class imgCIFAssembler:

    def __init__(self, filename, filetype, stem) -> None:

        if filetype in ['smv', 'cbf']:
            extractor_module = importlib.import_module(
                'imgCIF_Creator.information_extractors.cbf_smv')
        elif filetype == 'h5':
            extractor_module = importlib.import_module(
                'imgCIF_Creator.information_extractors.hdf5_NxMx')

        # try:
        #     extractor_module= importlib.import_module(extractor_filename)
        # except ImportError:
        #     pass
        #     # self.extractor = importlib.import_module(
        #     #     "onda.processing_layer.{0}".format(processing_layer_filename)
        #     # )

        self.extractor = extractor_module.extractor(filename, stem)
        # self.extractor = cbf_smv.extractor(filename, stem)
        self.cmd_parser = parser.CommandLineParser()


    # TODO add doi somewhere
    def create_imgCIF(self, cif_block, external_url, prepend_dir, filename, filetype):

        # cif_block['_database.dataset_doi'] =

        # # _diffrn_source block
        # # needs beamline name, facility or manufacturer, model opt location
        source_info = self.extractor.get_source_info()
        source_info = self.check_source_completeness(source_info)
        self.generate_facility(cif_block, source_info) #source
        # generate_facility(user_input, cif_block)

        misc_info = self.extractor.get_miscellaneous_info()
        misc_info = self.check_misc_completeness(misc_info)
        self.generate_misc(cif_block, misc_info)

        # self.generate _diffrn_wavelength block
        # needs rad_type, wavelength
        wavelength_info = self.extractor.get_wavelength_info()
        wavelength_info = self.check_wavelength_completeness(wavelength_info)
        # self.generate_wavelength(cif_block, scan_info)
        self.generate_wavelength(cif_block, wavelength_info)

        # # describe _axis block
        # # needs goniometer axes
        axes_info = self.extractor.get_axes_info()
        axes_info = self.check_axes_completeness(axes_info)
        self.generate_axes(cif_block, axes_info)
        # self.generate()...
        # self.generate_axes(user_input, cif_block)


        # # describe _array_structure_list_axis and _array_structure_list
        # # needs pixel size, array dimension, fast direction
        array_info = self.extractor.get_array_info()
        array_info = self.check_array_completeness(array_info)
        self.generate_array(cif_block, array_info)

        # # return

        # # describe _diffrn_detector and _diffrn_detector_axis
        # # needs detector axes,
        # this correlates with the detector axes in generate axes!
        detector_info = self.extractor.get_detector_info()
        detector_info = self.check_detector_completeness(detector_info)
        self.generate_detector(cif_block, detector_info)
        # self.generate_detector(user_input, cif_block)


        scan_setting_info = self.extractor.get_scan_settings_info()
        # this is not checking anything right now
        scan_setting_info = self.check_scan_settings_completeness(scan_setting_info)
        scan_list = self.get_scan_list(scan_setting_info)


        # self.generate _diffrn_scan_axis block
        # needs scan info
        self.generate_scan_settings(cif_block, scan_setting_info)

        # self.generate _diffrn_scan block
        # needs scan list
        self.generate_scan_info(cif_block, scan_list)

        # self.generate _diffrn_scan_frame block
        # needs scan_list, integration time
        self.generate_step_info(cif_block, scan_setting_info, scan_list)

        # self.generate _diffrn_data_frame block
        # needs scan list
        self.generate_array_info(cif_block, scan_list)

        # self.generate _array_data block
        # needs scan list
        self.generate_ids(cif_block, scan_list)


        archive, external_url = self.get_archive_type(external_url, filename)

        # self.generate _array_data_external_data
        # needs uri/zenod?, all frames, scan list, arch, prepend dir, file format
        self.generate_external_ids(
            cif_block, external_url, self.extractor.all_frames,
            scan_list, archive, prepend_dir, filetype)




    def get_scan_list(self, scan_info):

        # move to assebler

        # something like scan_list [('01', 325)]

        # Create scan list of (scanid, frame_no) where
        # frame_no is just the number of frames

        scans = sorted(scan_info)
        slist = [(s, scan_info[s][1]["frames"]) for s in scans]

        return slist

    def check_misc_completeness(self, misc_info):

        if self.param_is_none(misc_info['doi']):
            misc_info['doi'] = self.cmd_parser.request_input('doi')

        return self.lists_to_values(misc_info)


    def check_source_completeness(self, source_info):

        layout = ''
        if not any(source_info.values()):
            layout = self.cmd_parser.request_input('layout')

        if layout == 'Beamline' or \
            (source_info.get('beamline') is not None) or \
            (source_info.get('facility') is not None):
            required_info = ['facility', 'beamline']
            print('Creating a imgCIF file for a beamline.')

        if layout == 'Laboratory' or \
            (source_info.get('manufacturer') is not None) or \
            (source_info.get('model') is not None) or \
            (source_info.get('location') is not None):
            required_info = ['manufacturer', 'model', 'location']
            print('Creating a imgCIF file for a laboratory setup.')


        for info in required_info:
            # present_info = \
            #     [info for info in required_info if source_info.get(info) is not None]
#             for info in present_info:
#                 print(f'Found the following information about the {info}: \
# {source_info[info]}')

            if source_info.get(info) is None:
                source_info[info] = self.cmd_parser.request_input(info)
            else:
                print(f'\nFound the following information about the {info}: \
{source_info[info]}', end='')

        # base = "_diffrn_source."
        # if "Beamline name" in raw_info or "Facility name" in raw_info:
        #     imgblock[base + "beamline"] = [raw_info["Beamline name"]]
        #     imgblock[base + "facility"] = [raw_info["Facility name"]]
        # else:
        #     imgblock[base + "make"] = [raw_info["Name of manufacturer"] + "-" + \
        #         raw_info["Model"]]
        #     if "Location" in raw_info:
        #         imgblock[base + "details"] = [f"Located at {raw_info['Location']}"]

        return self.lists_to_values(source_info)


    def check_axes_completeness(self, axes_info):

        print('axinf', axes_info)
        lengths = [len(values) for values in axes_info.values() if values is not None]
        hast_same_len = all([length == lengths[0] for length in lengths])

        missing_information = False
        for prop, value in axes_info.items():

            if self.param_is_none(value):
                missing_information = True

            # if value is None:
            #     missing_information = True
            # elif any(elem is None for elem in value):
            #     missing_information = True


        if missing_information or not hast_same_len:
            # print(f'Missing information about the goniometer axes: {prop}.')
            print('Some information about the goiometer/detector is missing, please enter\
the missing information:')

            # if not make gonio axes and make detector axes
            # gonio axes needs
            # goniometer_axes
            # eventually kappa and chi
            print('Goniometer information: \
Answer the following questions for all axes in zero position.')
            goniometer_axes = self.cmd_parser.request_input('goniometer_axes')
            goniometer_axes = user_input.parse_axis_string(goniometer_axes)
            kappa_axis = self.cmd_parser.request_input('kappa_axis')
            chi_axis = self.cmd_parser.request_input('chi_axis')

            gon_axes = user_input.make_goniometer_axes(
                goniometer_axes, kappa_axis, chi_axis)


            # detector axes needs
            # goniometer_axes
            # principal_orientation
            # image_orientation
            # two_theta_axis
            print('Detector information: Answer the following questions assuming \
all detector positioning axes are at their home positions. If there is more than \
one detector, or detectors are non-rectangular, please describe in the comments below.')
            principal_orientation = self.cmd_parser.request_input('principal_orientation')
            image_orientation = self.cmd_parser.request_input('image_orientation')
            two_theta_axis = self.cmd_parser.request_input('two_theta_axis')

            det_axes = user_input.make_detector_axes(
                goniometer_axes, principal_orientation, image_orientation, two_theta_axis)


            # for i, _ in enumerate(gon_axes):
            #     gon_axes[i] += det_axes[i]

            for key in  gon_axes.keys():
                gon_axes[key] += det_axes[key]

            # print('gohomoni', gon_axes)
            axes_info = gon_axes
            # gohomoni [['chi', 'omega', 'beta', 'two_theta', 'trans', 'detx', 'dety'],
            #  ['rotation', 'rotation', 'rotation', 'translation', 'rotation',
            # 'translation', 'translation'], ['goniometer', 'goniometer', 'goniometer',
            # 'detector', 'detector', 'detector', 'detector'], ['omega', 'beta', 'none',
            # 'none', 'two_theta', 'trans', 'detx'], [[-1, 0, 0], [1, 0, 0], [1, 0, 0],
            # [-1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]], [[0, 0, 0], [0, 0, 0],
            # [0, 0, 0], [0, 0, 0], [0, 0, 0], ['none', 'none', 0], [0, 0, 0]]]


        # axes_info = {'axes' : gon_axes[0],
        #              'axis_type' : gon_axes[1],
        #              'equip' : gon_axes[2],
        #              'depends_on' : gon_axes[3],
        #              'vector' : gon_axes[4],
        #              'offset' : gon_axes[5]}

        return axes_info


    def check_array_completeness(self, array_info):


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

            array_info["array_id"] = ['array_x', 'array_y']
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

        # TODO does not include multiople detectors yet
        print('detrinf', detector_info)


        if self.param_is_none(detector_info["detector_id"]):
            detector_info["detector_id"] = ['det1']

        detector_axes_labels = ['number_of_axes', 'axis_id', 'detector_id',
            'detector_axis_id']

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
                detector_info["detector_id"] * len(det_axes)

        return detector_info


    def check_wavelength_completeness(self, wavelength_info):

        # {'rad_type' : rad_type,
        # 'wavelength' : wavelength}

        if self.param_is_none(wavelength_info['rad_type']):
            wavelength_info['rad_type'] = \
                self.cmd_parser.request_input('rad_type')
        if self.param_is_none(wavelength_info['wavelength']):
            wavelength_info['wavelength'] = \
                self.cmd_parser.request_input('wavelength')

        return self.lists_to_values(wavelength_info)


    def check_scan_settings_completeness(self, scan_settings_info):

        # does not need to ask for wl since this is done in wavelength_info
        # same for pixel size, which is done in array_info
        # the other parameters are scan dependent? doesn't make sense to request
        # that from the user?
        # TODO

        # for axes_single_frame, scan_details in scan_settings_info.values():
        #     missing = any(val is None for val in axes_single_frame.values())

        #     missing = any(val is None for val in scan_details.values())


        return scan_settings_info



    def check_completeness(info_block, info_block_name):

        return info_block


    def param_is_none(self, param):

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

        for key, value in param_dict.items():
            if (type(value) == list) and (len(value) == 1):
                param_dict[key] = value[0]
                print('set', value, 'to', param_dict[key])

        return param_dict


    def add_scan_info_to_block(scan_info, all_frames, cif_block, external_url,
                            prepend_dir, directory, file_format):
        """
        We identify each frame by its sequence number overall (not just within
        its own scan.)
        """

        # Output CIF fragment
        arch, external_url = get_archive_type(external_url, directory)
        print('new', external_url, arch)
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
            cif_block, external_url, all_frames, scan_list, arch, prepend_dir, file_format)









    def create_scan_list(scan_info):

        # Create scan list of (scanid, frame_no) where
        # frame_no is just the number of frames

        scans = sorted(scan_info)
        slist = [(s, scan_info[s][1]["frames"]) for s in scans]

        return slist

    def generate_misc(self, cifblock, misc_info):

        cifblock['_database.dataset_doi'] = misc_info['doi']
        if not self.param_is_none(misc_info['overload']):
            cifblock['_array_intensities.overload'] = misc_info['overload']



    def generate_wavelength(self, cifblock, wavelength_info):

        # print('scaninf', scan_info)
        cifblock["_diffrn_radiation.type"] = wavelength_info['rad_type']

        base = "_diffrn_radiation_wavelength"
        cifblock[base + ".id"] = [1]
        cifblock[base + ".value"] = \
            [wavelength_info['wavelength']]
        cifblock.CreateLoop([base + '.id', base + '.value'])



    def generate_scan_settings(self, cif_block, scan_info):
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



    def generate_ids(self, cif_block, scan_list):
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
                entries[base + ".id"].append("IMAGE01") #TODO is that always IMAGE?
                entries[base + ".binary_id"].append(counter)
                entries[base + ".external_data_id"].append(f'ext{counter}') # TODO put ext here?
                # println(op, "   IMAGE $ctr $ctr")

        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))



    # new url from get arch typem all frames from get all frames, scan linst same
    # prepend dir from creator call, file format from creator
    def generate_external_ids(self, cif_block, external_url, all_frames, scan_list,
                              arch, prepend_dir, file_format):
        """
        `external_url` is the location of a single archive file. The individual files
        on the local storage are assumed to be at the same locations relative to
        this top-level directory.
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
            print('sc fr scl', scan, frames, scan_list)
            for frame in range(1, frames + 1):
                # print('allfr, scan, fram', all_frames[(scan, frame)])
                counter += 1
                frame_name = f"/{all_frames[(scan, frame)]['filename']}"

                entries[base + ".id"].append(f"ext{counter}")
                file_format = 'HDF5' if file_format == 'h5' else file_format
                entries[base + ".format"].append(file_format.upper())
                entries[base + ".uri"].append(external_url) #TODO zenodo uri here?

                # TODO path and frame only for hdf5? repsectively containerized images?
                if file_format == 'h5':
                    entries[base + ".path"].append(all_frames[(scan, frame)]['path'])
                    entries[base + ".frame"].append(all_frames[(scan, frame)]['frame'])
                # entries[base + ".frame"].append()
                # println(op, "frm$ctr  ELEMENT  IMAGE $(ctr)")

                # A too-clever-by-half way of optionally live constructing a URL
                if arch is not None:
                    # print(op, "  $comp $prepend_dir")
                    entries[base + ".archive_format"].append(arch)
                    entries[base + ".archive_path"].append(prepend_dir + frame_name)


        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))



    def generate_scan_info(self, cif_block, scan_list):
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



    def generate_step_info(self, cif_block, scan_setting_info, scan_list):
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
                entries[base + ".integration_time"]\
                    .append(scan_setting_info[scan][1]['time'])

        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))



    def generate_array_info(self, cif_block, scan_list):

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
                # println(op, "frm$ctr  ELEMENT  IMAGE $(ctr)")

        for name, value in entries.items():
            cif_block[name] = value

        cif_block.CreateLoop(list(entries.keys()))




    def get_archive_type(self, external_url, frame_dir):

        archive = None

        if external_url == '':
            frame_dir = frame_dir[1:] if frame_dir.startswith('/') else frame_dir
            external_url = "file://" + frame_dir

        else:
            archives = {"TGZ" : r".*((\.tgz\Z)|(\.tar\.gz\Z|))",
                        "TBZ" : r".*((\.tbz\Z)|(\.tar\.bz2\Z|))",
                        "ZIP" : r".*(\.zip\Z)"}

            for archive, regex in archives.items():
                if re.match(regex, external_url):
                    return archive, external_url

        return archive, external_url


    def generate_facility(self, imgblock, source_info):

        regex = r"[^A-Za-z0-9_]"
        base = "_diffrn_source."
        if source_info.get("beamline") is not None or\
             source_info.get("facility") is not None:
            print('facilit', source_info["facility"])
            imgblock[base + "beamline"] = source_info["beamline"]
            imgblock[base + "facility"] = source_info["facility"]
            audit_id_1 = re.sub(regex, '_', source_info["beamline"])
            audit_id_2 = re.sub(regex, '_', source_info["facility"])
        else:
            imgblock[base + "make"] = source_info["manufacturer"] + "-" + \
                source_info["model"]
            imgblock[base + "details"] = \
                f"Located at {source_info['location']}"
            audit_id_1 = re.sub(regex, '_', source_info["manufacturer"])
            audit_id_2 = re.sub(regex, '_', source_info["model"])

        imgblock['_audit.block_code'] = audit_id_1 + '_' + audit_id_2.strip('_')





    def generate_array(self, imgblock, array_info):

        """
        describe_detector(raw_info)

        Produce the information required for array_structure_list. Here
        we assume a rectangular detector with x horizontal, y vertical
        """

        hor, vert = [float(element) for element in array_info['pixel_size']]
        # array structure list axis
        base = "_array_structure_list_axis."
        # TODO always detx dety?
        imgblock[base + "axis_id"] = array_info['axis_id']
        imgblock[base + "axis_set_id"] = array_info['axis_set_id']
        imgblock[base + "displacement"] = [hor / 2, vert / 2]   #half pixel size
        imgblock[base + "displacement_increment"] = [hor, vert]
        imgblock.CreateLoop([
            base + "axis_id",
            base + "axis_set_id",
            base + "displacement",
            base + "displacement_increment"])

        # array structure list
        base = "_array_structure_list."
        imgblock[base + "array_id"] = array_info['array_id']
        imgblock[base + "index"] = array_info['array_index']
        imgblock[base + "axis_set_id"] = array_info['axis_set_id']
        imgblock[base + "dimension"] = array_info['array_dimension']# ['none', 'none']   #number of elements in each direction
        imgblock[base + "direction"] = array_info['array_direction']
        imgblock[base + "precedece"] = array_info['array_precedence']

        print('arryinfo', array_info)

        imgblock.CreateLoop([
            base + "array_id",
            base + "index",
            base + "axis_set_id",
            base + "dimension",
            base + "direction",
            base + "precedece"])


    def generate_detector(self, imgblock, detector_info):

        # len detector axes
        # detector axes
        # rename to det1 etc

        print('detinf', detector_info)

        base = "_diffrn_detector."
        imgblock[base + "id"] = detector_info["detector_id"]
        imgblock[base + "number_of_axes"] = detector_info['number_of_axes']
        imgblock.CreateLoop([base + "id", base + "number_of_axes"])
        #
        base = "_diffrn_detector_axis."
        # print('olla', imgblock["_diffrn_detector.number_of_axes"])
        imgblock[base + "axis_id"] = detector_info['axis_id']
        imgblock[base + "detector_id"] = detector_info['detector_axis_id']
        imgblock.CreateLoop([base + "axis_id", base + "detector_id"])


    def generate_axes(self, imgblock, axes_info):
        """
        describe_axes(raw_info,imgblock)

        Create the goniometer axes corresponding to the data in `raw_info`,
        placing the result in CIF block `imgblock`.
        """


        # print("gonaxes2", gon_axes)
        base = "_axis."
        imgblock[base + "id"] = axes_info['axes']
        imgblock[base + "type"] = axes_info['axis_type']
        imgblock[base + "equipment"] = axes_info['equip']
        imgblock[base + "depends_on"] = axes_info['depends_on']
        user_input.split_vector(axes_info['vector'], base + "vector", imgblock)
        user_input.split_vector(axes_info['offset'], base + "offset", imgblock)
        imgblock.CreateLoop([base + "id", base + "type", base + "equipment",
                            base + "depends_on",
                            base + "vector[1]", base + "vector[2]", base + "vector[3]",
                            base + "offset[1]", base + "offset[2]", base + "offset[3]"])
        print("nnt")
        # print(imgblock.items())


        # [[1, 0, 0], [-1, 0, 0], [1, 0, 0], [[-1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]]], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [[0, 0, 0], [0, 0, 0], [None, None, 0], [0, 0, 0]]])
        # [[1, 0, 0], [-1, 0, 0], [1, 0, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [missing, missing, 0], [0, 0, 0]])
