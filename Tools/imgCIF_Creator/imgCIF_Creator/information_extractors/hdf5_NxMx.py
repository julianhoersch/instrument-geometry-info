import sys
import h5py as h5
import numpy as np
from collections import defaultdict
from imgCIF_Creator.output_assembler import imgCIF_assembler
from imgCIF_Creator.information_extractors import user_input
from imgCIF_Creator.command_line_interfaces import parser

from . import extractor_interface

class extractor(extractor_interface.ExtractorInterface):

    def __init__(self, filename, _) -> None:
        super().__init__()

        # TODO design choice: put scan info into the class? if we want to check
        # later that the number ov frames in all files from the user input is
        # correct this makes it easier re request the input
        frame_files, frames_per_file = scan_frames_from_master(filename)
        self.scan_info = get_scan_info(filename, frame_files)
        self.all_frames = get_all_frames(frame_files, frames_per_file)
        # print('all frames', self.all_frames)

        self.goniometer_axes, self.offsets, self.setup_info = \
            get_setup_info(filename, frame_files)

        self._first_scan = sorted(self.scan_info.keys())[0]
        print('frst sc', self._first_scan)


    def get_miscellaneous_info(self):

        doi = None
        overload = self.setup_info.get('overload')

        return {'doi' : doi,
                'overload' : overload,
                }


    def get_source_info(self):

        facility = self.setup_info['facility']
        beamline = self.setup_info['beamline']
        manufacturer = None
        model = self.setup_info['detector_model']
        location = None

        source_info = {'beamline' : beamline,
                       'facility' : facility,
                       'manufacturer' : manufacturer,
                       'model' : model,
                       'location' : location}

        return source_info


    def get_axes_info(self):

        axes, axis_type, equip, depends_on, vector, offset = \
            make_axes(self.goniometer_axes, self.offsets)

        axes_info = {'axes' : axes,
                     'axis_type' : axis_type,
                     'equip' : equip,
                     'depends_on' : depends_on,
                     'vector' : vector,
                     'offset' : offset}

        return axes_info


    def get_array_info(self):

        pixel_size = \
            (self.setup_info['x_pixel_size'], self.setup_info['y_pixel_size'])

        if self.setup_info['fast_direction'] == "horizontal":
            array_precedence = [1, 2]
        else:
            array_precedence = [2, 1]


        # can I find any of this information in the h5?
        array_info = {
            'axis_id' : None, # names as array_x, array_y
            'axis_set_id': None, # integers as 1, 2
            'pixel_size' : pixel_size,
            'array_id' : None,
            'array_index' : None,
            'array_dimension' : self.setup_info['array_dimensions'],
            'array_direction' : self.setup_info['fast_direction'],
            'array_precedence' : array_precedence,
        }

        return array_info


    def get_detector_info(self):

        #TODO what is detector axes? can I extract it from the hdf5 with
        # the number of nxpositioners? >> should be the number of detector axes
        # in the gonio settings, is hardcoded...

        #TODO multiple detectors are not supportet (yet?)
        #TODO nothing from hdf5?


        detector_info = {
            'detector_id' : None,
            'number_of_axes' : None,
            'axis_id' : None,
            'detector_axis_id' : None
        }

        return detector_info


    def get_wavelength_info(self):

        # TODO what do we actually do with additional information??
        # TODO not creating a list of wl id's here (yet)
        # TODO could infer rad type from wl...

        wl = self.scan_info[self._first_scan][1]['wavelength']
        rad_type = self.scan_info[self._first_scan][1]['rad_type']

        return {'rad_type' : rad_type,
                'wavelength' : wl}


    def get_scan_settings_info(self):

        # is already in a dict format
        return self.scan_info



def extract_hdf5_NxMx_scan(
    master_file, cif_block, new_url=[], prepend_dir=''):

    scan_info, all_frames = get_scan_info(master_file)
    goniometer_axes, offsets, setup_info = get_setup_info(master_file)

    # Rename axes
    # if len(axis) > 0:
    #     axis_to_rename = {axis[0].lower() : axis[1].lower()}
    #     # print("lo axes ", axis_to_rename)
    #     rename_axes(scan_info, axis_to_rename)

    # if include_archive_directory:
    #     prepend_dir = os.path.split(directory)[-1]
    # else:
    #     prepend_dir = ""


    info = {
            # "Data DOI" : "10.5281/zenodo.7155191",
            # "Comments" : "No response",
            # "Goniometer axes" : "Phi, a, Chi, c, Omega, a",
            # "Two theta axis" : "a",
            "fast_direction" : setup_info["fast_direction"],
            # "Other detector axes" : "No response",
            # "Principal axis orientation" : "270",
            # "chi" : "180",
            # "kappa" : "No response",
            # "Model" : "Stadivari",
            # "xds.inp" : "No response",
            # "Location" : "Benemerita Uiversidad Autonoma de Puebla",
            # "Name of manufacturer" : "Stoe",
            # "Image orientation" : "top left",
            "x_pixel_size" : setup_info["x_pixel_size"],
            "y_pixel_size" : setup_info["y_pixel_size"],
            "facility" : setup_info['facility'],
            "beamline" : setup_info['beamline'],
            "gonio_axes_nxmx" : goniometer_axes,
            "offsets_nxmx" : offsets,
            "array_dimensions" : setup_info['array_dimensions'],
            "number_of_detector_axes" : setup_info['number_of_detector_axes']
            }

    user_input.convert_user_input_to_imgcif(info, cif_block)

    imgCIF_assembler.add_scan_info_to_block(
        scan_info, all_frames, cif_block, new_url, prepend_dir, master_file, 'HDF5'
        )


def get_scan_info(master_file, frame_files):



    scan_info = {}
    #TODO recnum
    # rec_num = 5886687
    # axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
    # print('master', master_file)


    #TODO is input?
    # goniometer_rot_direction = 'clockwise'
    # print(f'Goniometer rotation direction set to: {goniometer_rot_direction}')
    goniometer_rot_direction =  \
        parser.CommandLineParser().request_input('goniometer_rot_direction')

    print('ff', frame_files)

    with h5.File(master_file, 'r') as h5_master:

         for scan_no, scan in frame_files.keys():
            scan_axis_found = False
            for axis in ['phi', 'chi', 'omega', 'fast', 'slow', 'trans']:
                path = f'{scan}/sample/sample_{axis}/{axis}'
                print('hdff5 path', path)
                h5item = h5_master.get(path)
                if h5item is not None:
                    print('Collect INFO for', axis)
                    scan_axis_found = True
                    try:
                        print('thisitem', axis, h5item[()])
                        if len(h5item[()]) > 1:
                            print(' ... identified as scan axis')
                            scan_axis_found = True
                            scan_axis = axis
                            n_frames = len(h5item[()])
                            scan_start = h5item[()][0]
                            scan_stop = h5item[()][-1]
                            print('thisitem', h5item[()])
                            scan_incr = h5_master[
                                f'{scan}/sample/transformations/{scan_axis}_increment_set'
                                ][()][0]
                            # if the rotation direction is the other way round,
                            # chage inrement
                            if goniometer_rot_direction == 'counter_clockwise':
                                scan_incr = scan_incr * -1
                            # total range
                            scan_range = n_frames * scan_incr
                    except TypeError:
                        # for fast and slow this is a scalar and no list
                        pass
            if not scan_axis_found:
                print('No scan axis found! Aborting!')
                sys.exit() #TODO sys exit

            # return 0, 0
            # # if mode == 'hdf5':
            # for scan_no, scan in frame_files.keys():
            # print(f'{n_frames} frames in {len(frame_files[(scan_no, scan)])} HDF5 data files.')
            # print('Please give a list of N(frames) per file (info missing in the master)')
            # user_str = input(' > ')
            # for pattern in [' ', ',', ', ']:
            #     try:
            #         frame_nums = [int(x) for x in user_str.split(pattern)]
            #     except:
            #         continue
            # frame_nums_in_file = [1000, 1000, 1000, 600]
            # for i, fn in enumerate(frame_files[(scan_no, scan)]):
            #     print(fn, frame_nums_in_file[i])

            axes_single_frame = {}
            for axis in ['phi', 'chi', 'omega', 'trans']:
                path = f'{scan}/sample/sample_{axis}/{axis}'
                if axis == 'trans':
                    path = f'{scan}/instrument/detector_z/det_z'

                # we take the last as this was done in the julia cbf code also
                # TODO convert units if neccessary
                # this are the settings of the axes, hence the value in the hdf5 file
                axes_single_frame[axis] =  h5_master[path][()][0]
                # print('thats what i want', h5_master[path][()][0])

            # Check increment and range match
            if not np.isclose(scan_start + scan_incr * (n_frames - 1), scan_stop, atol=1e-6):
                raise Exception(
                    f"Scan range does not match increment: \
{scan_start} to {scan_stop}, {n_frames} steps of {scan_incr}")

            exposure = h5_master[f'{scan}/instrument/detector/count_time'][()]
            wl = h5_master[f'{scan}/instrument/beam/incident_wavelength'][()]
            rad_type = h5_master[f'{scan}/instrument/source/type'][()].decode('utf-8')
            scan_details = {"frames" : n_frames,
                            "axis" : scan_axis,
                            "incr" : scan_incr,
                            "time" : exposure,
                            "start" : scan_start,
                            "range" : scan_range,
                            "wavelength" : wl,
                            "rad_type" : rad_type,
                            }

            scan_info[scan_no] = (axes_single_frame, scan_details)


    return scan_info



def get_all_frames(frame_files, frames_per_file):

    # frame_files, frames_per_file = scan_frames_from_master(master_file)

    # TODO prune scan info
    # obtain all frames
    all_frames = {}
    for scan_no, scan in frame_files:
        #TODO ensure ordering of files 1,2,3...
        counter = 1
        for i, entry in enumerate(frame_files[(scan_no, scan)]) :
            filename, dpath = entry
            # print('filen', filename)
            # print('filen', dpath)
            for frame_in_file in range(1, frames_per_file[i] + 1):
                all_frames[(scan_no, counter)] = {
                    'filename' : filename,
                    'path' : dpath,
                    'frame' : frame_in_file
                }
                counter += 1
                    # (f'https://zenodo.org/record/{rec_num}/files/{filename}',
                    # dpath,
                    # j
                    # )
                # frame_no += 1

    # print('all frames', all_frames)

    return all_frames


def scan_frames_from_master(file_name):
    """Obtain the frames from the master.

    Args:
        file_name (str): the file name

    Returns:
        list: a list of tuples containing the frame files
    """

    print('myfname', file_name)
    scan_frame_files = defaultdict(list)
    with h5.File(file_name, 'r') as h5_master:
        for scan_no, scan in enumerate(h5_master.keys(), 1):
            print(scan_no, scan)
            group = h5_master[f'{scan}/data']
            print('grp', group)
            for key in group:
                link = group.get(key, getlink=True)
                if isinstance(link, h5.ExternalLink):
                    print(f'  {key} -> {link.filename}/{link.path}')
                    scan_frame_files[(scan_no, scan)].append((link.filename, link.path))

        # if mode == 'hdf5':
        # for scan_no, scan in frame_files.keys():
        # print(f'{n_frames} frames in {len(frame_files[(scan_no, scan)])} HDF5 data files.')
        print('sccfr', scan_frame_files)
        first_key = list(scan_frame_files.keys())[0]
        n_files = len(scan_frame_files[first_key])
        frame_numbers = []
        while len(frame_numbers) != n_files:
            print(f'\nFound {n_files} files in the hdf5 master file.', end='')
            frame_numbers = parser.CommandLineParser().request_input('frame_numbers')
            frame_numbers = frame_numbers.replace(' ', '').split(',')
            frame_numbers = [int(number) for number in frame_numbers]
            print('frame nums', frame_numbers)

        # for pattern in [' ', ',', ', ']:
        #     try:
        #         frame_nums = [int(x) for x in user_str.split(pattern)]
        #     except:
        #         continue
        # # frame_nums_in_file = [1000, 1000, 1000, 600]

    return scan_frame_files, frame_numbers


def get_setup_info(master_file, frame_files):#, filename, fpath_stem, imgfiles):

    # scan_info = {}
    setup_info = {}
    #TODO recnum
    # rec_num = 5886687
    # axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
    # frame_files = scan_frames_from_master(master_file)

    # print('ff', frame_files)


    goniometer_axes = {}
    with h5.File(master_file, 'r') as h5_master:
        # TODO right now this does not work for multiple scans in a file
        for scan_no, scan in frame_files.keys():

            path = f'{scan}/instrument/source/name'
            setup_info['facility'] = get_hdf5_item(h5_master, path)

            path = f'{scan}/instrument/source/beamline' #TODO can this path in theory exist?
            setup_info['beamline'] = get_hdf5_item(h5_master, path)

            path = f'{scan}/entry/instrument/detector/description'
            setup_info['detector_model'] = get_hdf5_item(h5_master, path)

            path = f'{scan}/instrument/detector/saturation_value'
            setup_info['overload'] = get_hdf5_item(h5_master, path)

            path = f'{scan}/instrument/detector/x_pixel_size'
            setup_info['x_pixel_size'] = get_hdf5_item(h5_master, path)

            path = f'{scan}/instrument/detector/y_pixel_size'
            setup_info['y_pixel_size'] = get_hdf5_item(h5_master, path)

            path = f'{scan}/instrument/detector/module/data_size'
            setup_info['array_dimensions'] = get_hdf5_item(h5_master, path)

            # path = f'{scan}/instrument/detector/depends_on'
            # setup_info['detector_axis_id'] = get_hdf5_item(h5_master, path)

            path = f'{scan}/instrument/detector/module/fast_pixel_direction'
            fast_direction = get_hdf5_item(h5_master, path, 'vector')

            if fast_direction[0] !=0 and fast_direction[1]==fast_direction[2]==0:
                # print('is hor ')
                setup_info['fast_direction'] = 'horizontal'
            elif fast_direction[1] !=0 and fast_direction[0]==fast_direction[2]==0:
                setup_info['fast_direction'] = 'vertical'

            # setup_info['Start date'] = \
            #     h5_master[f'{scan}/start_time']
            # # TODO check if this always works
            # fast_direction = h5_master[
            #     f'{scan}/instrument/detector/module/fast_pixel_direction'].attrs['vector']



            # h5_master[f'{scan}/instrument/detector/module/data_size'][()]

            # TODO loop through hdf5?
            # for axis in ['phi', 'chi', 'omega']:
            path = f'{scan}/sample/transformations' #sample_{axis}/{axis}'
            h5item = h5_master.get(path)
            for axis in list(h5item.keys()):

                path = f'{scan}/sample/transformations/{axis}' #sample_{axis}/{axis}'
                axis = axis.lower()
                # print('hdff5 path', path)
                h5item = h5_master.get(path)
                # print('hdf5', list(h5item.keys()))
                if h5item is None or len(h5item.attrs) == 0:
                    continue
                else:
                    print('Collect INFO for', axis)

                axis_type = h5item.attrs['transformation_type'].decode('utf-8')
                depends_on = h5item.attrs['depends_on'].decode('utf-8').split('/')[-1]
                depends_on_h5item = h5_master[h5item.attrs['depends_on']]

                # while depends_on in ['sam_x', 'sam_y', 'sam_z', 'module_offset']:
                #     depends_on = depends_on_h5item.attrs['depends_on']\
                #         .decode('utf-8').split('/')[-1]
                    # print('depo while', depends_on)

                    # depends_on_h5item = h5_master[depends_on_h5item.attrs['depends_on']]
                    # depends_on_h5item = h5_master[list(depends_on_h5item.attrs.items())\
                    #     [0][1].decode('utf-8')]

                # the dependecy chain starts with detx for chi
                # rename det_z to trans
                if depends_on == 'sam_z':
                    depends_on = 'trans'
                elif depends_on == 'sam_x':
                    depends_on = 'detx'
                elif depends_on == 'sam_y':
                    depends_on = 'dety'

                # some manual changes in the depends_on fields
                if axis == 'sam_x': #'slow':
                    axis = 'detx'
                    # depends_on = 'detx'
                elif axis == 'sam_y':
                    axis = 'dety'
                elif axis == 'sam_z':
                    axis = 'trans'
                    # depends_on = 'two_theta'

                vector = h5item.attrs['vector']#.items())[3][1]

                goniometer_axes[axis] = \
                    {
                     'depends_on': depends_on,
                     'vector': vector,
                     'axis_type': axis_type
                     }

                # print('myaxt', axis)

                # sortorder = {"sun":0, "mon":1, "tue":2, "wed":3, "thu":4, "fri":5, "sat":6}

                # TODO is an omega axis always present?
                if axis == 'omega':
                    # two_theta axis is duplicate of omega
                    goniometer_axes['two_theta'] = \
                        {
                         'depends_on': depends_on,
                         'vector': vector,
                         'axis_type': axis_type
                         }

            print('goooni', goniometer_axes)
            # base = f'{scan}/instrument/detector/'

            # beam_center = (h5_master[f'{base}beam_center_x'][()] *
            #             h5_master[f'{base}x_pixel_size'][()],
            #             h5_master[f'{base}beam_center_y'][()] *
            #             h5_master[f'{base}y_pixel_size'][()])


            # setup_info['pixel_size'] = (h5_master[f'{base}x_pixel_size'][()],
            #             h5_master[f'{base}y_pixel_size'][()])

            # find detector axes
            number_of_detector_axes = 0
            instrument_group = h5_master[f'{scan}/instrument/']
            for sub_grp in instrument_group.keys():
                try:
                    nx_class = instrument_group[sub_grp].attrs["NX_class"].decode('utf-8')
                    if "NXpositioner" in nx_class:
                        number_of_detector_axes += 1
                except KeyError:
                    pass
            setup_info['number_of_detector_axes'] = number_of_detector_axes

        # conversion from NeXus/McStas to CBF coordinate system
        # the kind of transformations (rotation) which have to be performed depend
        # on the position and rotation direction of the goniometer (seen from
        # the source)
        # TODO is this correct?
        if goniometer_axes['omega']['vector'].all() == np.array([-1, 0, 0]).all():
            goniometer_pos = 'right'
        elif goniometer_axes['omega']['vector'].all() == np.array([1, 0, 0]).all():
            goniometer_pos = 'right'

        print(f'Identified goiniometer position (from source): {goniometer_pos}')

        offsets = {}
        for axis in goniometer_axes:
            offsets[axis] = {'vector' : np.array([0.0, 0.0, 0.0])}

        path = f'{scan}/instrument/detector/module/module_offset'
        offsets['detx'] = {'vector' : get_hdf5_item(h5_master, path, 'offset')}

        path = f'{scan}/instrument/detector_distance'
        offsets['trans'] = {'vector' :
            np.array([0.0, 0.0, get_hdf5_item(h5_master, path)])}

        print('ö+ffis', offsets)

    for axis, content in goniometer_axes.items():
        goniometer_axes[axis]['vector'] = \
            rotate_from_nexus_to_cbf(content['vector'], goniometer_pos)
    for axis, content in offsets.items():
        offsets[axis]['vector'] = \
            rotate_from_nexus_to_cbf(content['vector'], goniometer_pos)

    return goniometer_axes, offsets, setup_info



def rotate_from_nexus_to_cbf(vector, goniometer_pos):
    """Rotate the vectors from NeXus to CBF coordinate system.

    Args:
        vector (string): The vector that shall be rotated as string e.g. '1 0 0'
        goniometer_pos (string): The position of the goniometer seen from the
            source ('left' or 'right')

    Returns:
        rotated_vector (string): The rotated vector.
    """

    assert goniometer_pos in ['left', 'right'], \
        "Unknown goniometer position! Please choose between 'left' and 'right'."

    if goniometer_pos == 'right':
        rotated_vector = rotate_vector(vector, 'y')

    elif goniometer_pos == 'left':
        rotated_vector = rotate_vector(vector, 'x')

    return rotated_vector


def rotate_vector(vector, rotation_axis):
    """Rotate a vector around a rotation axis.

    Args:
        vector (string): The vector that shall be rotated as string e.g. '1 0 0'
        rotation_axis (string): The axis around which shall be rotated
            'x', 'y' or 'z')

    Returns:
        rotated_vector (string): The rotated vector.
    """
    print('vector to rotate', vector)

    if rotation_axis == 'x':
        rot = vector * np.array([1, -1, -1])
    elif rotation_axis == 'y':
        rot = vector * np.array([-1, 1, -1])
    elif rotation_axis == 'z':
        rot = vector * np.array([-1, -1, 1])
    else:
        rot = None

    return rot



            #  goniometer_axes_sorted = []
            # def recursive_sort(axis):

            #     if goniometer_axes[axis]['dep'] == '.':
            #         goniometer_axes_sorted

            #     for axis, attrs in goniometer_axes.items():
            #         if attrs['dep'] == '.':
            #             return

            #     # goniometer_axes_sorted.append(axis)
            #     goniometer_axes.keys()   [axis]['depends_on']
            #     if goniometer_axes[axis]['depends_on'] in goniometer_axes:
            #         dependent_axis = goniometer_axes[axis]['axis']
            #         goniometer_axes_sorted.append(dependent_axis)


            #         return recursive_sort(dependent_axis)


def get_hdf5_item(h5_master, path, attr=None):

    h5item = h5_master.get(path)
    if h5item is not None:
        if attr is not None:
            h5item = h5item.attrs[attr]
        h5item = h5item[()]
        try:
            h5item = h5item.decode('utf-8')
        except AttributeError:
            pass

    return h5item


# def get_hdf5_item(h5_master, path, required=False, message=''):

#     found = False
#     if required:
#         # while not found:
#         h5item = h5_master.get(path)
#         if h5item is None:
#             # found = False
#             print(f"Could not find item at {path}, but is required.")
#             print(f"Please enter the value for {message}:")
#             value = input(' >> ')
#         # else:
#         #     found = True
#     else:
#         h5item = h5_master.get(path)
#         if h5item is None:
#             print(f"Could not find item at {path}, skipping.")
#             value = None

#     if h5item is not None:
#         value = h5item[()].decode('utf-8')

#     return value


def get_attribute_hdf5(h5_master, path, attr, required=False, message=''):

    found = False
    if required:
        while not found:
            h5item = h5_master.get(path)

        if h5item is None:
            found = False
        else:
            if attr is not None:
                try:
                    attribute = h5_master[path].attrs[attr]
                    found = True
                except KeyError:
                    if required:
                        print(f"Could not find attribute {attr} in {path}, but is required.")
                        print(f"Please enter the right location for {message}:")
                        path = input(' >> ')
            else:
                return None
    else:
        pass # TODO

    return attribute



def make_axes(goniometer_axes, offsets):

    axes = []
    axis_type = []
    equip = []
    depends_on = []
    vector = []
    offset = []

    print('raw egg', goniometer_axes)
    # add two theta somewhered

    for axis in goniometer_axes:
        axes.append(axis)
        ax_type = goniometer_axes[axis]['axis_type']
        axis_type.append(ax_type)
        #TODO can i do that?
        if ax_type == 'rotation':
            equip.append('goniometer')
        elif ax_type == 'translation':
            equip.append('detector')
        else:
            equip.append('equipment')
        depends_on.append(goniometer_axes[axis]['depends_on'])
        vector.append(list(goniometer_axes[axis]['vector']))
        offset.append(list(offsets[axis]['vector']))

    print('offs', offset)

    return [axes, axis_type, equip, depends_on, vector, offset]