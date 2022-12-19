import h5py as h5
import numpy as np
from collections import defaultdict
from imgCIF_Creator.output_assembler import imgCIF_assembler
from imgCIF_Creator.information_extractors import user_input



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
            "Fast direction" : setup_info["Fast direction"],
            # "Other detector axes" : "No response",
            # "Principal axis orientation" : "270",
            # "chi" : "180",
            # "kappa" : "No response",
            # "Model" : "Stadivari",
            # "xds.inp" : "No response",
            # "Location" : "Benemerita Uiversidad Autonoma de Puebla",
            # "Name of manufacturer" : "Stoe",
            # "Image orientation" : "top left",
            "Pixel size" : setup_info["Pixel size"],
            "Facility name" : setup_info['Facility name'],
            "Beamline name" : setup_info['Beamline name'],
            "gonio_axes_nxmx" : goniometer_axes,
            "offsets_nxmx" : offsets,
            "Array dimension" : setup_info['Array dimension'],
            "detector_axes" : setup_info['detector_axes']
            }

    user_input.convert_user_input_to_imgcif(info, cif_block)

    imgCIF_assembler.add_scan_info_to_block(
        scan_info, all_frames, cif_block, new_url, prepend_dir, master_file, 'HDF5'
        )


def get_scan_info(master_file):

    scan_info = {}
    #TODO recnum
    rec_num = 5886687
    # axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
    print('master', master_file)
    frame_files = scan_frames_from_master(master_file)

    #TODO is input?
    goniometer_rot_direction = 'clockwise'
    print(f'Goniometer rotation direction set to: {goniometer_rot_direction}')

    # print('ff', frame_files)

    with h5.File(master_file, 'r') as h5_master:
        for scan_no, scan in frame_files.keys():
            for axis in ['phi', 'chi', 'omega', 'fast', 'slow', 'trans']:
                scan_axis_found = False
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
                    break #TODO sys exit

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
            frame_nums_in_file = [1000, 1000, 1000, 600]
            for i, fn in enumerate(frame_files[(scan_no, scan)]):
                print(fn, frame_nums_in_file[i])

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
                        "rad_type" : rad_type}

            scan_info[scan_no] = (axes_single_frame, scan_details)

    # TODO prune scan info
    # obtain all frames
    all_frames = {}
    for scan_no, scan in frame_files:
        frame_no = 1
        #TODO ensure ordering of files 1,2,3...
        for i, entry in enumerate(frame_files[(scan_no, scan)]) :
            filename, dpath = entry
            # print('filen', filename)
            # print('filen', dpath)
            for j in range(1, frame_nums_in_file[i] + 1):

                all_frames[(scan_no, frame_no)] = \
                    (f'https://zenodo.org/record/{rec_num}/files/{filename}',
                    dpath,
                    j
                    )
                frame_no += 1

    return scan_info, all_frames


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
            for key in group:
                link = group.get(key, getlink=True)
                if isinstance(link, h5.ExternalLink):
                    print(f'  {key} -> {link.filename}/{link.path}')
                    scan_frame_files[(scan_no, scan)].append((link.filename, link.path))
    return scan_frame_files


def get_setup_info(master_file):#, filename, fpath_stem, imgfiles):

    scan_info = {}
    setup_info = {}
    #TODO recnum
    rec_num = 5886687
    axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
    frame_files = scan_frames_from_master(master_file)

    # print('ff', frame_files)


    goniometer_axes = {}
    with h5.File(master_file, 'r') as h5_master:
        # TODO right now this does not work for multiple scans in a file
        for scan_no, scan in frame_files.keys():
            setup_info['Facility name'] = \
                h5_master[f'{scan}/instrument/source/name'][()].decode('utf-8')

            path = f'{scan}/instrument/source/name'
            setup_info['Facility name'] = get_item_hdf5(h5_master, path, True, 'facility name')

            path = f'{scan}/instrument/source/beamline' #TODO can this path in theory exist?
            setup_info['Beamline name'] = get_item_hdf5(h5_master, path, True, 'beamline name')

            # setup_info['Start date'] = \
            #     h5_master[f'{scan}/start_time']
            # TODO check if this always works
            fast_direction = h5_master[
                f'{scan}/instrument/detector/module/fast_pixel_direction'].attrs['vector']
            if fast_direction[0] !=0 and fast_direction[1]==fast_direction[2]==0:
                # print('is hor ')
                setup_info['Fast direction'] = 'horizontal'
            elif fast_direction[1] !=0 and fast_direction[0]==fast_direction[2]==0:
                setup_info['Fast direction'] = 'vertical'

            setup_info['Array dimension'] = \
                h5_master[f'{scan}/instrument/detector/module/data_size'][()]

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

                trafo = h5item.attrs['transformation_type'].decode('utf-8')
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
                    {'dep': depends_on, 'vec': vector, 'trafo': trafo}

                # print('myaxt', axis)

                # sortorder = {"sun":0, "mon":1, "tue":2, "wed":3, "thu":4, "fri":5, "sat":6}

                # TODO is an omega axis always present?
                if axis == 'omega':
                    # two_theta axis is duplicate of omega
                    goniometer_axes['two_theta'] = \
                    {'dep': depends_on, 'vec': vector, 'trafo': trafo}



            print('goooni', goniometer_axes)
            base = f'{scan}/instrument/detector/'

            # beam_center = (h5_master[f'{base}beam_center_x'][()] *
            #             h5_master[f'{base}x_pixel_size'][()],
            #             h5_master[f'{base}beam_center_y'][()] *
            #             h5_master[f'{base}y_pixel_size'][()])

            setup_info['Pixel size'] = (h5_master[f'{base}x_pixel_size'][()],
                        h5_master[f'{base}y_pixel_size'][()])

            # find detector axes
            detector_axes = 0
            instrument_group = h5_master[f'{scan}/instrument/']
            for sub_grp in instrument_group.keys():
                try:
                    nx_class = instrument_group[sub_grp].attrs["NX_class"].decode('utf-8')
                    if "NXpositioner" in nx_class:
                        detector_axes += 1
                except KeyError:
                    pass
            setup_info['detector_axes'] = [detector_axes]

        # conversion from NeXus/McStas to CBF coordinate system
        # the kind of transformations (rotation) which have to be performed depend
        # on the position and rotation direction of the goniometer (seen from
        # the source)
        # TODO is this correct?
        if goniometer_axes['omega']['vec'].all() == np.array([-1, 0, 0]).all():
            goniometer_pos = 'right'
        elif goniometer_axes['omega']['vec'].all() == np.array([1, 0, 0]).all():
            goniometer_pos = 'right'

        print(f'Identified goiniometer position (from source): {goniometer_pos}')

        offsets = {}
        for axis in goniometer_axes:
            offsets[axis] = {'vec' : np.array([0.0, 0.0, 0.0])}
        offsets['detx'] = { 'vec' : h5_master[
            f'{scan}/instrument/detector/module/module_offset'].attrs['offset']}
        offsets['trans'] = {'vec' : np.array(
            [0.0, 0.0, h5_master[f'{scan}/instrument/detector_distance'][0]])}

        print('ö+ffis', offsets)

    for axis, content in goniometer_axes.items():
        goniometer_axes[axis]['vec'] = \
            rotate_from_nexus_to_cbf(content['vec'], goniometer_pos)
    for axis, content in offsets.items():
        offsets[axis]['vec'] = \
            rotate_from_nexus_to_cbf(content['vec'], goniometer_pos)

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


def get_item_hdf5(h5_master, path, required=False, message=''):

    found = False
    if required:
        # while not found:
        h5item = h5_master.get(path)
        if h5item is None:
            # found = False
            print(f"Could not find item at {path}, but is required.")
            print(f"Please enter the value for {message}:")
            value = input(' >> ')
        # else:
        #     found = True
    else:
        h5item = h5_master.get(path)
        if h5item is None:
            print(f"Could not find item at {path}, skipping.")
            value = None

    if h5item is not None:
        value = h5item[()].decode('utf-8')

    return value


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





    # try



    return location