import h5py as h5
import numpy as np
from collections import defaultdict
import imgCIF_Creator.assembler.imgCIF_assembler as imgCIF_assembler
import imgCIF_Creator.information_extractors.user_input_extractor as input_extractor




# def add_scan_info_to_block(scan_info, all_frames, cif_block, new_url,
#                            prepend_dir, directory):

def extract_hdf5_NxMx_scan(
    master_file, cif_block, file_format, new_url=[], prepend_dir=''):

    scan_info, all_frames = get_scan_info(master_file)
    goniometer_axes, offsets = get_setup_info(master_file)

    # Rename axes
    # if len(axis) > 0:
    #     axis_to_rename = {axis[0].lower() : axis[1].lower()}
    #     # print("lo axes ", axis_to_rename)
    #     rename_axes(scan_info, axis_to_rename)

    # if include_archive_directory:
    #     prepend_dir = os.path.split(directory)[-1]
    # else:
    #     prepend_dir = ""


    user_input = {
            "Data DOI" : "10.5281/zenodo.7155191",
            "Comments" : "No response",
            "Goniometer axes" : "Phi, a, Chi, c, Omega, a",
            "Two theta axis" : "a",
            "Fast direction" : "horizontal",
            "Other detector axes" : "No response",
            "Principal axis orientation" : "270",
            "chi" : "180",
            "kappa" : "No response",
            "Model" : "Stadivari",
            "xds.inp" : "No response",
            "Location" : "Benemerita Uiversidad Autonoma de Puebla",
            "Name of manufacturer" : "Stoe",
            "Image orientation" : "top left",
            "Pixel size" : "0.172",
            "gonio_axes_nxmx" : goniometer_axes,
            "offsets_nxmx" : offsets,
            }

    input_extractor.convert_user_input_to_imgcif(user_input, cif_block)

    # if output is not N:
    imgCIF_assembler.add_scan_info_to_block(
        scan_info, all_frames, cif_block, new_url, prepend_dir, master_file, file_format
        )


def get_scan_info(master_file):

    scan_info = {}
    #TODO recnum
    rec_num = 5886687
    axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
    frame_files = scan_frames_from_master(master_file)

    goniometer_rot_direction = 'clockwise'
    print(f'Goniometer rotation direction set to: {goniometer_rot_direction}')


    print('ff', frame_files)

    with h5.File(master_file, 'r') as h5_master:
        for scan_no, scan in frame_files.keys():
            for axis in ['phi', 'chi', 'omega', 'fast', 'slow', 'trans']:

                path = f'{scan}/sample/sample_{axis}/{axis}'
                print('hdff5 path', path)
                h5item = h5_master.get(path)
                if h5item is not None:
                    print('Collect INFO for', axis)
                    try:
                        print('thisitem', axis, h5item[()])
                        if len(h5item[()]) > 1:
                            print(' ... identified as scan axis')
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

            # for scan in s     can_names.keys()[0]:

            # print('unic scan', scan)

            # scan_names = list(filter(lambda x : x[0] == scan, all_frames.keys()))
            # print('scnn', scan_names)
            # frames = [x[1] for x in scan_names]


            # Get information for first
            # file_name = os.path.join(frame_dir, all_frames[(scan, 1)])
            # vals, scan_ax, scan_incr, exposure, wl = get_frame_info(file_name, axes)

            axes_single_frame = {}

            for axis in ['phi', 'chi', 'omega', 'trans']:
                path = f'{scan}/sample/sample_{axis}/{axis}'
                if axis == 'trans':
                    path = f'{scan}/instrument/detector_z/det_z'

                # print('h5item', h5_master[path][()])
                # we take the last as this was done in the julia cbf code also
                # TODO convert units if neccessary
                axes_single_frame[axis] =  h5_master[path][()][0]

            # print('axsingfr', axes_single_frame)

            # return 1, 1

            # get_frame_info(frame_type, file_name, axes)
            # print("vals ", vals, "scanax ", scan_ax, "scaninc ", scan_incr, "exposure ", exposure, "wl ", wl)
            # start = vals[scan_ax]

            # Get information and last frame
            # file_name = os.path.join(frame_dir, all_frames[(scan, len(frames))])
            # vals, _, _, _, _ = get_frame_info(file_name, axes)
            # finish = vals[scan_ax]

            # Check increment and range match

            if not np.isclose(scan_start + scan_incr * (n_frames - 1), scan_stop, atol=1e-6):
                raise Exception(
                    f"Scan range does not match increment: \
{scan_start} to {scan_stop}, {n_frames} steps of {scan_incr}")

            exposure = h5_master[f'{scan}/instrument/detector/count_time'][()]
            wl = h5_master[f'{scan}/instrument/beam/incident_wavelength'][()]
            scan_details = {"frames" : n_frames,
                        "axis" : scan_axis,
                        "incr" : scan_incr,
                        "time" : exposure,
                        "start" : scan_start,
                        "range" : scan_range,
                        "wavelength" : wl}

            scan_info[scan_no] = (axes_single_frame, scan_details)

    # TODO prune scan info


    # print('printiprint')
    # frame_files = []
    # with h5.File(master_file, 'r') as h5_master:
    #     group = h5_master['entry/data']
    #     print('master', h5_master.keys())
    #     print('keykey', h5_master[list(h5_master.keys())[0]].keys())


    # for i, entry in enumerate(imgfiles):
    #     fn, dpath = entry
    #     for j in range(1, n_frames_per_file[i]+1):
    #          frame_links += f'        ext{k:<4} HDF5 https://zenodo.org/record/{rec_num}/files/{fn} {dpath} {j}\n'


    # return None, None

    # obtain all frames
    all_frames = {}

    # scan_frame_files[(scan_no, scan)].append((link.filename, link.path))
    for scan_no, scan in frame_files:
        frame_no = 1
        #TODO ensure ordering of files 1,2,3...
        for i, entry in enumerate(frame_files[(scan_no, scan)]) :
            filename, dpath = entry
            print('filen', filename)
            print('filen', dpath)
            for j in range(1, frame_nums_in_file[i] + 1):
                # array_links += f'        {k:<4} 1 ext{k:<4}\n'
                # frame_links += f'        ext{k:<4} HDF5 https://zenodo.org/record/{rec_num}/files/{fn} {dpath} {j}\n'
                # frame_ids += f'        {k:<4}  {k:<4} 1\n'
                # scan_frames += f'        {k:<4}  SCAN1 {k:<4}\n'

                all_frames[(scan_no, frame_no)] = \
                    (f'https://zenodo.org/record/{rec_num}/files/{filename}',
                    dpath,
                    j
                    )
                #TODO  {dpath} {j}\n'
                frame_no += 1

    # scan_info = {}
    # scan_info[scan] = (vals, details)

    # print("this is scanifo", scan_info)
    # print("this is all frames", list(all_frames)[-100:])
    # print("this is all frames", list(all_frames.values())[-100:])

    return scan_info, all_frames


def scan_frames_from_master(file_name):
    """Obtain the frames from the master.

    Args:
        file_name (str): the file name

    Returns:
        list: a list of tuples containing the frame files
    """

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


# def get_frame_info():


#     ax_vals = \
#     list(map(lambda ax : (ax.lower(), get_CBF_header_values(lines, ax)),
#                 axes))




def get_setup_info(master_file):#, filename, fpath_stem, imgfiles):

    scan_info = {}
    #TODO recnum
    rec_num = 5886687
    axes = imgCIF_assembler.ROT_AXES + imgCIF_assembler.TRANS_AXES
    frame_files = scan_frames_from_master(master_file)

    print('ff', frame_files)


    setup_info = {
        'layout': 'beamline',
        'Facility name': '',
        'Beamline name': 'p08',
        'Is updated?': 'no',
        'Start date': '1.1.1',
        'Finish date': '2.2.2',
        'Principal axis orientation':'270',
        'Goniometer axes': "Phi, a, Chi, c, Omega, a",
        'Rotation axis name': 'chi',
        'kappa': 'kappa 60',
        'Image orientation':'top left',
        'Fast direction': 'top right',
        'Pixel size': '0.12',
        'Detector axes': '',
        'xds.inp': '',
        'Data DOI': '',
        'Comments': ''}

    goniometer_axes = {}
    with h5.File(master_file, 'r') as h5_master:
        for scan_no, scan in frame_files.keys():
            setup_info['Beamline name'] = \
                h5_master[f'{scan}/instrument/source/name'].attrs['short_name']
            # TODO loop through hdf5?
            # for axis in ['phi', 'chi', 'omega']:
            path = f'{scan}/sample/transformations' #sample_{axis}/{axis}'
            h5item = h5_master.get(path)
            for axis in list(h5item.keys()):

                path = f'{scan}/sample/transformations/{axis}' #sample_{axis}/{axis}'
                print('hdff5 path', path)
                h5item = h5_master.get(path)
                # print('hdf5', list(h5item.keys()))
                if h5item is None or len(h5item.attrs) == 0:
                    continue
                else:
                    print('Collect INFO for', axis)


    # the transformation 'type' attribute is a bytestring tb. coverted to string
                # tag = self.items[f'{k}_info']['tftype_tag']
                # attr = self.items[f'{k}_info']['tftype_attr_itemno']
                # tags_dict[tag] = list(h5item.attrs.items())[attr][1].decode('utf-8')

                # attr = axis.upper() + '_TYPE'
                print('hdfattr', dict(h5item.attrs.items()))
                trafo = h5item.attrs['transformation_type'] #.items())[3][1]#.decode('utf-8')
                print('ttf', trafo)

                # same for the 'depends_on' attribute, but for imgCIF we need
                # only the path tail
                # tag = self.items[f'{k}_info']['depend_tag']
                # attr = self.items[f'{k}_info']['depend_attr_itemno']

                # iterate through the dependencies until a depndency on
                # '.', 'omega', 'chi' or 'phi' is found
                # depends_on = list(h5item.attrs.items())\
                #     [0][1].decode('utf-8').split('/')[-1]
                # depends_on_h5item = h5_master[list(h5item.attrs.items())\
                #     [0][1].decode('utf-8')]

                depends_on = h5item.attrs['depends_on'].decode('utf-8').split('/')[-1]
                depends_on_h5item = h5_master[h5item.attrs['depends_on']]

                print('depo', depends_on)

                while depends_on in ['sam_x', 'sam_y', 'sam_z', 'module_offset']:
                    depends_on = depends_on_h5item.attrs['depends_on']\
                        .decode('utf-8').split('/')[-1]
                    # print('depo while', depends_on)

                    depends_on_h5item = h5_master[depends_on_h5item.attrs['depends_on']]
                    # depends_on_h5item = h5_master[list(depends_on_h5item.attrs.items())\
                    #     [0][1].decode('utf-8')]

                # rename det_z to trans
                if depends_on == 'det_z':
                    depends_on = 'trans'

                # some manual changes in the depends_on fields
                if axis == 'sam_x': #'slow':
                    depends_on = 'detx'
                elif axis == 'trans':
                    depends_on = 'two_theta'

                # tags_dict[tag] = depends_on

                # the orientation vector attribute is an array of floats that
                # we convert to string and concatenate
                # tag = self.items[f'{k}_info']['vector_tag']
                # attr = self.items[f'{k}_info']['vector_attr_itemno']
                # separator = '  '
                # tags_dict[tag] = separator.join(
                #     [str(s) for s in list(h5item.attrs.items())[attr][1]])
                # # print(tags_dict[tag])

                vector = h5item.attrs['vector']#.items())[3][1]
                print('vecdorr', vector)
                print('dep', depends_on)
                goniometer_axes[axis] = \
                    {'dep': depends_on, 'vec': vector, 'trafo': trafo}

            print('goooni', goniometer_axes)
            # goniometer_axes_sorted = sorted(goniometer_axes, key=)

            # recursive_sort = lambda goniometer_axes['.'][0]
            # \

            goniometer_axes_sorted = []
            def recursive_sort(axis):

                if goniometer_axes[axis]['dep'] == '.':
                    goniometer_axes_sorted

                for axis, attrs in goniometer_axes.items():
                    if attrs['dep'] == '.':
                        return

                # goniometer_axes_sorted.append(axis)
                goniometer_axes.keys()   [axis]['depends_on']
                if goniometer_axes[axis]['depends_on'] in goniometer_axes:
                    dependent_axis = goniometer_axes[axis]['axis']
                    goniometer_axes_sorted.append(dependent_axis)


                    return recursive_sort(dependent_axis)
                # else:
                #     goniometer_axes_sorted.append()

            # for axis, attrs in goniometer_axes.items():
            #     if attrs['dep'] == '.':
            #         return


            # recursive_sort('.')
            # goniometer_axes_sorted.reverse()

            print('srt', goniometer_axes_sorted)

            base = f'{scan}/instrument/detector/'
            dims = h5_master[f'{base}module/data_size'][()]

            beam_center = (h5_master[f'{base}beam_center_x'][()] *
                        h5_master[f'{base}x_pixel_size'][()],
                        h5_master[f'{base}beam_center_y'][()] *
                        h5_master[f'{base}y_pixel_size'][()])

            setup_info['Pixel Size'] = (h5_master[f'{base}x_pixel_size'][()],
                        h5_master[f'{base}y_pixel_size'][()])

            # print('Collect INFO for detector/frame/array')
            # sdict = self.items['data_arr_info']
            # for tag, item in [
            #         (sdict['arrdim_tag'][0], f[sdict['dims_path']][()][0]),
            #         (sdict['arrdim_tag'][1], f[sdict['dims_path']][()][1]),
            #         (sdict['beamcent_tag'][0],
            #             f[sdict['xcent_path']][()] * f[sdict['xpxsz_path']][()]),
            #         (sdict['beamcent_tag'][1],
            #             f[sdict['ycent_path']][()] * f[sdict['ypxsz_path']][()]),

            #         (sdict['pixsize_tag'][0], f[sdict['xpxsz_path']][()]),
            #         (sdict['pixsize_tag'][1], f[sdict['ypxsz_path']][()]),
            #         (sdict['halfsize_tag'][0], f[sdict['xpxsz_path']][()]/2.0),
            #         (sdict['halfsize_tag'][1], f[sdict['ypxsz_path']][()]/2.0),
            #     ]:
            #     tags_dict[tag] = item


            base = f'{scan}/instrument/'
            rad_type = h5_master[f'{base}source/type'][()].decode('utf-8')
            wavelength = h5_master[f'{base}beam/incident_wavelength'][()]
            detector_distance = h5_master[f'{base}detector_distance'][()]
            integration_time = h5_master[f'{base}detector/count_time'][()]

            # print('Collect INFO for instrument setup')
            # sdict = self.items['instrument_info']
            # for tag, item in [(sdict['radtype_tag'], f[sdict['radtype_path']][()].decode('utf-8')),
            #                     (sdict['wavelen_tag'], f[sdict['wavelen_path']][()]),
            #                     (sdict['detdist_tag'], f[sdict['detdist_path']][()][0]),
            #                     (sdict['int_time_tag'], f[sdict['int_time_path']][()]),
            #                     ]:
            #     tags_dict[tag] = item



        # print('Collect INFO for scan axis')
        # sdict = self.items['scan_info']
        # for tag, item in [(sdict['nframes_tag'], sdict['nframes_value']),
        #                     (sdict['axisid_tag'], sdict['axisid_name']),
        #                     (sdict['axstart_tag'], sdict['axstart_value']),
        #                     (sdict['axrange_tag'], sdict['axrange_value']),
        #                     (sdict['axincr_tag'], sdict['axincr_value'])
        #                     ]:
        #     tags_dict[tag] = item

        # conversion from NeXus/McStas to CBF coordinate system
        # the kind of transformations (rotation) which have to be performed depend
        # on the position and rotation direction of the goniometer (seen from
        # the source)
        # if tags_dict['OMG_OVEC'] == separator.join(['-1.0', '0.0', '0.0']):
        #     self._goniometer_pos = 'right'
        # elif tags_dict['OMG_OVEC'] == separator.join(['1.0', '0.0', '0.0']):
        #     self._goniometer_pos = 'left'
        # else:
        #     self._goniometer_pos = 'undefined'

        # this is omega
        # TODO is this correct?
        print('arr', goniometer_axes['omega']['vec'])
        if goniometer_axes['omega']['vec'].all() == np.array([-1, 0, 0]).all():
            goniometer_pos = 'right'
        elif goniometer_axes['omega']['vec'].all() == np.array([1, 0, 0]).all():
            goniometer_pos = 'right'

        print(f'Identified goiniometer position (from source): {goniometer_pos}')


        # # create offsets
        # for offset in ['PHI_OFFSET', 'CHI_OFFSET', 'OMG_OFFSET', 'SLW_OFFSET']:
        #     tags_dict[offset] = separator.join(['0.0', '0.0', '0.0'])

        # tags_dict['TRS_OFFSET'] = \
        #     separator.join(['0.0', '0.0', f"{tags_dict['DET_ZDIST']}"])
        # tags_dict['FST_OFFSET'] = \
        #     separator.join([f"{tags_dict['DET_BMCENT_X']}",
        #                     f"{tags_dict['DET_BMCENT_Y']}", '0.0'])

        offsets = {}
        for axis in goniometer_axes:
            offsets[axis] = {'vec' : np.array([0.0, 0.0, 0.0])}
        offsets['detx'] = { 'vec' : h5_master[
            f'{scan}/instrument/detector/module/module_offset'].attrs['offset']}
        offsets['trans'] = {'vec' : np.array(
            [0.0, 0.0, h5_master[f'{scan}/instrument/detector_distance'][0]])}

        print('ö+ffis', offsets)

    # offsets[]

    for axis, content in goniometer_axes.items():
        goniometer_axes[axis]['vec'] = \
            rotate_from_nexus_to_cbf(content['vec'], goniometer_pos)
        print('axiiicont', axis)
    for axis, content in offsets.items():
        print('axiii', axis)
        print('goto', content['vec'])
        offsets[axis]['vec'] = \
            rotate_from_nexus_to_cbf(content['vec'], goniometer_pos)

    # for axis in offsets:
    #     axis['vec'] = rotate_from_nexus_to_cbf(content['vec'], goniometer_pos)

    # goniometer_axes_cbf, offsets_cbf = transform_nexus_to_cbf_coordinates(
    #     goniometer_axes, offsets, goniometer_pos)

    return goniometer_axes, offsets


# def transform_nexus_to_cbf_coordinates(axes, goniometer_pos):
#     """Transform the from NeXus to CBF coordinate system depending on the
#     position and rotation direction of the goniometer.

#     #TODO updata doc
#     Args:
#         tags_dict (dict) : The dictionary containing the tags and its values
#         goniometer_pos (string): The position of the goniometer seen from the
#             source ('left' or 'right')
#         goniometer_rot_direction (string): The rotation direction of the
#             goniometer seen from the tail of the vector ('clockwise' or
#             'counter_clockwise')

#     Returns:
#         tags_dict (dict): The transformed tags dictionary.
#     """

#     for axis in axes:
#         rotate_from_nexus_to_cbf(axis['vec'], goniometer_pos)


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
    # single_separator = ' '
    # vector_array = np.fromstring(vector, dtype=np.float64, sep=single_separator)
    # count = int(vector.count(single_separator) / (len(vector_array) - 1))
    # separator = single_separator*count

    if rotation_axis == 'x':
        rot = vector * np.array([1, -1, -1])
    elif rotation_axis == 'y':
        rot = vector * np.array([-1, 1, -1])
    elif rotation_axis == 'z':
        rot = vector * np.array([-1, -1, 1])
    else:
        rot = None
    # remove leading minus at zeros
    # rot[rot == 0] = 0.0
    # rot = separator.join(str(x) for x in rot)

    return rot








def extract_from_hdf5(self, filename, fpath_stem, imgfiles):
    """Extract information from hdf5.

    Args:
        filename (str): the filename from which to extract

    Returns:
        dict: a dictionary containing the tags
    """

    tags_dict = {}
    scan_axes = {'phi': False, 'chi': False, 'omega': False }

    with h5.File(filename) as f:
        for k in ['phi', 'chi', 'omega', 'fast', 'slow', 'trans']:

            print('Collect INFO for', k)
            h5item = f[self.items[f'{k}_info']['path']]
            try:
                if len(h5item[()]) > 1:
                    print(' ... identified as scan axis')
                    self.items['scan_info']['axisid_name'] = k
                    n_frames = len(h5item[()])
                    self.items['scan_info']['nframes_value'] = n_frames
                    self.items['scan_info']['axstart_value'] = f'{h5item[()][0]}'
                    incr_range = f[f'entry/sample/transformations/{k}_increment_set'][()]
                    self.items['scan_info']['axincr_value'] = f'{incr_range[0]}'
                    axis_angle_range = n_frames * incr_range[0]        # total range
                    self.items['scan_info']['axrange_value'] = f'{axis_angle_range}'
            except TypeError:
                # for fast and slow this is a scalar and no list
                pass

            # the transformation 'type' attribute is a bytestring tb. coverted to string
            tag = self.items[f'{k}_info']['tftype_tag']
            attr = self.items[f'{k}_info']['tftype_attr_itemno']
            tags_dict[tag] = list(h5item.attrs.items())[attr][1].decode('utf-8')

            # same for the 'depends_on' attribute, but for imgCIF we need only the path tail
            tag = self.items[f'{k}_info']['depend_tag']
            attr = self.items[f'{k}_info']['depend_attr_itemno']

            # iterate through the dependencies until a depndency on
            # '.', 'omega', 'chi' or 'phi' is found
            depends_on = list(h5item.attrs.items())\
                [attr][1].decode('utf-8').split('/')[-1]
            depends_on_h5item = f[list(h5item.attrs.items())\
                [attr][1].decode('utf-8')]

            while depends_on in ['sam_x', 'sam_y', 'sam_z', 'module_offset']:
                depends_on = list(depends_on_h5item.attrs.items())\
                    [attr][1].decode('utf-8').split('/')[-1]

                depends_on_h5item = f[list(depends_on_h5item.attrs.items())\
                    [attr][1].decode('utf-8')]

            # rename det_z to trans
            if depends_on == 'det_z':
                depends_on = 'trans'

            # some manual changes in the depends_on fields
            if k == 'slow':
                depends_on = 'detx'
            elif k == 'trans':
                depends_on = 'two_theta'

            tags_dict[tag] = depends_on


            # the orientation vector attribute is an array of floats that
            # we convert to string and concatenate
            tag = self.items[f'{k}_info']['vector_tag']
            attr = self.items[f'{k}_info']['vector_attr_itemno']
            separator = '  '
            tags_dict[tag] = separator.join(
                [str(s) for s in list(h5item.attrs.items())[attr][1]])
            # print(tags_dict[tag])

        print('Collect INFO for detector/frame/array')
        sdict = self.items['data_arr_info']
        for tag, item in [
                (sdict['arrdim_tag'][0], f[sdict['dims_path']][()][0]),
                (sdict['arrdim_tag'][1], f[sdict['dims_path']][()][1]),
                (sdict['beamcent_tag'][0],
                    f[sdict['xcent_path']][()] * f[sdict['xpxsz_path']][()]),
                (sdict['beamcent_tag'][1],
                    f[sdict['ycent_path']][()] * f[sdict['ypxsz_path']][()]),
                (sdict['pixsize_tag'][0], f[sdict['xpxsz_path']][()]),
                (sdict['pixsize_tag'][1], f[sdict['ypxsz_path']][()]),
                (sdict['halfsize_tag'][0], f[sdict['xpxsz_path']][()]/2.0),
                (sdict['halfsize_tag'][1], f[sdict['ypxsz_path']][()]/2.0),
            ]:
            tags_dict[tag] = item

        print('Collect INFO for instrument setup')
        sdict = self.items['instrument_info']
        for tag, item in [(sdict['radtype_tag'], f[sdict['radtype_path']][()].decode('utf-8')),
                            (sdict['wavelen_tag'], f[sdict['wavelen_path']][()]),
                            (sdict['detdist_tag'], f[sdict['detdist_path']][()][0]),
                            (sdict['int_time_tag'], f[sdict['int_time_path']][()]),
                            ]:
            tags_dict[tag] = item

    print('Collect INFO for scan axis')
    sdict = self.items['scan_info']
    for tag, item in [(sdict['nframes_tag'], sdict['nframes_value']),
                        (sdict['axisid_tag'], sdict['axisid_name']),
                        (sdict['axstart_tag'], sdict['axstart_value']),
                        (sdict['axrange_tag'], sdict['axrange_value']),
                        (sdict['axincr_tag'], sdict['axincr_value'])
                        ]:
        tags_dict[tag] = item

    # conversion from NeXus/McStas to CBF coordinate system
    # the kind of transformations (rotation) which have to be performed depend
    # on the position and rotation direction of the goniometer (seen from
    # the source)
    if tags_dict['OMG_OVEC'] == separator.join(['-1.0', '0.0', '0.0']):
        self._goniometer_pos = 'right'
    elif tags_dict['OMG_OVEC'] == separator.join(['1.0', '0.0', '0.0']):
        self._goniometer_pos = 'left'
    else:
        self._goniometer_pos = 'undefined'

    print(f'Identified goiniometer position (from source): {self._goniometer_pos}')
    print(f'Goniometer rotation direction set to: {self._goniometer_rot_direction}')

    # create offsets
    for offset in ['PHI_OFFSET', 'CHI_OFFSET', 'OMG_OFFSET', 'SLW_OFFSET']:
        tags_dict[offset] = separator.join(['0.0', '0.0', '0.0'])

    tags_dict['TRS_OFFSET'] = \
        separator.join(['0.0', '0.0', f"{tags_dict['DET_ZDIST']}"])
    tags_dict['FST_OFFSET'] = \
        separator.join([f"{tags_dict['DET_BMCENT_X']}",
                        f"{tags_dict['DET_BMCENT_Y']}", '0.0'])

    tags_dict = transform_nexus_to_cbf_coordinates(
        tags_dict, self._goniometer_pos, self._goniometer_rot_direction)

    return tags_dict