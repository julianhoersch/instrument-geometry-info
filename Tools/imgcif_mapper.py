'''
Conversion of metadata information from various input formats
[WIP: currently HDF5 NxMx]
to imgCIF output
'''

from argparse import ArgumentParser
from gettext import translation
from glob import glob
from tkinter.ttk import Separator
import h5py as h5
import os
import re
import requests
import sys
import numpy as np

TOKEN_NAME = 'fda_zen_01'
ACCESS_TOKEN = '73LIoratcpkJ9uRfRUTY4dOiZBWZMP8BoUfUyrcA2ySlIjBH3kCEGYPybhoI'

def get_files_online(rec_num):

    r = requests.get(f'https://zenodo.org/api/records/{rec_num}',
        params={'access_token': ACCESS_TOKEN})

    file_list = []
    for item in r.json()['files']:
        file_link = item['links']['self']
        file_name = os.path.split(file_link)[-1]
        file_list.append((file_link, file_name))

    return file_list


def download_master(file_url, file_name):

    r = requests.get(file_url, params={'access_token': ACCESS_TOKEN})

    with open(file_name, 'wb') as f:
        f.write(r.content)


def frames_from_master(file_name):

    frame_files = []
    with h5.File(file_name, 'r') as f:
        group = f['entry/data']
        for key in group:
            link = group.get(key, getlink=True)
            if isinstance(link, h5.ExternalLink):
                print(f'  {key} -> {link.filename}/{link.path}')
                frame_files.append((link.filename, link.path))
    return frame_files


class CBF_ARRAY_INFO:

    '''
    For CBF data: class to map data array properties
    from a miniCBF header
    '''

    def __init__(self):

        self.items = {
            'order_tag': 'BYTE_ORDER',
            'order_src': 'X-Binary-Element-Byte-Order',
            'compr_tag': 'BYTE_COMPRN',
            'compr_src': 'conversions=',
            'encod_tag': 'BYTE_ENCODE',
            'encod_src': 'Content-Transfer-Encoding'
        }

    def extract_from_header(self, fn):
        with open(fn, 'rb') as f:
            raw_header = f.read(4096).decode('latin-1')
            header_lines = raw_header.splitlines()
        tags_dict = {}
        for ln in header_lines:
            if self.items['order_src'] in ln:
                tags_dict[self.items['order_tag']] = ln.split()[-1]
            if self.items['encod_src'] in ln:
                tags_dict[self.items['encod_tag']] = ln.split()[-1]
            if self.items['compr_src'] in ln:
                tags_dict[self.items['compr_tag']] = ln.split('=')[-1]
        return tags_dict


    def set_void_tags(self):
        '''
        construct imgCIF definitions for HDF5-only use w/o content
        note that this block is not written at all currently, the
        method is for 'reserve'
        '''
        return {'BYTE_ORDER': 'n/a', 'BYTE_ENCODE': 'n/a', 'BYTE_COMPRN': 'n/a'}


    def fill_template(self, td):
        """Fill the arrayt info template.

        Args:
            td (dict): A dictionary containing the byte info.

        Returns:
            string : The byte info string.
        """
        BYTE_INFO_TEMPLATE = f"""\

_array_structure_byte_order         {td['BYTE_ORDER']}
_array_structure_compression_type   {td['BYTE_COMPRN']}
_array_structure.encoding_type      {td['BYTE_ENCODE']}
"""
        return BYTE_INFO_TEMPLATE


class DLS_I04_MAP:

    '''
    Class to map from a Diamond Light Source HDF5-NxMx master
    '''

    def __init__(self):

        # as seen from from the tail of the vector, either 'clockwise' or
        # 'counter_clockwise' -> counter_clockwise will lead to negative
        # increments for the scan
        self._goniometer_rot_direction = 'clockwise'

        self.items = {}
        '''
        for many of the required imgCIF items, correponding HDF5 information is
          at the actual dataset value, as scalar or vector
        '''
        self.items['data_arr_info'] = {
           'dims_path': 'entry/instrument/detector/module/data_size',
           'arrdim_tag': ['ARR_DIMN_1', 'ARR_DIMN_2'],
           'xcent_path': 'entry/instrument/detector/beam_center_x',
           'ycent_path': 'entry/instrument/detector/beam_center_y',
           'beamcent_tag': ['DET_BMCENT_X', 'DET_BMCENT_Y'],
           'xpxsz_path': 'entry/instrument/detector/x_pixel_size',
           'ypxsz_path': 'entry/instrument/detector/y_pixel_size',
           'pixsize_tag': ['DET_PXSIZE_X', 'DET_PXSIZE_Y'],
           'halfsize_tag': ['HALF_PXSIZE_X', 'HALF_PXSIZE_Y']
        }
        self.items['instrument_info'] = {
           'radtype_path': 'entry/instrument/source/type',
           'radtype_tag': 'RADN_TYPE',
           'wavelen_path': 'entry/instrument/beam/incident_wavelength',
           'wavelen_tag': 'RADN_WAVELEN',
           # for the detector distance we have different source items
           'detdist_path': 'entry/instrument/detector_distance',
           'detdist_tag': 'DET_ZDIST',
           'int_time_tag' : 'INT_TIME',
           'int_time_path' : '/entry/instrument/detector/count_time'
        }
        '''
        for axes (NXpositioner) the information sub-dictionary under the given
          key is structured like this:
        - path to HDF5 endpoint (dataset array plus attributes list)
        - content of dataset object:
          default/dummy value for rotation angle (single-element array)
        - attributes list index for 'depends_on'
          target template tag '<AX>_DEP'
        - attributes list index for 'transformation_type'
          target template tag '<AX>_TYPE')
        - attributes list index for 'vector'
          target template tag '<AX>_OVEC')
        '''
        self.items['phi_info'] = {
           'path': 'entry/sample/sample_phi/phi',
           'angle_tag': 'PHI_ANGLE',
           'depend_attr_itemno': 0,
           'depend_tag': 'PHI_DEP',
           'tftype_attr_itemno': 1,
           'tftype_tag': 'PHI_TYPE',
           'vector_attr_itemno': 3,
           'vector_tag': 'PHI_OVEC'
        }
        self.items['chi_info'] = {
           'path': 'entry/sample/sample_chi/chi',
           'angle_tag': 'CHI_ANGLE',
           'depend_attr_itemno': 0,
           'depend_tag': 'CHI_DEP',
           'tftype_attr_itemno': 1,
           'tftype_tag': 'CHI_TYPE',
           'vector_attr_itemno': 3,
           'vector_tag': 'CHI_OVEC'
        }
        self.items['omega_info'] = {
           'path': 'entry/sample/sample_omega/omega',
           'angle_tag': 'OMG_ANGLE',
           'depend_attr_itemno': 0,
           'depend_tag': 'OMG_DEP',
           'tftype_attr_itemno': 1,
           'tftype_tag': 'OMG_TYPE',
           'vector_attr_itemno': 3,
           'vector_tag': 'OMG_OVEC'
        }
        self.items['scan_info'] = {
           'nframes_tag': 'N_FRAMES',
           'nframes_value': 0,
           'axisid_tag': 'SCAN_AXIS',
           'axisid_name': '',
           'axstart_tag': 'SCAN_START',
           'axstart_value': '',
           'axincr_tag': 'SCAN_INCR',
           'axincr_value': '',
           'axrange_tag': 'SCAN_RANGE',
           'axrange_value': ''
        }
        self.items['fast_info'] = {
            'path': '/entry/instrument/detector/module/fast_pixel_direction',
            'angle_tag': 'FST_ANGLE',
            'depend_attr_itemno': 0,
            'depend_tag': 'FST_DEP',
            'offset_attr_itemno': 1,
            'offset_tag': 'FST_OFF',
            'tftype_attr_itemno': 2,
            'tftype_tag': 'FST_TYPE',
            'vector_attr_itemno': 4,
            'vector_tag': 'FST_OVEC'
        }
        self.items['slow_info'] = {
            'path': '/entry/instrument/detector/module/slow_pixel_direction',
            'angle_tag': 'SLW_ANGLE',
            'depend_attr_itemno': 0,
            'depend_tag': 'SLW_DEP',
            'offset_attr_itemno': 1,
            'offset_tag': 'SLW_OFF',
            'tftype_attr_itemno': 2,
            'tftype_tag': 'SLW_TYPE',
            'vector_attr_itemno': 4,
            'vector_tag': 'SLW_OVEC'
        }
        self.items['trans_info'] = {
            'path': '/entry/instrument/detector_z/det_z',
            'angle_tag': 'SLW_ANGLE',
            'depend_attr_itemno': 0,
            'depend_tag': 'TRS_DEP',
            'tftype_attr_itemno': 1,
            'tftype_tag': 'TRS_TYPE',
            'vector_attr_itemno': 3,
            'vector_tag': 'TRS_OVEC'
        }


    def extract_from_hdf5(self, fn, fpath_stem, imgfiles):

        tags_dict = {}
        # scan_axes = {'phi': False, 'chi': False, 'omega': False }

        with h5.File(fn) as f:
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

                '''
                the orientation vector attribute is an array of floats that we convert to string
                and concatenate
                '''
                tag = self.items[f'{k}_info']['vector_tag']
                attr = self.items[f'{k}_info']['vector_attr_itemno']
                separator = '  '
                tags_dict[tag] = separator.join([str(s) for s in list(h5item.attrs.items())[attr][1]])
                # print(tags_dict[tag])

            print('Collect INFO for detector/frame/array')
            sdict = self.items['data_arr_info']
            for tag, item in [(sdict['arrdim_tag'][0], f[sdict['dims_path']][()][0]),
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


    def fill_frame_info_cbf(self, n_frames, imgfiles):

        '''
        Verify that the number of external files matches
        the specified array size as in the metadata.
        '''
        if n_frames != len(imgfiles):
            print(' WARNING!\n'
                  ' N(scan frames) in master file differs from N(image files).\n'
                  ' Frame links will be truncated to available information.'
                  )
        '''
        Construct the required three loops related to scan frames,
        data/image files and the binary IDs that map the frames to files.
        '''
        tags_dict = {}
        tags_dict['ARRAY_COLUMNS_SPEC'] = '_array_data.external_path'
        frame_links = ''
        frame_ids = ''
        scan_frames = ''
        for i, fnpath in enumerate(imgfiles):
            array_links += f'        {(i+1):<4} 1 ext{(i+1):<4}\n'
            frame_links += f'        ext{(i+1):<4} CBF file://{fnpath}\n'
            frame_ids   += f'        {(i+1):<4}  {(i+1):<4} 1\n'
            scan_frames += f'        {(i+1):<4}  SCAN1 {(i+1):<4}\n'
        tags_dict['ARRAY_DATA_INFO'] = array_links
        tags_dict['DATA_EXT_LINKS'] = frame_links
        tags_dict['DATA_FRAME_IDS'] = frame_ids
        tags_dict['SCAN_FRAME_IDS'] = scan_frames

        return tags_dict

    def fill_frame_info_hdf5(self, rec_num, imgfiles, n_frames_per_file):

        tags_dict = {}
        tags_dict['ARRAY_COLUMNS_SPEC'] = \
      '    _array_data_external_data.id\n'\
      '    _array_data_external_data.format\n'\
      '    _array_data_external_data.uri\n'\
      '    _array_data_external_data.path\n'\
      '    _array_data_external_data.frame'
        array_links = ''
        frame_links = ''
        frame_ids = ''
        scan_frames = ''
        k = 1
        for i, entry in enumerate(imgfiles):
            fn, dpath = entry
            for j in range(1, n_frames_per_file[i]+1):
                array_links += f'        {k:<4} 1 ext{k:<4}\n'
                frame_links += f'        ext{k:<4} HDF5 https://zenodo.org/record/{rec_num}/files/{fn} {dpath} {j}\n'
                frame_ids   += f'        {k:<4}  {k:<4} 1\n'
                scan_frames += f'        {k:<4}  SCAN1 {k:<4}\n'
                k += 1
        tags_dict['ARRAY_DATA_INFO'] = array_links
        tags_dict['DATA_EXT_LINKS'] = frame_links
        tags_dict['DATA_FRAME_IDS'] = frame_ids
        tags_dict['SCAN_FRAME_IDS'] = scan_frames

        return tags_dict

    def find_max_len(self, col, dic):
        """Find the maximum string lenght of a column to align entries propertly.

        Args:
            col (string): The string that identifies entries belonging to the
                same column (from tags)..
            dic (dict): The dictionary containing all tags.

        Returns:
            max_len (int) : The maximum lenght.
        """
        max_len = 0
        for tag, value in dic.items():
            if 'RADN_TYPE' in tag:
                continue
            if col in tag:
                if len(value) > max_len:
                    max_len = len(value)

        return max_len


    def fill_template(self, td):
        """Fill the template for DLS_I04 with the values in the tags dictionary.

        Note: The depends_on attributes of the two_theta and dety axis are altered
        manually in the extract from hdf5 method.

        Args:
            td (dict): The tags dictionary.

        Returns:
            DLS_I04_TEMPLATE (string): The template filled with values.
        """

        stop = {}
        for col in ['TYPE', 'DEP', 'OVEC', 'OFFSET']:
            stop[col[:2]] = self.find_max_len(col, td)

        DLS_I04_TEMPLATE = f"""\
data_{td['RECORD']}
_database.dataset_doi       '10.5281/zenodo.5886687'
_audit.block_code           Diamond_I04
_diffrn_source.beamline     I04
_diffrn_source.facility     Diamond
_diffrn_radiation.type      '{td['RADN_TYPE']}'
{td['BYTE_INFO_BLOCK']}
loop_
    _diffrn_radiation_wavelength.id
    _diffrn_radiation_wavelength.value
        1          {td['RADN_WAVELEN']}

loop_
    _axis.id
    _axis.type
    _axis.equipment
    _axis.depends_on
    _axis.vector[1]
    _axis.vector[2]
    _axis.vector[3]
    _axis.offset[1]
    _axis.offset[2]
    _axis.offset[3]
        phi        {td['PHI_TYPE']:{stop['TY']}}  goniometer  {td['PHI_DEP']:{stop['DE']}}  {td['PHI_OVEC']:{stop['OV']}}  {td['PHI_OFFSET']}
        chi        {td['CHI_TYPE']:{stop['TY']}}  goniometer  {td['CHI_DEP']:{stop['DE']}}  {td['CHI_OVEC']:{stop['OV']}}  {td['CHI_OFFSET']}
        omega      {td['OMG_TYPE']:{stop['TY']}}  goniometer  {td['OMG_DEP']:{stop['DE']}}  {td['OMG_OVEC']:{stop['OV']}}  {td['OMG_OFFSET']}
        two_theta  {td['OMG_TYPE']:{stop['TY']}}  detector    {td['OMG_DEP']:{stop['DE']}}  {td['OMG_OVEC']:{stop['OV']}}  {td['OMG_OFFSET']}
        trans      {td['TRS_TYPE']:{stop['TY']}}  detector    {td['TRS_DEP']:{stop['DE']}}  {td['TRS_OVEC']:{stop['OV']}}  {td['TRS_OFFSET']}
        detx       {td['FST_TYPE']:{stop['TY']}}  detector    {td['FST_DEP']:{stop['DE']}}  {td['FST_OVEC']:{stop['OV']}}  {td['FST_OFFSET']}
        dety       {td['SLW_TYPE']:{stop['TY']}}  detector    {td['SLW_DEP']:{stop['DE']}}  {td['SLW_OVEC']:{stop['OV']}}  {td['SLW_OFFSET']}

loop_
    _array_structure_list_axis.axis_id
    _array_structure_list_axis.axis_set_id
    _array_structure_list_axis.displacement
    _array_structure_list_axis.displacement_increment
        detx        1         {td['HALF_PXSIZE_X']}        {td['DET_PXSIZE_X']}
        dety        2         {td['HALF_PXSIZE_Y']}        {td['DET_PXSIZE_Y']}

loop_
    _array_structure_list.array_id
    _array_structure_list.axis_set_id
    _array_structure_list.direction
    _array_structure_list.index
    _array_structure_list.precedence
    _array_structure_list.dimension
        1          1          increasing    1        1       {td['ARR_DIMN_1']}
        1          2          increasing    2        2       {td['ARR_DIMN_2']}

_diffrn_detector.id               det1
_diffrn_detector.number_of_axes   1

loop_
    _diffrn_detector_axis.axis_id
    _diffrn_detector_axis.detector_id
        trans      det1

loop_
    _array_data.binary_id
    _array_data.array_id
    _array_data.external_data_id
{td['ARRAY_DATA_INFO']}
loop_
{td['ARRAY_COLUMNS_SPEC']}
{td['DATA_EXT_LINKS']}
loop_
    _diffrn_data_frame.id
    _diffrn_data_frame.binary_id
    _diffrn_data_frame.array_id
{td['DATA_FRAME_IDS']}
_diffrn_scan.id                               SCAN1
_diffrn_scan.frames                           {td['N_FRAMES']}
_diffrn_scan.integration_time 	              {td['INT_TIME']}
_diffrn_scan_axis.axis_id                     {td['SCAN_AXIS']}
_diffrn_scan_axis.angle_start                 {td['SCAN_START']}
_diffrn_scan_axis.angle_range                 {td['SCAN_RANGE']}
_diffrn_scan_axis.angle_increment             {td['SCAN_INCR']}
_diffrn_scan_axis.displacement_start          0.0
_diffrn_scan_axis.displacement_range          0.0
_diffrn_scan_axis.displacement_increment      0.0

loop_
    _diffrn_scan_frame.frame_id
    _diffrn_scan_frame.scan_id
    _diffrn_scan_frame.frame_number
{td['SCAN_FRAME_IDS']}"""

        return DLS_I04_TEMPLATE

    def write_imgcif(self, tags):

        with open('cbf_metadata.cif', 'w') as f:
            f.write(tags)


def check_frames_location(fpath_stem):
    '''
    Only applicable in case of a local folder with file set
    '''
    fnames = sorted(glob(f'{fpath_stem}*.cbf'))
    print(f'found N = {len(fnames)} files under the path/stem name')
    if len(fnames) == 0:
        print(' WARNING!\n'
              ' No files found matching the pattern;'
              ' Check correctness of path and stem name'
             )
    return fnames


def transform_nexus_to_cbf_coordinates(tags_dict, goniometer_pos,
                                       goniometer_rot_direction):
    """Transform the from NeXus to CBF coordinate system depending on the
    position and rotation direction of the goniometer.

    Args:
        tags_dict (dict) : The dictionary containing the tags and its values
        goniometer_pos (string): The position of the goniometer seen from the
            source ('left' or 'right')
        goniometer_rot_direction (string): The rotation direction of the
            goniometer seen from the tail of the vector ('clockwise' or
            'counter_clockwise')

    Returns:
        tags_dict (dict): The transformed tags dictionary.
    """
    for tag, value in tags_dict.items():
        if 'OVEC' in tag or 'OFFSET' in tag:
            tags_dict[tag] = rotate_from_nexus_to_cbf(value, goniometer_pos)
        elif 'SCAN_INCR' in tag:
            assert goniometer_rot_direction in ['clockwise', 'counter_clockwise'], \
                """Unknown goniometer position! Please choose between
                'clockwise' and 'counter_clockwise.'"""
            if goniometer_rot_direction == 'counter_clockwise':
                tags_dict[tag] = f'-{value}'
                tags_dict['SCAN_RANGE'] = f"-{tags_dict['SCAN_RANGE']}"

    return tags_dict


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
    single_separator = ' '
    vector_array = np.fromstring(vector, dtype=np.float64, sep=single_separator)
    count = int(vector.count(single_separator) / (len(vector_array) - 1))
    separator = single_separator*count

    if rotation_axis == 'x':
        rot = vector_array * np.array([1, -1, -1])
    elif rotation_axis == 'y':
        rot = vector_array * np.array([-1, 1, -1])
    elif rotation_axis == 'z':
        rot = vector_array * np.array([-1, -1, 1])
    else:
        rot = None
    # remove leading minus at zeros
    rot[rot==0] = 0.0
    rot = separator.join(str(x) for x in rot)

    return rot


def main():

    print('imgCIF mapping tool')

    ap = ArgumentParser(prog='imgcif_mapper')
    ap.add_argument('record', type=str, help='record # of the Zenodo repository')
    ap.add_argument('-i', '--input-file', type=str, help='path/name of input HDF5 file with metadata')
    ap.add_argument('-f', '--frames', type=str, help='path and name stem of external image frames')
    ap.add_argument('-a', '--archive', type=str, help='path and name stem of external image frames')
    args = ap.parse_args(sys.argv[1:])

    print('record #:',args.record)

    if args.input_file is None:
        # Default
        print('taking metadata from online HDF5 master')
        file_list = get_files_online(args.record)
        for flink, fname in file_list:
            print(f'{flink:84s} {fname:20s}')
            if 'master' in fname:
                meta_file = fname
                meta_link = flink
        print('master found:', meta_file)
        download_master(meta_link, meta_file)
    else:
        # Overrides online mode if present
        print('metadata source is input file:', args.input_file)
        if not os.path.exists(args.input_file):
            print('local input file not found')
            exit(0)
        meta_file = args.input_file

    data_array_info = CBF_ARRAY_INFO()

    if args.frames is None:
        print('get image frames from master/repository')
        imgfiles = frames_from_master(meta_file)
        byte_tags = data_array_info.set_void_tags()  # not used at all under this condition
        byte_block_str = ''
        mode = 'hdf5'
    else:
        print('external image file location:', '/'.join(args.frames.split('/')[:-1]))
        imgfiles = check_frames_location(args.frames)
        byte_tags = data_array_info.extract_from_header(imgfiles[0])
        byte_block_str = data_array_info.fill_template(byte_tags)
        mode = 'hybrid'

    metadata_map = DLS_I04_MAP()
    metadata_tags = metadata_map.extract_from_hdf5(meta_file, args.frames, imgfiles)

    if mode == 'hdf5':
        n_frames = metadata_map.items['scan_info']['nframes_value']
        print(f'{n_frames} frames in {len(imgfiles)} HDF5 data files.')
        print('Please give a list of N(frames) per file (info missing in the master)')
        user_str = input(' > ')
        for pattern in [' ', ',', ', ']:
            try:
                frame_nums = [int(x) for x in user_str.split(pattern)]
            except:
                continue
        for i, fn in enumerate(imgfiles):
            print(fn, frame_nums[i])
        framelink_tags = metadata_map.fill_frame_info_hdf5(args.record, imgfiles, frame_nums)
    else:
        framelink_tags = metadata_map.fill_frame_info_cbf(len(imgfiles), imgfiles)

    tags = {'RECORD': args.record, 'BYTE_INFO_BLOCK': byte_block_str, **metadata_tags, **framelink_tags}
    template_filled = metadata_map.fill_template(tags)


    metadata_map.write_imgcif(template_filled)


if __name__ == '__main__':

    main()

