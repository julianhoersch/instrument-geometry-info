'''
Conversion of metadata information from various input formats
[WIP: currently HDF5 NxMx]
to imgCIF output
'''

from argparse import ArgumentParser
from glob import glob
import h5py as h5
import os
import re
import requests
import sys

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

        # construct for HDF5-only use 

        return {'BYTE_ORDER': 'n/a', 'BYTE_ENCODE': 'n/a', 'BYTE_COMPRN': 'n/a'} 
    

DLS_I04_TEMPLATE = """\
_audit.block_id	Diamond_I04
_diffrn_source.beamline	I04
_diffrn_source.facility	Diamond

    _array_structure_byte_order         %(BYTE_ORDER)s
    _array_structure_compression_type   %(BYTE_COMPRN)s
    _array_structure.encoding_type      %(BYTE_ENCODE)s

   _diffrn_radiation.type     '%(RADN_TYPE)s'
 
 loop_
      _diffrn_radiation_wavelength.id
      _diffrn_radiation_wavelength.value
       1          %(RADN_WAVELEN)s

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
         phi        %(PHI_TYPE)s  goniometer  %(PHI_DEP)s  %(PHI_OVEC)s  0  0  0
         chi        %(CHI_TYPE)s  goniometer  %(CHI_DEP)s  %(CHI_OVEC)s  0  0  0
         omega      %(OMG_TYPE)s  goniometer  %(OMG_DEP)s  %(OMG_OVEC)s  0  0  0
         two_theta  rotation     detector    .          ?   0  0   0  0  0
         trans      translation  detector    two_theta  0   0  1   0  0  %(DET_ZDIST)s
         detx       translation  detector    trans      1   0  0  -%(DET_BMCENT_X)s  0  0
         dety       translation  detector    trans      0  -1  0   0  %(DET_BMCENT_Y)s  0

    loop_
      _array_structure_list_axis.axis_id
      _array_structure_list_axis.axis_set_id
      _array_structure_list_axis.displacement
      _array_structure_list_axis.displacement_increment
         detx                    1                    0                  %(DET_PXSIZE_X)s
         dety                    2                    0                  %(DET_PXSIZE_Y)s

    loop_
      _array_structure_list.array_id
      _array_structure_list.axis_set_id
      _array_structure_list.direction
      _array_structure_list.index
      _array_structure_list.precedence
      _array_structure_list.dimension
         1             1             increasing             1             1       %(ARR_DIMN_1)s
         1             2             increasing             2             2       %(ARR_DIMN_2)s

    loop_
      _diffrn_detector.id
      _diffrn_detector.number_of_axes
         det1                        1

    loop_
      _diffrn_detector_axis.axis_id
      _diffrn_detector_axis.detector_id
         trans                    det1

    loop_
      _array_data.array_id
      _array_data.external_format
      %(MODE_DEPENDENT_INFO)s
      _array_data.external_path
      _array_data.external_location_uri
      _array_data.external_archive_path
      _array_data.external_frame
%(DATA_EXT_LINKS)s

    loop_
      _diffrn_data_frame.id
      _diffrn_data_frame.binary_id
      _diffrn_data_frame.array_id
%(DATA_FRAME_IDS)s

    _diffrn_scan.id SCAN1
    _diffrn_scan.frames                      %(N_FRAMES)s
    _diffrn_scan_axis.axis_id                %(SCAN_AXIS)s
    _diffrn_scan_axis.angle_start            %(SCAN_START)s
    _diffrn_scan_axis.displacement_start     0
    _diffrn_scan_axis.angle_increment        %(SCAN_INCR)s

    loop_
      _diffrn_scan_frame.frame_id
      _diffrn_scan_frame.scan_id
      _diffrn_scan_frame.frame_number
%(SCAN_FRAME_IDS)s      
"""

class DLS_I04_MAP:

    '''
    Class to map from a Diamond Light Source HDF5-NxMx master 
    '''

    def __init__(self):

        self.template = DLS_I04_TEMPLATE

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
           'pixsize_tag': ['DET_PXSIZE_X', 'DET_PXSIZE_Y']
        }
        self.items['instrument_info'] = {
           'radtype_path': 'entry/instrument/source/type',
           'radtype_tag': 'RADN_TYPE',
           'wavelen_path': 'entry/instrument/beam/incident_wavelength',
           'wavelen_tag': 'RADN_WAVELEN',
           # for the detector distance we have different source items
           'detdist_path': 'entry/instrument/detector_distance',
           'detdist_tag': 'DET_ZDIST'
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
           'axincr_value': ''
        }

    def extract_from_hdf5(self, fn, fpath_stem, imgfiles):

        tags_dict = {}
        scan_axes = {'phi': False, 'chi': False, 'omega': False }

        with h5.File(fn) as f:
            for k in ['phi', 'chi', 'omega']:

                print('Collect INFO for', k)
                h5item = f[self.items[f'{k}_info']['path']]
                if len(h5item[()]) > 1:
                    print(' ... identified as scan axis')
                    self.items['scan_info']['axisid_name'] = k
                    n_frames = len(h5item[()])
                    self.items['scan_info']['nframes_value'] = n_frames
                    self.items['scan_info']['axstart_value'] = f'{h5item[()][0]}'
                    incr_range = f[f'entry/sample/transformations/{k}_increment_set'][()]
                    self.items['scan_info']['axincr_value'] = f'{incr_range[0]}'


                # the transformation 'type' attribute is a bytestring tb. coverted to string
                tag = self.items[f'{k}_info']['tftype_tag']
                attr = self.items[f'{k}_info']['tftype_attr_itemno']
                tags_dict[tag] = list(h5item.attrs.items())[attr][1].decode('utf-8')

                # same for the 'depends_on' attribute, but for imgCIF we need only the path tail
                tag = self.items[f'{k}_info']['depend_tag']
                attr = self.items[f'{k}_info']['depend_attr_itemno']
                tags_dict[tag] = list(h5item.attrs.items())[attr][1].decode('utf-8').split('/')[-1]

                '''
                the orientation vector attribute is an array of floats that we convert to string
                and concatenate
                '''
                tag = self.items[f'{k}_info']['vector_tag']
                attr = self.items[f'{k}_info']['vector_attr_itemno']
                tags_dict[tag] = ' '.join([str(s) for s in list(h5item.attrs.items())[attr][1]])

            print('Collect INFO for detector/frame/array')
            sdict = self.items['data_arr_info']
            for tag, item in [(sdict['arrdim_tag'][0], f[sdict['dims_path']][()][0]),
                              (sdict['arrdim_tag'][1], f[sdict['dims_path']][()][1]),
                              (sdict['beamcent_tag'][0], f[sdict['xcent_path']][()]),
                              (sdict['beamcent_tag'][1], f[sdict['ycent_path']][()]),
                              (sdict['pixsize_tag'][0], f[sdict['xpxsz_path']][()]),
                              (sdict['pixsize_tag'][1], f[sdict['ypxsz_path']][()]),
                             ]:
                tags_dict[tag] = item

            print('Collect INFO for instrument setup')
            sdict = self.items['instrument_info']
            for tag, item in [(sdict['radtype_tag'], f[sdict['radtype_path']][()].decode('utf-8')),
                              (sdict['wavelen_tag'], f[sdict['wavelen_path']][()]),
                              (sdict['detdist_tag'], f[sdict['detdist_path']][()][0])
                             ]:
                tags_dict[tag] = item

        print('Collect INFO for scan axis')
        sdict = self.items['scan_info']
        for tag, item in [(sdict['nframes_tag'], sdict['nframes_value']),
                          (sdict['axisid_tag'], sdict['axisid_name']),
                          (sdict['axstart_tag'], sdict['axstart_value']),
                          (sdict['axincr_tag'], sdict['axincr_value'])
                         ]:
            tags_dict[tag] = item
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
        tags_dict['MODE_DEPENDENT_INFO'] = '_array_data.external_path'
        frame_links = ''
        frame_ids = ''
        scan_frames = '' 
        for i, fnpath in enumerate(imgfiles):
            frame_links += f'        ext{(i+1):<4} CBF file://{fnpath}\n'
            frame_ids   += f'        {(i+1):4}  ext{(i+1):<4} 1\n'
            scan_frames += f'        {(i+1):4}  SCAN1 {(i+1):4}\n'
        tags_dict['DATA_EXT_LINKS'] = frame_links
        tags_dict['DATA_FRAME_IDS'] = frame_ids
        tags_dict['SCAN_FRAME_IDS'] = scan_frames

        return tags_dict

    def fill_frame_info_hdf5(self, rec_num, imgfiles, n_frames_per_file):

        tags_dict = {}
        tags_dict['MODE_DEPENDENT_INFO'] = '_array_data.external_location_uri\n
_array_data.external_archive_path\n
_array_data.external_frame'
        frame_links = ''
        frame_ids = ''
        scan_frames = '' 
        k = 1
        for i, entry in enumerate(imgfiles):
            fn, dpath = entry
            for j in range(1, n_frames_per_file[i]+1):
                frame_links += f'        ext{k:<4} HDF5 https://zenodo.org/api/records/{rec_num}/{fn} {dpath} {j}\n'
                frame_ids   += f'        {k:4}  ext{k:<4} 1\n'
                scan_frames += f'        {k:4}  SCAN1 {k:4}\n'
                k += 1
        tags_dict['DATA_EXT_LINKS'] = frame_links
        tags_dict['DATA_FRAME_IDS'] = frame_ids
        tags_dict['SCAN_FRAME_IDS'] = scan_frames

        return tags_dict
        
    def write_imgcif(self, tags):

        with open('cbf_metadata.cif', 'w') as f:
            f.write(self.template % tags)


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
        byte_tags = data_array_info.set_void_tags()
        mode = 'hdf5'
    else:
        print('external image file location:', '/'.join(args.frames.split('/')[:-1]))
        imgfiles = check_frames_location(args.frames)
        byte_tags = data_array_info.extract_from_header(imgfiles[0])
        mode = 'hybrid'
    
    metadata_map = DLS_I04_MAP()
    metadata_tags = metadata_map.extract_from_hdf5(meta_file, args.frames, imgfiles)

    if args.frames is None:
        n_frames = metadata_map.items['scan_info']['nframes_value']
        print(f'{n_frames} frames in {len(imgfiles)} HDF5 data files.')
        print('Please give a list of frames per file (info not available from the master)')
        user_str = input(' > ')
        for pattern in [' ', ',', ', ']:
            try:
                frame_nums = [int(x) for x in user_str.split(pattern)]
            except:
                continue
        for i, fn in enumerate(imgfiles):
            print(fn, frame_nums[i])
        #framelink_tags = metadata_map.fill_frame_info_cbf(n_frames, imgfiles)
        framelink_tags = metadata_map.fill_frame_info_hdf5(args.record, imgfiles, frame_nums)
    else:
        framelink_tags = metadata_map.fill_frame_info_cbf(len(imgfiles), imgfiles)

    tags = {**byte_tags, **metadata_tags, **framelink_tags}
    metadata_map.write_imgcif(tags)

if __name__ == '__main__':

    main()

