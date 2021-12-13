from argparse import ArgumentParser
import h5py as h5
import sys

DLS_I04_TEMPLATE = """\
_audit.block_id	Diamond_I04
_diffrn_source.beamline	I04
_diffrn_source.facility	Diamond

    _array_structure_byte_order         ?
    _array_structure_compression_type   ?
    _array_structure.encoding_type      ?

    loop_
      _diffraction_radiation.type
      _diffraction_radiation_wavelength.value
         %(RADN_TYPE)s          %(RADN_WAVELEN)s

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
         two_theta  rotation     detector    .          ?   0  0  0  0  0
         trans      translation  detector    two_theta  0   0  1  0  0  %(DET_ZDIST)s
         detx       translation  detector    trans      0   1  0  %(DET_BMCENT_X)s  0  0
         dety       translation  detector    trans      -1  0  0  0  %(DET_BMCENT_Y)s  0

    loop_
      _array_structure_list_axis.axis_id
      _array_structure_list_axis.axis_set_id
      _array_structure_list_axis.start
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
         1                        2

    loop_
      _diffrn_detector_axis.axis_id
      _diffrn_detector_axis.detector_id
         detx                     1
         dety                     1

    _diffrn_scan_axis.axis_id                %(SCAN_AXIS)s
    _diffrn_scan_axis.angle_start            %(SCAN_START)s
    _diffrn_scan_axis.displacement_start     0  0
    _diffrn_scan_axis.angle_increment        %(SCAN_INCR)s
    _diffrn_scan_frame_axis.angle_increment  %(SCAN_INCR)s
"""

class DLS_I04_MAP:

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
           'axisid_tag': 'SCAN_AXIS',
           'axisid_name': '',
           'axstart_tag': 'SCAN_START',
           'axstart_range': '',
           'axincr_tag': 'SCAN_INCR',
           'axincr_range': ''
        }
    

    def extract_from_hdf5(self, fn):

        tags_dict = {}
        scan_axes = {'phi': False, 'chi': False, 'omega': False }

        with h5.File(fn) as f:
            for k in ['phi', 'chi', 'omega']:

                print('Collect INFO for', k)
                h5item = f[self.items[f'{k}_info']['path']]
                if len(h5item[()]) > 1:
                    print(' ... identified as scan axis')
                    self.items['scan_info']['axisid_name'] = k
                    self.items['scan_info']['axstart_range'] = \
                      f'{h5item[()][0]} {h5item[()][-1]}'
                    incr_range = f[f'entry/sample/transformations/{k}_increment_set'][()]
                    self.items['scan_info']['axincr_range'] = \
                      f'{incr_range[0]} {incr_range[-1]}'


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
            for tag, item in [(sdict['axisid_tag'], sdict['axisid_name']),
                              (sdict['axstart_tag'], sdict['axstart_range']),
                              (sdict['axincr_tag'], sdict['axincr_range'])
                             ]:
                tags_dict[tag] = item

        return tags_dict

    def write_imgcif(self, tags):

        with open('cbf_metadata.cif', 'w') as f:
            f.write(self.template % tags)


def main():

    print('imgCIF mapping tool')
    ap = ArgumentParser(prog='imgcif_mapper')
    ap.add_argument('infile', type=str, help='path/name of input HDF5 file with metadata')
    args = ap.parse_args(sys.argv[1:])
    print('metadata source is input file:', args.infile)
    metadata_map = DLS_I04_MAP()
    tags = metadata_map.extract_from_hdf5(args.infile)
    print(tags)
    metadata_map.write_imgcif(tags)

if __name__ == '__main__':
    main()
