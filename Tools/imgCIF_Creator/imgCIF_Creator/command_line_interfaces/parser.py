
import os
import yaml
import re
import numpy as np
from inspect import getsourcefile


class CommandLineParser():

    def __init__(self) -> None:

        self.parsed = {}
        self._resources_path = os.path.abspath(getsourcefile(lambda: 0)).replace(
            'command_line_interfaces/parser.py', 'resources/')

        # opening a file
        with open(f"{self._resources_path}user_input.yaml", 'r') as stream:
            try:
                self.input_options = yaml.safe_load(stream)
            except yaml.YAMLError as e:
                print(e)


        self.validation_regex = {
            'doi': None, # not checking doi format yet
            'layout': r'Beamline|Laboratory',
            'facility': None,
            'beamline': None,
            'rad_type': None,
            "model" : None,
            "location" : None,
            "manufacturer" : None,
            # 'Is updated?': r'yes|no\Z', # no effect?
            # 'Start date': r'\d{4}-\d{2}-\d{2}\Z', # does not prevent the user from putting
            # in unrealistic values, no effect?
            # 'Finish date': r'\d{4}-\d{2}-\d{2}\Z', # no effect?
            'principal_orientation': r'\d{1,3}\Z',
            'goniometer_axes': r'(?:.+,\s*(a|c))+\Z', # matches n patterns xxxxx, a|c
            'rotation_axis': None, # no effect?
            "two_theta_axis" : r'(a|c|anticlockwise|clockwise)\Z', # no effect?
            # "detector_repositioning" : None, # no effect?
            'detector_repositioning': None, # no effect? almost the same as the above
            "detector_axes" : None, #r'(?:\S*,\s*)+(\S*)\Z', # this is not exact
            # r"(?:\d+,\s*)+(?:\d+)\s*\Z",
            "chi_axis" : r'.*(?:\s+\d{1,3})\Z',
            'kappa_axis': r'.*(?:\s+\d{1,3}){1,2}\Z', # capture something like kappa, 50, 50
            'image_orientation':r'(top left|top right|bottom left|bottom right)\Z',
            'fast_direction': r'(horizontal|vertical)\Z',
            'pixel_size': r'(\d+\.?\d*\Z)|(\d+\.?\d*,\s*\d+\.?\d*\Z)', # matches either
            # one pixel dimension or both separated by a comma
            # 'pixel_number': r'(\d+,{0,1}\s*\d+\Z)', # matches two numbers separated by
            'pixel_number': r'(\d+(,\s*)\d+\Z)',
            # a comma or a whitespace
            # 'xds.inp': None, # no effect
            'doi': None, # no effect?
            'comments': None, # no effect?
            'filename' : r'.*((\.h5)\Z|(\.cbf)\Z|(\.smv)\Z)',
            'goniometer_rot_direction' : r'(clockwise|anti-clockwise|c|a)\Z',
            'frame_numbers' : r'^\d+(,\s*\d+)*$'
        }

        self.required_information = {
            'doi': True,
            'layout': True,
            'facility': True,
            'beamline': True,
            "model" : True,
            "location" : True,
            "manufacturer" : True,
            # 'Is updated?': r'yes|no\Z', # no effect?
            # 'Start date': r'\d{4}-\d{2}-\d{2}\Z', # does not prevent the user from putting
            # in unrealistic values, no effect?
            # 'Finish date': r'\d{4}-\d{2}-\d{2}\Z', # no effect?
            'principal_orientation': True,
            'goniometer_axes': True, # matches n patterns xxxxx, a|c
            'rotation_axis': False, # no effect?
            "two_theta_axis" : False, # no effect?
            # "detector_repositioning" : None, # no effect?
            'detector_repositioning': False, # no effect? almost the same as the above
            'detector_axes' : True,
            "chi_axis" : False,
            'kappa_axis': False, # capture something like kappa, 50, 50
            'image_orientation': True,
            'fast_direction': True,
            'pixel_size': True, # matches either
            # one pixel dimension or both separated by a comma
            # 'xds.inp': None, # no effect
            'pixel_number' : True,
            'doi': True, # no effect?
            'comments': False, # no effect?
            'filename' : True,
            'goniometer_rot_direction' : True,
            'frame_numbers': True
        }


    def request_input(self, label):

        while self.parsed.get(label) is None:
            required = 'required' if self.required_information.get(label) \
                else 'not required'
            print(f"\n{self.input_options[label]['label']} ({required}):")
            choices = ''
            if self.input_options[label].get('options') is not None:
                choices = ' (choices: ' +\
                    ', '.join(self.input_options[label]['options']) + ')'
            print(f"{self.input_options[label]['description']}".strip('\n') + choices)
            self.parsed[label] = self.validated_user_input(label)

        return self.parsed[label]


    def validated_user_input(self, input_label):

        user_input = input(' >> ')
        # validation pattern
        if self.validation_regex.get(input_label) is not None:
            pattern = re.compile(
                self.validation_regex[input_label])
            parsed_input = pattern.match(user_input)
            if parsed_input:
                parsed_input = parsed_input.group(0)
        else:
            parsed_input = user_input

        # special cases
        # TODO still needed?
        if input_label == 'rotation_axis' and user_input != '':
            axes, _ = user_input.parse_axis_string(self.parsed['goniometer_axes'])
            parsed_input = parsed_input if parsed_input in axes else None

        # if input_label == 'detector_axes':
        #     parsed_input = parsed_input.split(',')

        # required parameter, but either regex failed or no input was made if no regex
        # is defined
        if self.required_information.get(input_label) and parsed_input in [None, '']:
            print(' ==> Could not interpret your input correctly! This input is required, \
please try again.')
            parsed_input = None
        # not required with defined regex, but no user input
        elif not self.required_information.get(input_label) and user_input in ['']:
            parsed_input = ''
        # not required with defined regex, but user input
        elif not self.required_information.get(input_label) and parsed_input is None:
            print(' ==> Could not interpret your input correctly! Please try again.')

        if parsed_input is not None:
            print(f" ==> Your input was: {parsed_input}")

        return parsed_input


    def parse_axis_string(self, x):
        """
        Parse a string of form axis,sense,axis,sense...
        """
        s = [sub.strip() for sub in x.split(',')] #strip.(split(x, ","))
        # print('mysub', s)
        if len(s) % 2 != 0:
            raise Exception("Axis string is incorrect: %x")

        axes = []
        senses = []
        for i in range(0, len(s), 2): #1:2:length(s)
            axes.append(s[i])
            a = s[i+1].lower()
            # print('my a', a)
            if a == "clockwise":
                a = "c"
            elif a == "anticlockwise":
                a = "a"
            elif a not in ["a", "c"]:
                raise Exception(f"Unrecognised rotation sense {a}")

            senses.append(a)

        return axes, senses


    def make_goniometer_axes(self, goniometer_axes, kappa_axis, chi_axis):

        """
        make_gonio_axes(raw_info)

        Given a list of gonio axes, create their representation in imgCIF. The list
        of gonio axes goes in order from top to bottom, meaning that the first
        "depends on" the second and so forth. We assume a two theta axis.
        The items we have to fill in are:
        1. type -> rotation
        2. depends_on -> next in list
        3. equipment -> goniometer
        4. vector -> almost always 1 0 0 (rotation about principal axis)
        5. offset -> always [0 0 0] but needed for loop integrity

        Note that our questions assume looking from above whereas imgCIF is looking from
        below, so the sense of rotation is reversed.

        """
        print('gonx', goniometer_axes)

        axes, senses = goniometer_axes
        # n = len(axes)
        axis_type = ["rotation" for _ in axes]
        equip = ["goniometer" for _ in axes]

        if kappa_axis is None:
            kappa_axis = "- - -"
        kappa_axis = re.split(r',| ', kappa_axis)[0]

        if chi_axis is None:
            chi_axis = "- - -"
        chi_axis = re.split(r',| ', chi_axis)[0]

        # Construct axis dependency chain
        depends_on = []
        depends_on += axes[1:]
        depends_on.append('none')

        # Remember the principal axis direction
        principal = senses[-1]

        # Create direction vectors
        vector = []
        offset = []
        for (a, d) in zip(axes, senses):
            print(a, d)
            rotfac = 1 if d == principal else -1
            if a.lower() == kappa_axis.lower():
                kv = self.get_kappa_info(goniometer_axes, kappa_axis)
                kv[1] *= rotfac
                vector.append(kv)
            elif a.lower() == chi_axis.lower():
                vector.append(self.get_chi_info(goniometer_axes, chi_axis))
            else:
                vector.append([i * rotfac for i in [1, 0, 0]])
            print('vec', vector)
            offset.append([0, 0, 0])

        axes_dict = {
            'axes' : goniometer_axes[0],
            'axis_type' : axis_type,
            'equip' : equip,
            'depends_on' : depends_on,
            'vector' : vector,
            'offset' : offset,
        }
        return axes_dict


    def make_detector_axes(self,
        goniometer_axes, principal_orientation, image_orientation, two_theta_axis):
        """
            make_detector_axes(raw_info)

        Add information concerning the detector axes. We define our own axis names,
        with the detector distance being inserted when the data file is read. We
        choose det_x to be in the horizontal direction, and det_y to be vertical.
        We need to add:
        1. type -> translation
        2. depends_on -> x,y depend on translation
        3. equipment -> detector
        4. vector -> worked out from user-provided info
        5. offset -> beam centre
        """

        # Insert assumed axis names and orientations
        axis_id = ["two_theta", "trans", "detx", "dety"]
        axis_type = ["translation" for _ in axis_id]
        axis_type[1] = "rotation"
        equip = ["detector" for _ in axis_id]
        depends_on = ['none', "two_theta", "trans", "detx"]

        # Read necessary information
        print('prinsens', goniometer_axes)
        principal_sense = goniometer_axes[1][-1]
        principal = principal_orientation
        corner = image_orientation

        # Adjust two theta direction #TODO ensure twotheta exists
        rotsense = 1 if two_theta_axis == principal_sense else -1
        vector = [[rotsense, 0, 0]]

        # Detector translation always opposite to beam direction
        vector.append([0, 0, -1])

        # Work out det_x and det_y
        x_d, y_d = self.determine_detx_dety(principal, principal_sense, corner)
        vector.append(x_d)
        vector.append(y_d)

        # Beam centre is unknown for now
        # offset = [[0, 0, 0], [0, 0, 0], [np.nan, np.nan, 0], [0, 0, 0]]
        # TODO this should probably not be none
        offset = [[0, 0, 0], [0, 0, 0], ['none', 'none', 0], [0, 0, 0]]

        axes_dict = {
            'axes' : axis_id,
            'axis_type' : axis_type,
            'equip' : equip,
            'depends_on' : depends_on,
            'vector' : vector,
            'offset' : offset,
        }
        return axes_dict


    def get_kappa_info(self, goinometer_axes, kappa_axis):

        # extract kappa information
        axes, senses = goinometer_axes
        # TODO keepempty=false is ok with python?
        kinfo = kappa_axis.split([' ',','])

        if len(kinfo) == 2:
            kinfo.append("0")

        axname = kinfo[0]
        if not axname in axes:
            # throw(error("Kappa axis $axname not listed as a goniometer axis"))
            # TODO reaise correct here
            raise Exception("Kappa axis $axname not listed as a goniometer axis")

        kapang = float(kinfo[1])
        kappoise = float(kinfo[2])

        # Now calculate direction of axis in X-Z plane assuming
        # rotation of kappa is same as principal axis. There is no
        # Y component as home position is always under
        # beam path.
        up_comp = np.cos(kapang)
        across_comp = np.sin(kapang)
        if kappoise == 0:
            #is under incident beam collimator in X-Z plane
            z_comp = -1 * across_comp
        elif kappoise == 180:
            z_comp = across_comp

        return [up_comp, 0.0, z_comp]


    def get_chi_info(self, goniometer_axes, chi_axis):

        # Extract chi information
        axes, senses = goniometer_axes
        # TODO keepempty?
        cinfo = chi_axis.split([' ',',']) #,keepempty=false)
        axname = cinfo[0].lower() #lowercase( first( cinfo ) )
        if not axname in axes.lower():
            raise Exception("Chi axis name $axname not listed as a goniometer axis")
        r = -1 * float(cinfo[2])

        # Now turn this into an axis direction
        # At the provided rotation, the chi axis is parallel to z. It is rotated
        # by -chi_rot about omega to bring it to the start position. The sense
        # of omega rotation is always about 1 0 0 by definition
        # chi_sense = indexin([axname], lowercase.(axes))[]
        #TODO correct?
        axes_lower = axes.lower()
        chi_sense = [axes_lower.index(val) for val in [axname]]
        chi_sense = senses[chi_sense]

        # TODO arrays?
        chi_beam_dir = [0, 0, 1] if chi_sense == "a" else [0, 0, -1]
        chi_rot = [
            [1, 0, 0],
            [0, np.cos(r), -np.sin(r)],
            [0, np.sin(r), np.cos(r)]
        ]

        return np.array(chi_rot) * np.array(chi_beam_dir)


    def determine_detx_dety(self, principal_angle, principal_sense, corner):

        # Determine direction of detx (horizontal) and dety (vertical) in
        # imgCIF coordinates.
        # Start with basic value and then flip as necessary

        x_direction = [-1, 0, 0]        # spindle rotates anticlockwise at 0, top_left origin
        y_direction = [0, 1, 0]       #

        #TODO go for np arrays?
        if corner == "top right":
            x_direction = [i * -1 for i in x_direction]

        elif corner == "bottom right":
            x_direction = [i * -1 for i in x_direction]
            y_direction = [i * -1 for i in y_direction]

        elif corner == "bottom left":
            y_direction = [i * -1 for i in y_direction]

        # The direction of the principal axis flips by 180 if the sense
        # changes
        pa = int(principal_angle) if principal_sense == "a" else int(principal_angle) + 180
        if pa >= 360:
            pa = pa - 360

        if pa == 90:
            temp = x_direction
            x_direction = y_direction
            y_direction = [i * -1 for i in temp] #-1*temp
        elif pa == 180:
            x_direction = [i * -1 for i in x_direction]
            y_direction = [i * -1 for i in y_direction]
        elif pa == 270:
            temp = x_direction
            x_direction = [i * -1 for i in y_direction] #-1*y_direction
            y_direction = temp

        return x_direction,y_direction

