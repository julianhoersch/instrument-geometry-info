

# Program to convert yaml file to dictionary
import os
import yaml
import re
from inspect import getsourcefile
from imgCIF_Creator.information_extractors import user_input


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

        # print('input opt', self.input_options)


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



        # self.validation_regex = {
        #     'layout': r'Beamline|Laboratory',
        #     'Facility name': None,
        #     'Beamline name': None,
        #     "Model" : None,
        #     "Location" : None,
        #     "Name of manufacturer" : None,
        #     'Is updated?': r'yes|no\Z', # no effect?
        #     'Start date': r'\d{4}-\d{2}-\d{2}\Z', # does not prevent the user from putting
        #     # in unrealistic values, no effect?
        #     'Finish date': r'\d{4}-\d{2}-\d{2}\Z', # no effect?
        #     'Principal axis orientation': r'\d{1,3}\Z',
        #     'Goniometer axes': r'(?:.+,\s*(a|c))+\Z', # matches n patterns xxxxx, a|c
        #     'Rotation axis name': None, # no effect?
        #     "Two theta axis" : r'(a|c|anticlockwise|clockwise)\Z', # no effect?
        #     "Other detector axes" : None, # no effect?
        #     'Detector axes': None, # no effect? almost the same as the above
        #     "chi" : r'.*(?:\s+\d{1,3})\Z',
        #     'kappa': r'.*(?:\s+\d{1,3}){1,2}\Z', # capture something like kappa, 50, 50
        #     'Image orientation':r'(top left|top right|bottom left|bottom right)\Z',
        #     'Fast direction': r'(horizontal|vertical)\Z',
        #     'Pixel size': r'(\d+\.?\d*\Z)|(\d+\.?\d*,\s*\d+\.?\d*\Z)', # matches either
        #     # one pixel dimension or both separated by a comma
        #     'xds.inp': None, # no effect
        #     'Data DOI': None, # no effect?
        #     'Comments': None, # no effect?
        #     'Filename' : r'.*((\.h5)\Z|(\.cbf)\Z|(\.smv)\Z)',
        # }

        # self.required_information = {
        #     'layout': True,
        #     'Facility name': True,
        #     'Beamline name': True,
        #     "Model" : True,
        #     "Location" : True,
        #     "Name of manufacturer" : True,
        #     'Is updated?': False, # no effect
        #     'Start date': False, # no effect
        #     'Finish date': False, # no effect
        #     'Principal axis orientation': True,
        #     'Goniometer axes': True,
        #     'Rotation axis name': False, # no effect
        #     "Two theta axis" : False, # no effect
        #     "Other detector axes" : False, # no effect
        #     'Detector axes': False, # no effect
        #     "chi" : False,
        #     'kappa': False,
        #     'Image orientation': True,
        #     'Fast direction': True,
        #     'Pixel size': True,
        #     'xds.inp': False, # no effect
        #     'Data DOI': False, # no effect
        #     'Comments': False, # no effect
        #     'Filename' : True,
        # }


    def request_input(self, label):


        while self.parsed.get(label) is None:
            required = 'required' if self.required_information.get(label) \
                else 'not required'
            print(f"\n{self.input_options[label]['label']} ({required}):")
            choices = ''
            if self.input_options[label].get('options') is not None:
                choices = ' (choices: ' + ', '.join(self.input_options[label]['options']) + ')'
            print(f"{self.input_options[label]['description']}".strip('\n') + choices)
            # print(f"{content['attributes']['description']}".strip('\n'))
            self.parsed[label] = self.validated_user_input(label)

        return self.parsed[label]



    def get_beamline_lab_choice(self):

        while self.parsed.get('layout') not in ['beamline', 'laboratory']:
            print('\nDo you want to create an imgCIF for beamline or laboratory data (choices: beamline or laboratory)?')
            self.parsed['layout'] = self.validated_user_input('layout')

        # choice = {'beamline' : 'beamline-layout.yaml',
        #         'laboratory' : 'lab-layout.yaml',
        #         }


        # return choice


    def retrieve_user_information(self):

#         print('\n--------------------------- imgCIF Creator ---------------------------\n' )
#         print("""This is an interactive command line interface to collect the necessary \
# information to create an imgCIF file out of HDF5, full CBF and some common subset \
# of miniCBF. You can skip parameters that are not required with an empty input,\
# if you however provide an input that will be checked against the required format.""")

        while self.parsed.get('layout') not in ['beamline', 'laboratory']:
            print('\nDo you want to create an imgCIF for beamline or laboratory data (choices: beamline or laboratory)?')
            self.parsed['layout'] = self.validated_user_input('layout')

        choice = {'beamline' : 'beamline-layout.yaml',
                'laboratory' : 'lab-layout.yaml',
                }

        # opening a file
        with open(f"{self._resources_path}{choice[self.parsed['layout']]}", 'r') as stream:
            try:
                user_requested_input = yaml.safe_load(stream)
            except yaml.YAMLError as e:
                print(e)

        for content in user_requested_input['body']:
            if content['attributes'].get('value') is not None:
                print('\n', content['attributes']['value'])
            if content['attributes'].get('description') is not None:

                while self.parsed.get(content['attributes']['label']) is None:
                    required = 'required' if self.required_information[content['attributes']['label']] \
                        else 'not required'
                    print(f"\n{content['attributes']['label']} ({required}):")
                    choices = ''
                    if content['attributes'].get('options') is not None:
                        choices = ' (choices: ' + ', '.join(content['attributes']['options']) + ')'
                    print(f"{content['attributes']['description']}".strip('\n') + choices)
                    # print(f"{content['attributes']['description']}".strip('\n'))
                    self.parsed[content['attributes']['label']] = \
                        self.validated_user_input(content['attributes']['label'])

        # print(self.parsed)
        return self.parsed


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


# if __name__ == '__main__':
    # retrieve_user_information()
