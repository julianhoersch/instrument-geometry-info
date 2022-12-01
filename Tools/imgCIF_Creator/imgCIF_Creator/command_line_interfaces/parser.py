

# Program to convert yaml file to dictionary
import os
import yaml
import re
from inspect import getsourcefile
from imgCIF_Creator.information_extractors import user_input_extractor


class CommandLineParser():

    def __init__(self) -> None:

        self.validation_regex = {
            'layout': r'beamline|laboratory',
            'Facility name': None,
            'Beamline name': None,
            "Model" : None,
            "Location" : None,
            "Name of manufacturer" : None,
            'Is updated?': r'yes|no\Z', # no effect?
            'Start date': r'\d{4}-\d{2}-\d{2}\Z', # does not prevent the user from putting
            # in unrealistic values, no effect?
            'Finish date': r'\d{4}-\d{2}-\d{2}\Z', # no effect?
            'Principal axis orientation': r'\d{1,3}\Z',
            'Goniometer axes': r'(?:.+,\s*(a|c))+\Z', # matches n patterns xxxxx, a|c
            'Rotation axis name': None, # no effect?
            "Two theta axis" : r'(a|c|anticlockwise|clockwise)\Z', # no effect?
            "Other detector axes" : None, # no effect?
            'Detector axes': None, # no effect? almost the same as the above
            "chi" : r'.*(?:\s+\d{1,3})\Z',
            'kappa': r'.*(?:\s+\d{1,3}){1,2}\Z', # capture something like kappa, 50, 50
            'Image orientation':r'(top left|top right|bottom left|bottom right)\Z',
            'Fast direction': r'(horizontal|vertical)\Z',
            'Pixel size': r'(\d+\.?\d*\Z)|(\d+\.?\d*,\s*\d+\.?\d*\Z)', # matches either
            # one pixel dimension or both separated by a comma
            'xds.inp': None, # no effect
            'Data DOI': None, # no effect?
            'Comments': None, # no effect?
        }

        self.required_information = {
            'layout': True,
            'Facility name': True,
            'Beamline name': True,
            "Model" : True,
            "Location" : True,
            "Name of manufacturer" : True,
            'Is updated?': False, # no effect
            'Start date': False, # no effect
            'Finish date': False, # no effect
            'Principal axis orientation': True,
            'Goniometer axes': True,
            'Rotation axis name': False, # no effect
            "Two theta axis" : False, # no effect
            "Other detector axes" : False, # no effect
            'Detector axes': False, # no effect
            "chi" : False,
            'kappa': False,
            'Image orientation': True,
            'Fast direction': True,
            'Pixel size': True,
            'xds.inp': False, # no effect
            'Data DOI': False, # no effect
            'Comments': False, # no effect
        }

        self.parsed = {}
        self._resources_path = os.path.abspath(getsourcefile(lambda: 0)).replace(
            'command_line_interfaces/parser.py', 'resources/')


    def retrieve_user_information(self):

        print('\n--------------------------- imgCIF Creator ---------------------------\n' )
        print("""This is an interactive command line interface to collect the necessary \
information to create an imgCIF file out of HDF5, full CBF and some common subset \
of miniCBF. You can skip parameters that are not required with an empty input,\
if you however provide an input that will be checked against the required format.""")

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
        if self.validation_regex[input_label] is not None:
            pattern = re.compile(
                self.validation_regex[input_label])
            parsed_input = pattern.match(user_input)
            if parsed_input:
                parsed_input = parsed_input.group(0)
        else:
            parsed_input = user_input

        # special cases
        if input_label == 'Rotation axis name' and user_input != '':
            axes, _ = user_input_extractor.parse_axis_string(self.parsed['Goniometer axes'])
            parsed_input = parsed_input if parsed_input in axes else None

        # required parameter, but either regex failed or no input was made if no regex
        # is defined
        if self.required_information[input_label] and parsed_input in [None, '']:
            print(' ==> Could not interpret your input correctly! This input is required, \
please try again.')
            parsed_input = None
        # not required with defined regex, but no user input
        elif not self.required_information[input_label] and user_input in ['']:
            parsed_input = ''
        # not required with defined regex, but user input
        elif not self.required_information[input_label] and parsed_input is None:
            print(' ==> Could not interpret your input correctly! Please try again.')

        if parsed_input is not None:
            print(f" ==> Your input was: {parsed_input}")

        return parsed_input


# if __name__ == '__main__':
    # retrieve_user_information()
