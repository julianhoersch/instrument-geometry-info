
import sys
import click
import os
import re
import CifFile
from imgCIF_Creator.command_line_interfaces import parser
from imgCIF_Creator.information_extractors import user_input
from imgCIF_Creator.information_extractors import cbf_smv
from imgCIF_Creator.information_extractors import full_cbf
from imgCIF_Creator.information_extractors import hdf5_NxMx
from imgCIF_Creator.output_assembler import imgCIF_assembler

# how generic can it be? catch errors if entries are not found or provide different
# ways for different facilities etc?


# def validate_filename(context, param, filename):
def validate_filename(filename):
    """_summary_

    A callback is a function that is invoked with two parameters: the current
    Context and the value. The context provides some useful features such as
    quitting the application and gives access to other already processed parameters.

    Args:
        ctx (_type_): _description_
        param (_type_): _description_
        value (_type_): _description_
    """

    if os.path.isdir(filename):
        for _, _, files in os.walk(filename + os.sep):
            pattern = re.compile(r'.*((?P<cbf>\.cbf)\Z|(?P<smv>\.smv)\Z)')
            matches = [bool(pattern.match(file)) for file in files]
            occurences = len([match for match in matches if match == True])
            if occurences > 0:
                print('Found a folder with .cbf or .smv files.')
                filetype = 'cbf'
            else:
                print('Could not find .cbf or .smv files in directory.')
                sys.exit()
    elif os.path.isfile(filename):
        regex = r'.*((?P<h5>\.h5)\Z)'
        match = re.match(regex, filename)
        filetype = 'h5'
        if not match:
            print('Only h5 (NxMx) files are supported! If you want to convert \
.cbf and .smv files please provide a directory. Exiting.')
            sys.exit()

    # regex = r'.*((?P<h5>\.h5)\Z|(?P<cbf>\.cbf)\Z|(?P<smv>\.smv)\Z)'
    # match = False
    # while not match:
#     print('\nPlease enter the filename of the file you want to convert into \
# imgCIF. Accepted files are .cbf, .smv and h5 (NxMx).')
#     filename = input(' >> ')
    # match = re.match(regex, filename)
    # if not match:
    #     print('Only .cbf, .smv and h5 (NxMx) files are supported! Exiting.')
    #     sys.exit()
    # else:
    # for group in ['smv', 'cbf', 'h5']:
    #     if match.group(group) is not None:
    #         print(f'\nIdentified input file as .{group} file.')
    #         return (match.group(), group)

    #     print(match.group('smv'))


    return filename, filetype



@click.command()
@click.option(
    "--gui",
    "-g",
    default=False,
    type=bool,
    is_flag=True,
    help="Start the imgCIF_creator with a graphical user interface.",
)

@click.argument(
    "filename",
    type=str,
    # help="The filename of the file that should be converted",
    # callback=validate_filename,
)

def main(filename, gui):
    """A tool to do things with a FILENAME

    Args:
        gui (_type_): _description_
    """


    print('\n--------------------------- imgCIF Creator ---------------------------\n' )
    print("""This is an interactive command line interface to collect the necessary \
information to create an imgCIF file out of HDF5, full CBF and some common subset \
of miniCBF. You can skip parameters that are not required with an empty input,\
if you however provide an input that will be checked against the required format.""")

    filename, filetype = validate_filename(filename)


    if gui:
        graphical_user_interface()
    else:
        command_line_interface(filename, filetype)



def command_line_interface(filename, filetype):


    print(f'working on {filename}')

    # return

    # cmd_parser = parser.CommandLineParser()

    # filename
    # cmd_sparser.validated_user_input('Filename')


    cif_file = CifFile.CifFile()
    cif_block = CifFile.CifBlock()
    cif_file['imgCIF'] = cif_block


    # filetype = 'cbf'


    ## mini CBF
    if filetype in ['smv', 'cbf']:
        # filename = '/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/cbf_cyclohexane_crystal2/CBF_crystal_2'

        # filename = '/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/010_Ni_dppe_Cl_2_150K01'
        # cif_file = full_cbf.extract_full_cbf_header_information(filename)


        imgCIF_assembler.create_imgCIF(filename)

        sys.exit()

        parsed_input = {
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
        "Array dimension" : ['none', 'none']}

        # uncomment for input parser
        parsed_input = parser.CommandLineParser().retrieve_user_information()
        parsed_input["Array dimension"] = [100, 100]

        # convert parsed input
        user_input.convert_user_input_to_imgcif(parsed_input, cif_block)


        scan_info, all_frames = cbf_smv.extract_cbf_scan(
            cif_block,
            filename,
            axis=['omega', 'my_new_name'],
            # file_stem='010_Ni_dppe_Cl_2_150K',
            )

        include_archive_directory = True
        if include_archive_directory:
            prepend_dir = os.path.split(filename)[-1]
        else:
            prepend_dir = ""

        imgCIF_assembler.add_scan_info_to_block(
            scan_info, all_frames, cif_block, 'my_new_name.tgz',
            prepend_dir, filename, 'CBF')

    elif filetype == 'h5':

        # Nexus NxMx
        # path = "/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/_b4_1_master.h5"
        hdf5_NxMx.extract_hdf5_NxMx_scan(filename, cif_block)



    # needed for hdf5 from user
    # _database.dataset_doi        '10.5281/zenodo.5886687' # always zenodo?
    # _audit.block_code            Diamond_I04
    # _diffrn_source.beamline      I04
    # _diffrn_source.facility      Diamond
    # _diffrn_radiation.type       '{td['RADN_TYPE']}'
    # _array_intensities.overload  65535 #TODO do I need this?

    parsed_input = {
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

    parsed_input = {
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
        "Array dimension" : ['none', 'none']}

    # uncomment for input parser
    # parsed_input = parser.CommandLineParser().retrieve_user_information()

    # convert parsed input
    # user_input.convert_user_input_to_imgcif(parsed_input, cif_block)


    # change order
    # for idx, item in enumerate(OUTPUT_ORDER):
    #     try:
    #         cif_file['imgCIF'].ChangeItemOrder(item, idx)
    #     except ValueError:
    #         continue

    # print('cblovck', cif_block)
    # print(cif_file.WriteOut(wraplength=1000))#, maxoutlength=2048))
    # with open("imgCIF_input.cif", 'w') as file:
    #     file.write(cif_file.WriteOut())

    fname = os.getcwd() + os.sep + "imgCIF_test.cif"
    with open(fname, 'w') as file:
        file.write(cif_file.WriteOut())
        print(f'File saved to: {fname}')



#     _item_description.description
# ;             The value of _diffrn_detector.number_of_axes  gives the
#             number of axes of the positioner for the detector identified
#             by _diffrn_detector.id.


# def choose_filetype():

#     print("Please give the filename of the file you want to convert:")
#     ino

#     return filename, filetype



def graphical_user_interface():

    print('GUI not implemented yet!')



