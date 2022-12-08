
import click
import os
import CifFile
from imgCIF_Creator.command_line_interfaces import parser
from imgCIF_Creator.information_extractors import user_input
from imgCIF_Creator.information_extractors import cbf_smv
from imgCIF_Creator.information_extractors import hdf5_NxMx


@click.command()
@click.option(
    "--gui",
    "-g",
    default=False,
    type=bool,
    is_flag=True,
    help="Start the imgCIF_creator with a graphical user interface.",
)

def main(gui):

    if gui:
        graphical_user_interface()
    else:
        command_line_interface()



def command_line_interface():

    print('\n--------------------------- imgCIF Creator ---------------------------\n' )
    print("""This is an interactive command line interface to collect the necessary \
information to create an imgCIF file out of HDF5, full CBF and some common subset \
of miniCBF. You can skip parameters that are not required with an empty input,\
if you however provide an input that will be checked against the required format.""")

    # choose_filetype()

    # cmd_parser = parser.CommandLineParser()

    # filename
    # cmd_sparser.validated_user_input('Filename')


    cif_file = CifFile.CifFile()
    cif_block = CifFile.CifBlock()
    cif_file['imgCIF'] = cif_block



    ## mini CBF
    path = '/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/cbf_cyclohexane_crystal2/CBF_crystal_2'
    cbf_smv.extract_cbf_scan(
        cif_block,
        path,
        axis=['omega', 'my_new_name'],
        new_url='my_new_name.tgz',
        include_archive_directory=True
        )



    ## Nexus NxMx
    # path = "/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/_b4_1_master.h5"
    # hdf5_NxMx.extract_hdf5_NxMx_scan(path, cif_block)



    # needed for hdf5 from user
    # _database.dataset_doi        '10.5281/zenodo.5886687'
    # _audit.block_code            Diamond_I04
    # _diffrn_source.beamline      I04
    # _diffrn_source.facility      Diamond
    # _diffrn_radiation.type       '{td['RADN_TYPE']}'
    # _array_intensities.overload  65535

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
    user_input.convert_user_input_to_imgcif(parsed_input, cif_block)


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