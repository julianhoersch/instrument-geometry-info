
import click
import CifFile
import imgCIF_Creator.command_line_interfaces.parser as cmd_parser
import imgCIF_Creator.information_extractors.user_input_extractor as user_input_extractor
import imgCIF_Creator.information_extractors.cbf_smv_extractor as cbf_smv_extractor
import imgCIF_Creator.information_extractors.hdf5_NxMx_extractor as hdf5_NxMx


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

    file_format = 'HDF5'

    cif_file = CifFile.CifFile()
    cif_block = CifFile.CifBlock()
    cif_file['imgCIF'] = cif_block

    # path = '/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/cbf_cyclohexane_crystal2/CBF_crystal_2'
    # cbf_smv_extractor.extract_cbf_scan(
    #     cif_block,
    #     path,
    #     file_format,
    #     axis=['omega', 'my_new_name'],
    #     new_url='my_new_name.tgz',
    #     include_archive_directory=True
    #     )

    path = "/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/_b4_1_master.h5"
    hdf5_NxMx.extract_hdf5_NxMx_scan(path, cif_block, file_format)

    # return


    if gui == True:
        print('GUI not implemented yet!')
    else:

        # parsed_input = cmd_parser.CommandLineParser().retrieve_user_information()

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
            "Pixel size" : "0.172"}

        # user_input_extractor.convert_user_input_to_imgcif(parsed_input)

        # describe_axes(user_input, cif_block)
        # print('citems', cif_block.items())

        # for idx, item in enumerate(OUTPUT_ORDER):
        #     try:
        #         cif_file['imgCIF'].ChangeItemOrder(item, idx)
        #     except ValueError:
        #         continue

        # print('cblovck', cif_block)
        # print(cif_file.WriteOut(wraplength=1000))#, maxoutlength=2048))
        with open("mycif.cif", 'w') as file:
            file.write(cif_file.WriteOut())


