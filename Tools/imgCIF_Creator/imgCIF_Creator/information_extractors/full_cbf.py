import CifFile
import cbf


from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull

@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


class text_container:

    def __init__(self, text) -> None:

        self.text = text

    def remove_char_from_tail(self):
        self.text = self.text[:-1]

    def read(self):
        return self.text


def extract_full_cbf_header_information(filename):

    # filename = '/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/010_Ni_dppe_Cl_2_150K01/010_Ni_dppe_Cl_2_150K01_00001.cbf'
    # print('fname', fname)
    with open(filename, 'rb') as file:
        # lines = file.readlines()

        # text = file.read()#.split("-BINARY-FORMAT-SECTION-")
        file_content = file.read()
        header_end_mark = b'\x0C\x1A\x04\xD5'
        header_end_index = file_content.find(header_end_mark)
        text = file_content[:header_end_index].decode('utf-8')
        text = text.split("--CIF-BINARY-FORMAT-SECTION--")[0]
        # print('hakunamatata', text)
        # print('hkanmstuta')

        # lines = []

        # for line in file:
        #     # print('my lane', line)
        #     if b"--CIF-BINARY-FORMAT-SECTION--" in line:
        #         break
        #     else:
        #         lines.append(line.decode().strip('\n').lower())

    # cf = CifFile.ReadStar('/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/010_Ni_dppe_Cl_2_150K01/010_Ni_dppe_Cl_2_150K01_00001.cbf',
    #     CBF=True)

    container = text_container(text)
    # filename = '/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/010_Ni_dppe_Cl_2_150K01/010_Ni_dppe_Cl_2_150K01_00001.cbf'
    # filename = '/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/010_Ni_dppe_Cl_2_150K01_00325.cbf'
    # cf = CifFile.ReadCif('/gpfs/exfel/data/scratch/jhoersch/instrument-geometry-info/Tools/010_Ni_dppe_Cl_2_150K01/010_Ni_dppe_Cl_2_150K01_00001.cbf')

    could_read = False
    while not could_read:
        try:
            with suppress_stdout_stderr():
                full_cbf_info = CifFile.ReadCif(container)
            could_read = True
        except CifFile.StarFile.StarError:
            container.remove_char_from_tail()

    # finalcif = CifFile.CifFile(scoping='instance',standard='CIF')
    # finalcif = CifFile.CifFile()
    # with open(text, 'r') as texi:
    #     print('its typ', type(texi))

    # cf = CifFile.StarFile.ReadStar(container, prepared=finalcif,grammar='auto',
    # scantype='standard', permissive=False, CBF=False)

    # print('cf', cf)

    # content = cbf.read(filename, parse_miniheader=True)
    # print('cont', content)

    # print('my meta: \n', content.metadata)
    # print('my mini: \n', content.miniheader)

    return full_cbf_info

