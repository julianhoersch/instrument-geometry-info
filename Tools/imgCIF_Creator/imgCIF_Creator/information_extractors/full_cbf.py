import CifFile
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

    with open(filename, 'rb') as file:

        file_content = file.read()
        header_end_mark = b'\x0C\x1A\x04\xD5'
        header_end_index = file_content.find(header_end_mark)
        text = file_content[:header_end_index].decode('utf-8')
        text = text.split("--CIF-BINARY-FORMAT-SECTION--")[0]

    container = text_container(text)
    could_read = False
    while not could_read:
        try:
            with suppress_stdout_stderr():
                full_cbf_info = CifFile.ReadCif(container)
            could_read = True
        except CifFile.StarFile.StarError:
            container.remove_char_from_tail()

    return full_cbf_info
