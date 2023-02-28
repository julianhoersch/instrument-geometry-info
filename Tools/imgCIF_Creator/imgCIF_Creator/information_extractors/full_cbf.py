"""Functionalities to extract information from full cbf headers.
"""

from os import devnull
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import CifFile

@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull. No output."""

    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


class TextContainer:
    """See documentation of the init method.
    """

    def __init__(self, text) -> None:
        """A container class for text that provides a read method.

        Args:
            text (str): the text the class should hold
        """
        self.text = text

    def remove_char_from_tail(self):
        """Remove the last char in the text and set as new text.
        """
        self.text = self.text[:-1]

    def read(self):
        """return the text.

        Returns:
            str: the text
        """
        return self.text


def extract_full_cbf_header_information(filename):
    """Extract the information from the full cbf header using the CifFile module.

    Args:
        filename (str): the filename of the cbf file

    Returns:
        CifFile.CifFile_module.CifFile: the full header information in the format
            of the CifFile module. Entries can be accessed like in a dictionary
    """

    with open(filename, 'rb') as file:
        file_content = file.read()
        header_end_mark = b'\x0C\x1A\x04\xD5'
        header_end_index = file_content.find(header_end_mark)
        text = file_content[:header_end_index].decode('utf-8')
        text = text.split("--CIF-BINARY-FORMAT-SECTION--")[0]

    container = TextContainer(text)
    could_read = False
    while not could_read:
        try:
            with suppress_stdout_stderr():
                full_cbf_info = CifFile.ReadCif(container)
            could_read = True
        except CifFile.StarFile.StarError:
            # remove the last character of the text until the CifFile module can
            # read the file, this is not super fast...
            container.remove_char_from_tail()

    return full_cbf_info
