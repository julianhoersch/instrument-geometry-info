import abc

class ExtractorInterface(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get_source_info') and
                callable(subclass.get_facility_info) and
                hasattr(subclass, 'get_axes_info') and
                callable(subclass.get_axes_info) and
                # hasattr(subclass, 'get_array_info') and
                # callable(subclass.get_array_info) and
                hasattr(subclass, 'get_detector_info') and
                callable(subclass.get_detector_info) and
                hasattr(subclass, 'get_wavelength_info') and
                callable(subclass.get_wavelength_info) and
                hasattr(subclass, 'get_scan_setting_info') and
                callable(subclass.get_scan_setting_info) #and
                # hasattr(subclass, 'get_scan_info') and
                # callable(subclass.get_scan_info) and
                # hasattr(subclass, 'get_step_info') and
                # callable(subclass.get_step_info) and
                # hasattr(subclass, 'get_array_info') and
                # callable(subclass.get_array_info) and
                # hasattr(subclass, 'get_ids_info') and
                # callable(subclass.get_ids_info) and
                # hasattr(subclass, 'get_external_ids_info') and
                # callable(subclass.get_external_ids_info)
                or NotImplemented)


    @abc.abstractmethod
    def get_source_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_axes_info():

        raise NotImplementedError

    # @abc.abstractmethod
    # def get_array_info():

    #     raise NotImplementedError

    @abc.abstractmethod
    def get_detector_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_wavelength_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_scan_settings_info():

        raise NotImplementedError

    # @abc.abstractmethod
    # def get_scan_info():

    #     raise NotImplementedError

    # @abc.abstractmethod
    # def get_step_info():

    #     raise NotImplementedError

    # @abc.abstractmethod
    # def get_array_info():

    #     raise NotImplementedError

    # @abc.abstractmethod
    # def get_ids_info():

    #     raise NotImplementedError

    # @abc.abstractmethod
    # def get_external_ids_info():

    #     raise NotImplementedError


    # @abc.abstractmethod
    # def load_data_source(self, path: str, file_name: str):
    #     """Load in the data set"""
    #     raise NotImplementedError

    # @abc.abstractmethod
    # def extract_text(self, full_file_path: str):
    #     """Extract text from the data set"""
    #     raise NotImplementedError








# class PdfParserNew(FormalParserInterface):
#     """Extract text from a PDF."""
#     def load_data_source(self, path: str, file_name: str) -> str:
#         """Overrides FormalParserInterface.load_data_source()"""
#         pass

#     def extract_text(self, full_file_path: str) -> dict:
#         """Overrides FormalParserInterface.extract_text()"""
#         pass

# class EmlParserNew(FormalParserInterface):
#     """Extract text from an email."""
#     def load_data_source(self, path: str, file_name: str) -> str:
#         """Overrides FormalParserInterface.load_data_source()"""
#         pass

#     def extract_text_from_email(self, full_file_path: str) -> dict:
#         """A method defined only in EmlParser.
#         Does not override FormalParserInterface.extract_text()
#         """
#         pass
