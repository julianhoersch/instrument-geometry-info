import abc

class ExtractorInterface(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get_source_info') and
                callable(subclass.get_facility_info) and
                hasattr(subclass, 'get_axes_info') and
                callable(subclass.get_axes_info) and
                hasattr(subclass, 'get_array_info') and
                callable(subclass.get_array_info) and
                hasattr(subclass, 'get_detector_info') and
                callable(subclass.get_detector_info) and
                hasattr(subclass, 'get_radiation_info') and
                callable(subclass.get_radiation_info) and
                hasattr(subclass, 'get_scan_setting_info') and
                callable(subclass.get_scan_setting_info) and
                hasattr(subclass, 'get_scan_settings_info') and
                callable(subclass.get_scan_settings_info) and
                hasattr(subclass, 'get_uncategorized_info') and
                callable(subclass.get_uncategorized_info) and
                hasattr(subclass, 'get_array_info') and
                callable(subclass.get_array_info) or
                NotImplemented)


    @abc.abstractmethod
    def get_source_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_axes_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_array_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_detector_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_radiation_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_scan_settings_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_scan_settings_info():

        raise NotImplementedError

    @abc.abstractmethod
    def get_uncategorized_info():

        raise NotImplementedError
