import abc

class ExtractorInterface(metaclass=abc.ABCMeta):
    """The interface an extractor has to implement.

    Args:
        metaclass (abc.ABCMeta, optional): Abstact base class. Defaults to abc.ABCMeta.

    Raises:
        NotImplementedError: required method not implemented
    """

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

    # mot working
    # @property
    # @abc.abstractproperty
    # def all_frames(self):
    #     """the scan name / scan frame file mapping
    #     Example:
    #     {('scan name', frame number): {'filename': 'ciclohexano3_010001.cbf'}, ...}
    #     """

    # @abc.abstractmethod
    # def __init__(directory_or_filename, stem, palim):
    #     """Return the information about the facility and beamline or the instrument,
    #     model and location. Cif block: _diffrn_source

    #     Returns:
    #         dict: a dictionary containing the information about the source
    #     """

    #     raise NotImplementedError


    @abc.abstractmethod
    def get_source_info():
        """Return the information about the facility and beamline or the instrument,
        model and location. Cif block: _diffrn_source

        Returns:
            dict: a dictionary containing the information about the source
        """

        raise NotImplementedError


    @abc.abstractmethod
    def get_axes_info():
        """Return the information about the axes settings. Cif block: _axis

        Returns:
            dict: a dictionary containing the information about the axes settings
        """

        raise NotImplementedError


    @abc.abstractmethod
    def get_array_info():
        """Return the information about the array. Cif block: _array_structure_list_axis
        and _array_structure_list

        Returns:
            dict: a dictionary containing the information about the array
        """

        raise NotImplementedError


    @abc.abstractmethod
    def get_detector_info():
        """Return the information about the detector. Cif block: _diffrn_detector
        and _diffrn_detector_axis.

        Returns:
            dict: a dictionary containing the information about the detector
        """

        raise NotImplementedError


    @abc.abstractmethod
    def get_radiation_info():
        """Return the information about the wavelength an type of radiation.
        Cif block: _diffrn_radiation and _diffrn_radiation_wavelength

        Returns:
           dict: a dictionary containing the information about the radiation
        """

        raise NotImplementedError


    @abc.abstractmethod
    def get_scan_settings_info():
        """Return the information about the scans, this is a dictionary containing
        the starting point settings of the axes and the details of each scan.

        For example for scan '08':
        {'08': ({'chi': -60.991, 'phi': 110.0, 'detector_2theta': -12.4,
        'omega': -18.679, 'distance': 40.0}, {'frames': 12, 'axis': 'omega',
        'incr': 2.0, 'time': 1800.0, 'start': -40.679, 'range': 24.0,
        'wavelength': 0.560834, 'x_pixel_size': 0.172, 'y_pixel_size': 0.172,
        'mini_header': ['# detector: pilatus100k',.... ],
        ...}

        Returns:
            dict: a dictionary containing the information about the scans
        """
        raise NotImplementedError


    @abc.abstractmethod
    def get_uncategorized_info():
        """Return the information that was found about the doi and the array
        intensities overload.

        Returns:
            dict: a dictionary containing the doi and the array intensities overload
        """

        raise NotImplementedError
