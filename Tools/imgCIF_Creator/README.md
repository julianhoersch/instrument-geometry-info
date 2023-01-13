# The imgCIF_Creator

The creator is a tool to create imgCIF files out of HDF5, full CBF and some common
subset of miniCBF. An interactive command line interface can be used to
collect the information that is missing from the dataset under investigation.

## Installation

The creator does not yet exist as a package and to install it your can clone this repository
and create a python virtual environment in which you install the pacakge:

`python -m venv my_environment_name`

source it with

`source my_environment_name/bin/activate`.

and install it into the environment by navigating your terminal to the package
directory with the `setup.py` file. There type:

`pip install .`

This will install the imgCIF_Creator in `my_environment_name`.

## Usage

If the virtual environment is sourced the imgCIF_Creator can be used by typing:

`creator.py path/to/a/file/or/directory/you/want/to/convert`

You always need to provide the path to the directory containing the scan files, or
in case it is an hdf5 file, to the master file.

You can provide more optional arguments to the call, for example to provide the
constant portion of the file name (`creator.py path --stem my/stem`) if the
programm can not deduce the scan/frame naming convention of the files. All options
are found by typing `creator.py --help`.

The programm will guide you through the process of information collection and if
information is missing in the provided file/directory or if it could not be automatically
retrieved, the input will be requested.

Missing information will be requested by an user prompt and non-required/optional
information can be skipped by an empty input. For some information a specific format
of the input is required and the program will repeatedly ask until the input matches
the format.

## Development

The programm is still in development and might produce faulty output, or crash.

The general structure of the package is build around the `imgCIF_creator.py` module
which contains classes to obtain the required information and generate the imgCIF
file out of that. Depending on the filetype different information extractors
from the `information_extractors` directory are choosen. These must implement the
`extractor_interface.py` and provide the output expected by the `imgCIF_creator.py`.
The `imgCIF_creator.py` checks if the required information is present and will request
an user input if that is not the case.

To generate the imgCIF the PyCIFRW package (https://github.com/jamesrhester/pycifrw)
is used.

Some open questions or not implemented features are marked with #TODO in the code.

## Limitations

Some current limitations are:

- All frames must be in a single local directory
- Frame filenames (for example for cbf) must be of the form `<xxxxx>_ss(0|_)ffff.cbf`
where `<xxxxx>` is arbitrary, `ss` is an arbitrary-length integer
labelling the scan, `ffff` is an arbitrary-length integer
sequentially numbering the frames within the scan starting from
1, and the scan label is followed either by the digit `0` or
an underscore.
- Many variants of axis names may not be recognised. Please
report any that you think should be recognised.

## Testing

There are no tests implemented yet.
