
from setuptools import setup, find_packages
# Warning : do not import the distutils extension before setuptools
# It does break the cythonize function calls
from distutils.extension import Extension
from distutils.sysconfig import get_config_vars

import glob
import os
import platform
import sys

from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy

DEBUG = False

SUPPORT_CODE_INCLUDE = './wrapper'

AT_LIBRARY = 'ActiveTickServerAPI'
AT_LIBRARY_DIR = './sdk/bin'
AT_INCLUDE_DIR = './sdk/include'

# FIXME: would be good to be able to customize the path with environment
# variables in place of hardcoded paths ...
if sys.platform == 'darwin':
    INCLUDE_DIRS = ['/usr/local/include', '.', '../sources/boost_1_55_0']
    LIBRARY_DIRS = ["/usr/local/lib"]
    ## From SO: hack to remove warning about strict prototypes
    ## http://stackoverflow.com/questions/8106258/cc1plus-warning-command-line-option-wstrict-prototypes-is-valid-for-ada-c-o
    (opt,) = get_config_vars('OPT')
    os.environ['OPT'] = " ".join(
        flag for flag in opt.split() if flag != '-Wstrict-prototypes')
elif sys.platform == 'win32':
    INCLUDE_DIRS = [ r'c:\dev\boost_1_56_0', '.' ]
    LIBRARY_DIRS = [ '.', r'.\dll', ]
elif sys.platform == 'linux2':
    INCLUDE_DIRS = ['/usr/local/include', '/usr/include', '.']
    LIBRARY_DIRS = ['/usr/local/lib', '/usr/lib', ]

INCLUDE_DIRS.append(SUPPORT_CODE_INCLUDE)
# INCLUDE_DIRS.append(AT_INCLUDE_DIR)
INCLUDE_DIRS.append(numpy.get_include())
# LIBRARY_DIRS.append(AT_LIBRARY_DIR)

def get_define_macros():
    #defines = [ ('HAVE_CONFIG_H', None)]
    defines = []
    if sys.platform == 'win32':
        # based on the SWIG wrappers
        defines += [
            (name, None) for name in [
                '__WIN32__', 'WIN32', 'NDEBUG', '_WINDOWS', 'NOMINMAX', 'WINNT',
                '_WINDLL', '_SCL_SECURE_NO_DEPRECATE', '_CRT_SECURE_NO_DEPRECATE',
                '_SCL_SECURE_NO_WARNINGS'
            ]
        ]
    return defines

def get_extra_compile_args():
    if sys.platform == 'win32':
        args = ['/GR', '/FD', '/Zm250', '/EHsc']
        if DEBUG:
            args.append('/Z7')
    else:
        args = ['-std=c++11', '-Wno-unused-function']

    return args

def get_extra_link_args():
    if sys.platform == 'win32':
        args = ['/subsystem:windows', '/machine:I386']
        if DEBUG:
            args.append('/DEBUG')
    elif sys.platform == 'darwin':
        major, minor = [
            int(item) for item in platform.mac_ver()[0].split('.')[:2]]
        if major == 10 and minor >= 9:
            # On Mac OS 10.9 we link against the libstdc++ library.
            args = ['-stdlib=libstdc++', '-mmacosx-version-min=10.6']
        else:
            args = []
    else:
        args = []

    return args

# http://docs.cython.org/src/reference/compilation.html#compiler-directives
CYTHON_DIRECTIVES = {"embedsignature": True, "gdb_debug": DEBUG}

def collect_extensions():
    """ Collect all the directories with Cython extensions and return the list
    of Extension.

    Th function combines static Extension declaration and calls to cythonize
    to build the list of extenions.
    """

    kwargs = {
        'language':'c++',
        'include_dirs':INCLUDE_DIRS,
        'library_dirs':LIBRARY_DIRS,
        'define_macros':get_define_macros(),
        'extra_compile_args':get_extra_compile_args(),
        'extra_link_args':get_extra_link_args(),
        'libraries':[AT_LIBRARY],
        'cython_directives':CYTHON_DIRECTIVES
    }

    test_extension = Extension('pet.pet',
        ['pet/pet.pyx', 'lib/patch.cxx'],
        **kwargs
    )

    manual_extensions = [
        test_extension,
    ]

    cython_extension_directories = []
    collected_extensions = []
    # for dirpath, directories, files in os.walk('activetick'):

    #     # skip the settings package
    #     if dirpath.find('settings') > -1 or dirpath.find('test') > -1:
    #         continue

    #     # if the directory contains pyx files, cythonise it
    #     if len(glob.glob('{0}/*.pyx'.format(dirpath))) > 0:
    #         cython_extension_directories.append(dirpath)

    # collected_extensions = cythonize(
    #     [
    #         Extension('*', ['{0}/*.pyx'.format(dirpath)], **kwargs)
    #         for dirpath in cython_extension_directories
    #     ]
    # )

    # # remove  all the manual extensions from the collected ones
    # names = [extension.name for extension in manual_extensions]
    # for ext in collected_extensions:
    #     if ext.name in names:
    #         collected_extensions.remove(ext)
    #         continue

    extensions = collected_extensions + manual_extensions

    return cythonize(extensions, gdb_debug=DEBUG)

from pet import version as petversion
setup(
    name = 'pet',
    version = petversion.full,
    author = 'Jan Paral',
    license = 'GNU',
    packages = find_packages(),
    ext_modules = collect_extensions(),
    cmdclass = {'build_ext': build_ext},
    install_requires = ['distribute', 'cython'],
    zip_safe = False
)
