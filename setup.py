import os
import pathlib
from setuptools import setup,find_packages,find_namespace_packages
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()




extensions = [
    Extension(
        name="convert_to_exr",
        sources=["./src/ToBlender/convert_to_exr.pyx"],
        include_dirs=["./src/common", numpy.get_include()],
        #define_macros=[('NPY_NO_DEPRECATED_API','NPY_1_7_API_VERSION')],
        extra_link_args=["-fopenmp"],
    ),
    Extension(
        "combine",
        ["./src/KGasSimulation/combine.pyx"],
        #define_macros=[('NPY_NO_DEPRECATED_API','NPY_1_7_API_VERSION')],
        include_dirs=["./src/common", numpy.get_include()],
    ),
    Extension(
        "gas_sim",
        ["./src/KGasSimulation/gas_sim.pyx"], 
        #define_macros=[('NPY_NO_DEPRECATED_API','NPY_1_7_API_VERSION')],
        include_dirs=["./src/common", numpy.get_include()],
    ),

    Extension(
        "particle_sim",
        ["./src/KParticleSimulation/particle_sim.pyx","./src/common/cached_array.c"],
        #define_macros=[('NPY_NO_DEPRECATED_API','NPY_1_7_API_VERSION')],
        include_dirs=["./src/common", numpy.get_include()],
    ),

    Extension(
        "combine",
        ["./src/KParticleSimulation/combine.pyx"],
        #define_macros=[('NPY_NO_DEPRECATED_API','NPY_1_7_API_VERSION')],
        include_dirs=["./src/common", numpy.get_include()],
    ),
]



setup(
    name = "JetCreator",
    version = "0.0.1",
    author = "Timo Väisänen, Aalto university",
    author_email = "",
    description = ("JetCreator for Sispo"),
    license = "BSD 2-Clause Simplified License",
    long_description=read('README.md'),
    python_requires='>=3.6',
    packages=find_packages(),
    #package_dir={"": "subprograms"},
    data_files=[
        ('', [
            "./src/KGasSimulation/defaults.in",
            "./src/KParticleSimulation/defaults.in",
            "./src/NSDSGenerator/defaults.in",
            "./src/ToBlender/defaults.in",
        ]),
    ],
    zip_safe=False,
    # Check dependencies
    install_requires=[
        "orekit",
        "numpy",
        "opencv-contrib-python",
        "pyembree",
        "cython",
        "trimesh",
        "noise",
        "scipy",
    ],

    include_package_data=True,
    ext_modules = cythonize(extensions, 
                                compiler_directives={'language_level' : 3}),
    cmdclass={'build_ext': build_ext},
)
