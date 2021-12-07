import os
from setuptools import setup, Extension, find_packages
from distutils.sysconfig import customize_compiler
import Cython.Build
from Cython.Build import cythonize


def CppExtension(name, sources):
    return Extension(
        name,
        sources=sources,
        language="c++",
        extra_compile_args=["-std=c++11", "-Werror=return-type", "-Werror=narrowing"],
        undef_macros=["NDEBUG"],
    )


extensions = [
    CppExtension(
        "whatshap.core",
        sources=[
            "whatshap/core.pyx",
            "src/pedigree.cpp",
            "src/pedigreedptable.cpp",
            "src/pedigreecolumncostcomputer.cpp",
            "src/columnindexingiterator.cpp",
            "src/columnindexingscheme.cpp",
            "src/entry.cpp",
            "src/graycodes.cpp",
            "src/read.cpp",
            "src/readset.cpp",
            "src/columniterator.cpp",
            "src/indexset.cpp",
            "src/genotype.cpp",
            "src/binomial.cpp",
            "src/pedigreepartitions.cpp",
            "src/phredgenotypelikelihoods.cpp",
            "src/genotyper.cpp",
            "src/genotypedistribution.cpp",
            "src/genotypedptable.cpp",
            "src/genotypecolumncostcomputer.cpp",
            "src/backwardcolumniterator.cpp",
            "src/transitionprobabilitycomputer.cpp",
            "src/hapchat/basictypes.cpp",
            "src/hapchat/balancedcombinations.cpp",
            "src/hapchat/binomialcoefficient.cpp",
            "src/hapchat/hapchatcore.cpp",
            "src/hapchat/hapchatcolumniterator.cpp",
            "src/polyphase/clustereditingsolution.cpp",
            "src/polyphase/clustereditingsolver.cpp",
            "src/polyphase/edgeheap.cpp",
            "src/polyphase/inducedcostheuristic.cpp",
            "src/polyphase/staticsparsegraph.cpp",
            "src/polyphase/switchflipcalculator.cpp",
            "src/polyphase/trianglesparsematrix.cpp",
            "src/polyphase/readscoring.cpp",
            "src/polyphase/haplothreader.cpp",
        ],
    ),
    CppExtension("whatshap.priorityqueue", sources=["whatshap/priorityqueue.pyx"]),
    CppExtension("whatshap.align", sources=["whatshap/align.pyx"]),
    CppExtension("whatshap._variants", sources=["whatshap/_variants.pyx"]),
]


class BuildExt(Cython.Build.build_ext):
    def build_extensions(self):
        # Remove the warning about “-Wstrict-prototypes” not being valid for C++,
        # see http://stackoverflow.com/a/36293331/715090
        customize_compiler(self.compiler)
        if self.compiler.compiler_so[0].endswith("clang"):
            # Clang needs this option in order to find the unordered_set header
            print("detected clang, using option -stdlib=libc++")
            self.compiler.compiler_so.append("-stdlib=libc++")
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        super().build_extensions()


with open("doc/README.rst", encoding="utf-8") as f:
    long_description = f.read()


# Avoid compilation if we are being installed within Read The Docs
if os.environ.get("READTHEDOCS") == "True":
    cmdclass = {}
    ext_modules = []
    install_requires = []
else:
    cmdclass = {"build_ext": BuildExt}
    ext_modules = extensions
    install_requires = [
        "pysam>=0.18.0",
        "pyfaidx>=0.5.5.2",
        "networkx",
        "biopython>=1.73",  # pyfaidx needs this for reading bgzipped FASTA files
        "scipy",
        "xopen>=1.2.0",
        "dataclasses; python_version<'3.7'",
    ]

setup(
    name="whatshap",
    use_scm_version={"write_to": "whatshap/_version.py"},
    author="WhatsHap authors",
    author_email="whatshap@cwi.nl",
    url="https://github.com/whatshap/whatshap/",
    description="phase genomic variants using DNA sequencing reads",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    license="MIT",
    cmdclass=cmdclass,
    ext_modules=cythonize(ext_modules),
    packages=find_packages(),
    entry_points={"console_scripts": ["whatshap = whatshap.__main__:main"]},
    install_requires=install_requires,
    extras_require={"dev": ["Cython", "pytest", "sphinx", "sphinx_issues", "pysam-stubs"]},
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Cython",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
