# coding utf8
import setuptools
from rnaseqtools.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="RNASeqTools",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="Some small tools to help with RNA-seq analysis",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/RNASeqTools",

    entry_points={
        "console_scripts": ["RNASeqTools = rnaseqtools.cli:main"]
    },  

    packages=setuptools.find_packages(),

    install_requires=[
        "yxutil",
        "yxseq",
        "yxalign",
        "taxontools",
        "interlap>=0.2.6",
        "numpy>=1.18.1",
        "bcbio-gff>=0.6.6",
        "biopython<=1.80",
        "pyfaidx>=0.5.5.2",        
    ],

    include_package_data=True,

    python_requires='>=3.5',

)