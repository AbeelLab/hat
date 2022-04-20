from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='HAT_phasing',
    version='0.0.2',
    packages=['hat'],
    url='',
    license='GPL-3.0-or-later',
    author='Ramin Shirali Hossein Zade',
    author_email='r.shirali.hz@gmail.com',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['pysam', 'numpy', 'seaborn', 'matplotlib', 'biopython'],
    entry_points={'console_scripts': ['HAT = hat.HAT_main:main']},
    python_requires=">=3.6",
    description='HAT:‌  Haplotype assembly tool that use both long and short reads to reconstruct haplotypes'
)
