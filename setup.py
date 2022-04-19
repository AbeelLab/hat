from setuptools import setup, find_packages

setup(
    name='HAT',
    version='0.0.1',
    packages=['hat'],
    url='',
    license='GPL-3.0-or-later',
    author='Ramin Shirali Hossein Zade',
    author_email='r.shirali.hz@gmail.com',
    entry_points={'console_scripts': ['HAT = hat.HAT_main:main']},
    description='HAT:‌  Haplotype assembly tool that use both long and short reads to reconstruct haplotypes'
)
