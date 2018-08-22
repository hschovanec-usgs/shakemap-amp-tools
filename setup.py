from distutils.core import setup
from Cython.Build import cythonize
import glob

setup(name='amptools',
      version='0.1dev',
      description='USGS ShakeMap Strong Motion Data Tools',
      include_package_data=True,
      author='Mike Hearne',
      ext_modules=cythonize("pgm/oscillators.pyx"),
      author_email='mhearne@usgs.gov',
      url='',
      packages=['amptools',
                'amptools.io',
                'amptools.io.cwb',
                'amptools.io.dmg',
                'amptools.io.geonet',
                'amptools.io.knet',
                'amptools.io.cosmos',
                'amptools.io.smc',
                'amptools.io.obspy',
                'amptools.io.usc',
                'pgm',
                'pgm.imt',
                'pgm.imc'],
      package_data={
          'amptools': glob.glob('amptools/io/*.csv') + glob.glob('tests/data/*/*'),
          'pgm': glob.glob('amptools/io/*.csv') + glob.glob('tests/data/*/*')
      },
      scripts=['bin/amps2xml',
               'bin/ftpfetch',
               'bin/sm2xml',
               'bin/fdsnfetch',
               'bin/ingvfetch']
      )
