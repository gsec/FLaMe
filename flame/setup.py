from setuptools import setup

setup(name='flame',
      version='0.1',
      description='A flake lattice modeler and growth simulation.',
      url='https://github.com/gsec/thesis',
      author='Guilherme Stein',
      author_email='gui@posteo.net',
      license='GPL',
      packages=['flame'],
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'],
      entry_points={
          'console_scripts': ['flame=flame.cli:main']}
      )