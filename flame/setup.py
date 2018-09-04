from setuptools import setup

setup(name='flame',
      version='0.2',
      description='A flake lattice modeler and growth simulation.',
      url='https://github.com/gsec/flame',
      author='Guilherme Stein',
      author_email='gui@posteo.net',
      license='GPL',
      packages=['flame'],
      zip_safe=False,
      test_suite='py.test',
      tests_require=['py.test'],
      entry_points={
          'console_scripts': ['flame=flame.cli:main']}
      )
