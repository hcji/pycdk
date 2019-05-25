from setuptools import setup, find_packages

setup(name='pycdk',
      version='0.0.1_alpha',
      description='A Python wrapper for the CDK',
      license='AGPLv3',
      author='Ji Hongchao',
      author_email='ji.hongchao@foxmail.com',
      url='https://github.com/hcji/PyFingerprint',
	  include_package_data = True,
      packages=find_packages()
     )