from distutils.core import setup
setup(
  name = 'Chemical Pressure Plotting',         # How you named your package folder (MyLib)
  packages = ['Chemical Pressure Plotting'],   # Chose the same as "name"
  version = '0.1',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Uses python to plot chemical pressures of atomic compounds',   # Give a short description about your library
  author = 'Alexander Smits',                   # Type in your name
  author_email = 'ajsmits2@wisc.edu',      # Type in your E-Mail
  url = 'https://github.com/ajsmits2/Chemical-Pressure/',   # Provide either the link to your github or to your website
  download_url = https://github.com/ajsmits2/Chemical-Pressure/archive/refs/tags/Chemistry.tar.gz,    # I explain this later on
  keywords = ['Chemistry', 'Pressure', 'Chemical'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'numpy',
          'pyvista',
      ],
  classifiers=[
    'Development Status :: 5 - Production',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package

    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',

    'License :: OSI Approved :: MIT License',   # Again, pick a license

    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
