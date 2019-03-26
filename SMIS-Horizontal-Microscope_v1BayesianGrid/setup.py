import sys

from cx_Freeze import setup, Executable

print(sys.path)
base = None
if sys.platform == 'win32':
    base = 'Win32GUI'

includefiles = ['Logo']

options = {
    'build_exe': {
        'includes': 'atexit',
        'packages': 'serial',
        "include_files":includefiles
    }
}

executables = [
    Executable('main.py', base=base)
]

setup(name='simple_Horizontal',
      version='0.1',
      description='Sample cx_Freeze Horizontal script',
      options=options,
      executables=executables
      )
