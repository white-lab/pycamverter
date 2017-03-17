# -*- mode: python -*-

import pymzml

block_cipher = None

def get_proteowizard_version():
    with open(os.path.join("pycamv", "proteowizard.py")) as f:
        for line in f:
            if line.startswith("PROTEOWIZARD_VERSION"):
                return line.split("=")[1].strip('" ')

def dir_files(path, rel):
    ret = []
    for p,d,f in os.walk(path):
        relpath = p.replace(path, '')[1:]
        for fname in f:
            ret.append((os.path.join(rel, relpath, fname),
                        os.path.join(p, fname), 'DATA'))
    return ret

a = Analysis(
        [os.path.join('pycamv', 'main.py')],
        pathex=[],
        binaries=None,
        datas=[
            (
                'pycamv/ProteoWizard/*',
                'pycamv/ProteoWizard/ProteoWizard ' +
                get_proteowizard_version() +
                '/',
            ),
        ],
        hiddenimports=[],
        hookspath=[],
        runtime_hooks=[],
        excludes=[],
        win_no_prefer_redirects=False,
        win_private_assemblies=False,
        cipher=block_cipher,
)

a.datas.extend(dir_files(os.path.join(os.path.dirname(pymzml.__file__),
    'obo'), os.path.join('pymzml', 'obo')))


pyz = PYZ(
        a.pure,
        a.zipped_data,
        cipher=block_cipher,
)

exe = EXE(
        pyz,
        a.scripts,
        a.binaries,
        a.zipfiles,
        a.datas,
        name='PyCAMVerter',
        debug=True,
        strip=False,
        upx=True,
        console=True,
        exclude_binaries=False,
)

coll = COLLECT(
        exe,
        strip=False,
        upx=True,
        name='PyCAMVerter',
)
