# -*- mode: python -*-

block_cipher = None

a = Analysis(
        [os.path.join('pycamv', 'main.py')],
        pathex=[],
        binaries=None,
        datas=[('pycamv/ProteoWizard/*', 'pycamv/ProteoWizard/ProteoWizard 3.0.10505')],
        hiddenimports=[],
        hookspath=[],
        runtime_hooks=[],
        excludes=[],
        win_no_prefer_redirects=False,
        win_private_assemblies=False,
        cipher=block_cipher,
)

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
