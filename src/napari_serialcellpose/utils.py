from importlib.metadata import PackageNotFoundError, version

def get_cp_version():
    try: 
        __cp_version__ = int(version("cellpose").split('.')[0])
        if __cp_version__ > 3:
            print(f'>> Initializing with CellposeSAM:')
        else:
            print(f'>> Initializing with Cellpose legacy release (v{__cp_version__}).')
    except PackageNotFoundError: 
        __cp_version__ = "uninstalled" 
        print('>> Cellpose not found. Please install Cellpose to use napari-serialcellpose.')
    
    return __cp_version__