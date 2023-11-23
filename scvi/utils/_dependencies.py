import importlib


def error_on_missing_dependencies(*modules):
    missing_modules = []
    for module in modules:
        try:
            importlib.import_module(module)
        except ImportError:
            missing_modules.append(module)
    if len(missing_modules) > 0:
        raise ModuleNotFoundError(
            f"Please install {missing_modules} to use this functionality."
        )
