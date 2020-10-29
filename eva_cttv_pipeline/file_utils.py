import argparse
import gzip
import sys
import importlib
import importlib._bootstrap
import importlib.util
import os


def open_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def get_resource_file(package, resource):
    spec = importlib.util.find_spec(package)
    if spec is None:
        return None
    mod = (sys.modules.get(package) or importlib._bootstrap._load(spec))
    if mod is None or not hasattr(mod, '__file__'):
        return None

    parts = resource.split('/')
    parts.insert(0, os.path.dirname(mod.__file__))
    resource_name = os.path.join(*parts)
    return resource_name
