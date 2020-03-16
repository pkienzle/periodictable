import sys

# Doctest needs doc/sphinx on the path so that discoverer package can be found.
# Using the pytest-pythonpath package, this would be:
#    setup.py: test_require=[..., 'pytest-pythonpath']
#    pytest.ini:  python_paths = "doc/sphinx"
# Not doing so because it introduces another dependency.
def pytest_configure(config):
    sys.path.insert(0, "doc/sphinx")
