from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os

os.environ["CC"] = "g++-4.7"
os.environ["CXX"] = "g++-4.7"

srcDir = "../src"
src = ["NetworKit.pyx"]	# list of source files
# include all .cpp files
# for (root, dirs, files) in os.walk(srcDir, topdown=False):
# 	for name in files:
# 		if name.endswith(".cpp") and (not root.endswith("/test")):
# 			src.append(os.path.join(root, name))
			
print("source files: {0}".format(src))

modules = [Extension("NetworKit",
					src,
					language = "c++",
					extra_compile_args=["-fopenmp", "-std=c++11", "-DNOLOG4CXX", "-DNOGTEST"],
					extra_link_args=["-fopenmp", "-std=c++11"],
					libraries=["NetworKit-Core-O"],
					library_dirs=["../NetworKit-Core-O/"])]

for e in modules:
	e.cython_directives = {"embedsignature" : True}

setup(name="NetworKit",
	 cmdclass={"build_ext": build_ext},
	 ext_modules=modules)

print("setup script terminated")

