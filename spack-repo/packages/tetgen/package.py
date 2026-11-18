import os

from spack.package import *


class Tetgen(Package):
    """TetGen executable and static library, using the project's pinned source."""

    homepage = "https://github.com/TetGen/TetGen"
    url = "https://codeload.github.com/TetGen/TetGen/tar.gz/refs/tags/v1.6.0"

    license("AGPL-3.0-or-later")
    version(
        "1.6.0",
        sha256="524bddc38afc205db075b74b67da424b87e676cfd12eccf7915d2e945afcac10",
        extension="tar.gz",
    )

    depends_on("cxx", type="build")
    depends_on("gmake", type="build")

    def install(self, spec, prefix):
        # Use Spack's compiler wrapper instead of upstream's hardcoded g++.
        # Keep robust predicates at upstream's -O0; both objects need PIC.
        make(
            "tetgen", "tetlib",
            "CXX=" + os.environ["CXX"],
            "CXXFLAGS=-O3 -fPIC",
            "PREDCXXFLAGS=-O0 -fPIC",
        )
        mkdirp(prefix.bin, prefix.include, prefix.lib, prefix.share.licenses.tetgen)
        install("tetgen", prefix.bin)
        install("tetgen.h", prefix.include)
        install("libtet.a", prefix.lib)
        install("LICENSE", prefix.share.licenses.tetgen)
