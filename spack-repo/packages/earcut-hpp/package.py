from spack.package import *


class EarcutHpp(Package):
    """Mapbox's header-only polygon triangulation library."""

    homepage = "https://github.com/mapbox/earcut.hpp"
    url = "https://codeload.github.com/mapbox/earcut.hpp/tar.gz/refs/tags/v2.2.4"

    license("ISC")
    version(
        "2.2.4",
        sha256="fcfa6a47a52d4c94dc960bdb747f17e077609235517b0bb5ce8097d6b747695a",
        extension="tar.gz",
    )

    def install(self, spec, prefix):
        install_tree("include", prefix.include)
        mkdirp(prefix.share.licenses.earcut_hpp)
        install("LICENSE", prefix.share.licenses.earcut_hpp)
