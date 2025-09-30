### Set LibTorch in Spack.
I'm using the corresponding version to CUDA 12.6.

1. Download the `libtorch` package.
```bash
cd ~
wget "https://download.pytorch.org/libtorch/cu126/libtorch-shared-with-deps-2.8.0%2Bcu126.zip"
unzip libtorch-shared-with-deps-2.8.0%2Bcu126.zip
rm libtorch-shared-with-deps-2.8.0%2Bcu126.zip
```
2. Convert the `/libtorch` to a Spack repo.
```bash
cd libtorch
echo -e "repo:\n  namespace: libtorch" > repo.yaml
```
```bash
mkdir packages/libtorch
cd packages/libtorch
```
```bash
cat > package.py << 'EOF'
from spack.package import *
import shutil
import os

class Libtorch(BundlePackage):
    homepage = "https://pytorch.org"
    version('2.8.0_cu126')
    
    def install(self, spec, prefix):
        src_path = os.path.expanduser('~/libtorch')
        if not os.path.exists(src_path):
            raise InstallError(f"Local libtorch folder not found: {src_path}")
        
        os.makedirs(prefix, exist_ok=True)
        for item in os.listdir(src_path):
            s = os.path.join(src_path, item)
            d = os.path.join(prefix, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)
EOF
```
3. Link the repo
```bash
spack repo add libtorch
```
4. Run `create-env.sh`