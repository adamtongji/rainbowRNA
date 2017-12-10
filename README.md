# rainbowRNA
## RNAseq网页自动化流程

**ONLY Browser Support: Firefox** <br>
**唯一支持浏览器: 火狐**<br>

### Installation
1. Download and install miniconda, and configure the miniconda environment.<br>
```Shell
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
# then add conda environ to you local enviroment.
# Please add:
# export PATH=PATH_TO_CONDA_BIN:$PATH
# to ~/.bashrc, for example

(conda config --add channels r)
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```
2. Download the source code to local disk from the git directory <br>
```Shell
git clone https://github.com/adamtongji/rainbowRNA
```
3. Install CIRI2  if you want to run circRNA module <br>

- CIRI2: ["https://sourceforge.net/projects/ciri/files/CIRI2/"](https://sourceforge.net/projects/ciri/files/CIRI2/)

4. Install the rainbowRNA <br>
```Shell
bash install.sh
```
If occurs with following error: <br>
```Python
if not line or line.startswith('#'):
   UnicodeDecodeError:
   'ascii' codec can't decode byte 0xe7 in position 50: ordinal not in range(128)
```
Please add following code to the script file `"YOUPATH"/miniconda2/lib/python2.7/site-packages/conda/cli/common.py` <br>
```Python
reload(sys)
sys.setdefaultencoding('utf8')
```
Just add these code below  `"import sys"` to avoid encoding error.<br>

5. (Recommend) Install the genome and annotation files from Illumina iGenomes <br>
The igenome website: **[iGenoms download page](https://support.illumina.com/sequencing/sequencing_software/igenome.html)**. <br>