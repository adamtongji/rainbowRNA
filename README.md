# rainbowRNA
## RNAseq网页自动化流程

**ONLY Browser Support: Firefox**
**唯一支持浏览器: 火狐**

### Installation
1. Download and install miniconda, and configure the miniconda environment.
`
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
`
2. Download the source code to local disk from the git directory
`
git clone https://github.com/adamtongji/rainbowRNA
`
3. Install CIRI2 and targetscan(including lib) if you want to run circRNA module

CIRI2: ["https://sourceforge.net/projects/ciri/files/CIRI2/"](https://sourceforge.net/projects/ciri/files/CIRI2/)
Targetscan: ["http://www.targetscan.org/vert_71/vert_71_data_download/targetscan_70.zip"](http://www.targetscan.org/vert_71/vert_71_data_download/targetscan_70.zip)

4. Install the rainbowRNA
`
bash install.sh
`
If occurs with following error:
`
if not line or line.startswith('#'):
   UnicodeDecodeError:
   'ascii' codec can't decode byte 0xe7 in position 50: ordinal not in range(128)
`
Please add following code to the script file `"YOUPATH"/miniconda2/lib/python2.7/site-packages/conda/cli/common.py`
`
reload(sys)
sys.setdefaultencoding('utf8')
`
Just add these code below  `import sys` to avoid encoding error.

5. (Recommend) Install the genome and annotation files from Illumina iGenomes
The igenome website: **[links](https://support.illumina.com/sequencing/sequencing_software/igenome.html)**.