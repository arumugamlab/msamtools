# Third-party software

msamtools uses third-party components from **samtools**, **HTSlib** and **ZOE/SNAP**.

The msamtools source distribution bundles modified versions of two files
from ZOE source code. The upstream license text is distributed as well as
installed alongside msamtools as:
- `zoe-LICENSE`

The msamtools source distribution does **not** bundle the samtools or HTSlib
source code. During the build process, the required upstream versions are
downloaded and used to build msamtools.

Binary distributions of msamtools may therefore contain code derived from or
statically linked with samtools and HTSlib. For such distributions, the
corresponding upstream license texts are installed alongside msamtools as:

- `samtools-LICENSE`
- `htslib-LICENSE`

These files are installed together with:

- `LICENSE`
- `THIRD_PARTY_NOTICES.md`

typically under:

- `share/licenses/msamtools/`

The msamtools Docker image additionally includes the `samtools` executable.
Its license is covered by the accompanying `samtools-LICENSE` file.

## samtools (v1.9)

samtools is developed by Genome Research Ltd. and the samtools authors and
contributors.

samtools is available at
[https://github.com/samtools/samtools](https://github.com/samtools/samtools)
and distributed under the MIT/Expat License.

The complete copyright and license notice for the version used to build
msamtools is provided in `samtools-LICENSE`.

## HTSlib (v1.9)

HTSlib is developed by Genome Research Ltd. and the HTSlib authors and
contributors.

HTSlib is available at
[https://github.com/samtools/htslib](https://github.com/samtools/htslib)
distributed under the MIT/Expat License.

The complete copyright and license notice for the version used to build
msamtools is provided in `htslib-LICENSE`.

## ZOE / SNAP

ZOE / SNAP is developed by Ian Korf.

`zoeTools.c` and `zoeTools.h` are modified, reduced versions of files from
the ZOE library distributed with SNAP.

Upstream source:
[https://github.com/KorfLab/SNAP](https://github.com/KorfLab/SNAP)

The versions included in msamtools were re-derived from SNAP commit
`4ad1e957cd8e68b63857cc1cb3380d39a7b518b1`.

The upstream ZOE/SNAP code is distributed under the MIT License. The complete
upstream copyright and license notice is provided in `zoe-LICENSE`.
