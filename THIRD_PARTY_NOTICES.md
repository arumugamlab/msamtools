# Third-party software

msamtools uses third-party components from **HTSlib** and **ZOE/SNAP**.

The msamtools source distribution bundles modified versions of two files
from ZOE source code. The upstream license text is distributed as well as
installed alongside msamtools as:

- `zoe-LICENSE`

The msamtools source distribution does **not** bundle HTSlib source code.
During the build process, the required upstream HTSlib version is downloaded
and statically linked into the msamtools executable. The corresponding
upstream license text is installed alongside msamtools as:

- `htslib-LICENSE`

These files are installed together with:

- `LICENSE`
- `THIRD_PARTY_NOTICES.md`

typically under:

- `share/licenses/msamtools/`

The msamtools Docker image additionally includes the **samtools** executable.
Samtools is not used to build or link the msamtools executable. The Docker
image uses the samtools release corresponding to the HTSlib version used by
msamtools. Its license is provided as:

- `samtools-LICENSE`

## HTSlib (v1.24)

HTSlib is developed by Genome Research Ltd. and the HTSlib authors and
contributors.

HTSlib is available at
[https://github.com/samtools/htslib](https://github.com/samtools/htslib)
and distributed under the MIT/Expat License.

HTSlib is downloaded at build time and statically linked into msamtools.
The complete copyright and license notice for the version used to build
msamtools is provided in `htslib-LICENSE`.

## samtools (v1.24; Docker image only)

Samtools is developed by Genome Research Ltd. and the samtools authors and
contributors.

Samtools is available at
[https://github.com/samtools/samtools](https://github.com/samtools/samtools)
and distributed under the MIT/Expat License.

Samtools is not a build or API dependency of msamtools. The msamtools Docker
image additionally distributes the samtools executable for user convenience.
The Docker image uses the samtools release corresponding to the HTSlib version
used by msamtools. Its complete copyright and license notice is provided in
`samtools-LICENSE`.

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
