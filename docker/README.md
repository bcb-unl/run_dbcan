# run_dbcan Docker Image

Container image for [run_dbcan](https://github.com/bcb-unl/run_dbcan), built with **micromamba** (conda) for system tools and **pip** for the `dbcan` package from the **same Git ref** as the release. This avoids bioconda indexing delay and keeps the image aligned with PyPI/GitHub tags.

## What gets installed

| Component | Source |
|-----------|--------|
| Python 3.10, DIAMOND | `docker/environment.deps.yml` (conda-forge / bioconda) |
| `dbcan` and Python dependencies | `pip install .` from repository checkout at build time |

`docker/environment.yml` is kept for backward compatibility and matches the conda-only deps file (no `dbcan` conda package).

## Build locally

From the repository root:

```bash
git fetch --tags
docker build -f docker/Dockerfile \
  --build-arg BUILD_VERSION=5.2.9 \
  -t run_dbcan:5.2.9 .
```

`BUILD_VERSION` must match the release you are building (without the `v` prefix). It is passed to `SETUPTOOLS_SCM_PRETEND_VERSION` when `.git` is unavailable in the build context.

Verify:

```bash
docker run --rm run_dbcan:5.2.9 run_dbcan --help
docker run --rm run_dbcan:5.2.9 python -c "import importlib.metadata as m; print(m.version('dbcan'))"
```

## CI / GHCR

- **Release published**: builds and pushes `ghcr.io/bcb-unl/run_dbcan:<version>` and `latest`.
- **workflow_dispatch**: Actions → *Build and Publish dbCAN Docker Image* → Run workflow (choose tag/branch and whether to push).

After build, CI checks that `importlib.metadata.version('dbcan')` equals the release version.

## Usage

```bash
docker pull ghcr.io/bcb-unl/run_dbcan:5.2.9
docker run --rm -it ghcr.io/bcb-unl/run_dbcan:5.2.9 run_dbcan --help
```

## Features

- Based on `micromamba:2-ubuntu22.04`
- `linux/amd64` and `linux/arm64`

## License

MIT License.
