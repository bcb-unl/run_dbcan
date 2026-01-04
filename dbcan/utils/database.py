import requests
from pathlib import Path
import logging
import tarfile
from typing import Dict, Optional
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm
from dbcan.configs.database_config import DBDownloaderConfig
from dbcan.constants.databases_constants import (DATABASES_CAZYME, DATABASES_CGC, COMPRESSED_DBCAN_PUL, DATABASES_CAZYME_S3, DATABASES_CGC_S3)


logger = logging.getLogger(__name__)

class DBDownloader:
    """Download dbCAN databases from HTTP or AWS S3"""

    def __init__(self, config: DBDownloaderConfig):
        """Initialize the database downloader

        Args:
            config: DBDownloaderConfig parameter
        """
        self.config = config
        # Ensure output directory exists
        self.db_dir_path.mkdir(parents=True, exist_ok=True)

    @property
    def db_dir_path(self) -> Path:
        return Path(self.config.db_dir).expanduser().resolve()

    @property
    def databases(self) -> Dict[str, str]:
        """Get database URLs based on source (HTTP or S3)"""
        if self.config.aws_s3:
            logger.debug("Using AWS S3 as download source")
            if self.config.cgc:
                logger.debug("Including CGC-related databases in download list.")
                return {**DATABASES_CAZYME_S3, **DATABASES_CGC_S3}
            logger.debug("CGC-related databases are disabled (cgc=False).")
            return dict(DATABASES_CAZYME_S3)
        else:
            logger.debug("Using HTTP as download source")
            if self.config.cgc:
                logger.debug("Including CGC-related databases in download list.")
                return {**DATABASES_CAZYME, **DATABASES_CGC}
            logger.debug("CGC-related databases are disabled (cgc=False).")
            return dict(DATABASES_CAZYME)

    def download_file(self):
        """Download all databases from HTTP or AWS S3"""
        session = self._prepare_session()
        try:
            for filename, url in self.databases.items():
                output_path = self.db_dir_path / filename

                # If file exists and no_overwrite is set, skip download
                if output_path.exists():
                    # skip if no_overwrite
                    if self.config.no_overwrite:
                        logger.info(f"File {filename} already exists, skipping download (no_overwrite).")
                        continue
                    # skip if resume and size matches
                    remote_size = self._head_content_length(session, url)
                    local_size = output_path.stat().st_size
                    if remote_size is not None and local_size == remote_size:
                        logger.info(f"File {filename} already complete (size match), skipping.")
                        continue

                logger.info(f"Downloading {filename} from {'S3' if self.config.aws_s3 else 'HTTP'}: {url}")
                try:
                    self._download_single_file(session, url, output_path)
                    logger.info(f"{filename} successfully downloaded")

                    # Extract PUL_ALL if it's the tar.gz file
                    if filename == COMPRESSED_DBCAN_PUL:
                        self._extract_tar_file(output_path, self.db_dir_path)

                except Exception as e:
                    logger.error(f"{filename} download error: {e}", exc_info=True)
        finally:
            session.close()

    def _prepare_session(self) -> requests.Session:
        """Prepare a requests session with retry policy."""
        session = requests.Session()
        retry = Retry(
            total=self.config.retries,
            connect=self.config.retries,
            read=self.config.retries,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=frozenset(["HEAD", "GET", "OPTIONS"])
        )
        adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        session.headers.update({"User-Agent": "dbCAN-DBDownloader/1.0"})
        return session

    def _head_content_length(self, session: requests.Session, url: str) -> Optional[int]:
        """adjust content-length by HEAD request"""
        try:
            resp = session.head(url, timeout=self.config.timeout, allow_redirects=True, verify=self.config.verify_ssl)
            resp.raise_for_status()
            cl = resp.headers.get("content-length")
            return int(cl) if cl is not None else None
        except Exception:
            return None

    def _download_single_file(self, session: requests.Session, url: str, output_path: Path):
        """Download a single file with progress bar (HTTP or S3)

        Args:
            session: requests Session object
            url: The URL to download from (HTTP or S3 URL)
            output_path: The path to save the file to
        """
        tmp_path = Path(str(output_path) + ".part")
        resume_pos = 0
        headers = {}
        if self.config.resume and tmp_path.exists():
            resume_pos = tmp_path.stat().st_size
            if resume_pos > 0:
                headers["Range"] = f"bytes={resume_pos}-"

        with session.get(url, stream=True, timeout=self.config.timeout, headers=headers, verify=self.config.verify_ssl) as response:
            # Check for HTTP errors
            if response.status_code not in (200, 206):
                response.raise_for_status()

            total_size = response.headers.get('content-length')
            # when resuming, total size is the remaining bytes
            total = int(total_size) if total_size is not None else None

            mode = 'ab' if "Range" in headers and response.status_code == 206 else 'wb'
            if mode == 'wb':
                resume_pos = 0  # reset if not resuming

            desc = f"Downloading {output_path.name}"
            with tmp_path.open(mode) as f, tqdm(
                desc=desc,
                total=(resume_pos + (total or 0)) if total is not None else None,
                initial=resume_pos,
                unit='iB',
                unit_scale=True,
                leave=True
            ) as bar:
                for chunk in response.iter_content(chunk_size=1024 * 256):
                    if not chunk:
                        continue
                    f.write(chunk)
                    bar.update(len(chunk))

        # Rename temp file to final output path
        tmp_path.replace(output_path)

    def _extract_tar_file(self, tar_path: Path, destination: Path):
        """Extract a tar archive

        Args:
            tar_path: Path to the tar archive
            destination: Destination directory
        """
        if not tarfile.is_tarfile(str(tar_path)):
            logger.error(f"File is not a valid tar archive: {tar_path}")
            return

        with tarfile.open(str(tar_path)) as tar:
            members = tar.getmembers()

            def _is_within_directory(directory: Path, target: Path) -> bool:
                abs_directory = directory.resolve()
                abs_target = target.resolve()
                return abs_target.is_relative_to(abs_directory) if hasattr(Path, "is_relative_to") else str(abs_target).startswith(str(abs_directory))

            safe_members = []
            for m in members:
                target_path = destination / m.name
                if _is_within_directory(destination, target_path):
                    safe_members.append(m)
                else:
                    logger.warning(f"Skipping suspicious tar member path: {m.name}")

            for member in tqdm(safe_members, desc="Extracting files"):
                tar.extract(member, path=str(destination))

            logger.info(f"Successfully extracted {len(safe_members)} files to {destination}")
            try:
                tar_path.unlink(missing_ok=True)
                logger.info(f"Removed tar file: {tar_path}")
            except Exception as e:
                logger.warning(f"Failed to remove tar file {tar_path}: {e}")


# if __name__ == "__main__":
#     config = DBDownloaderConfig(db_dir="db")
#     downloader = DBDownloader(config)
#     downloader.download_file()
