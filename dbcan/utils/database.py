import requests
import os
import logging
import tarfile
from tqdm import tqdm
from dbcan.parameter import DBDownloaderConfig
from dbcan.constants import (CAZY_DB_URL, HMMER_DB_URL, DBCAN_SUB_DB_URL,
                            DBCAN_SUB_MAP_URL, TCDB_DB_URL, TF_DB_URL,
                            STP_DB_URL, PUL_DB_URL, PUL_MAP_URL, PUL_ALL_URL,
                            PEPTIDASE_DB_URL, SULFATLAS_DB_URL)

class DBDownloader:
    """Download dbCAN databases"""

    def __init__(self, config: DBDownloaderConfig):
        """Initialize the database downloader

        Args:
            config: DBDownloaderConfig parameter
        """
        self.config = config
        self.db_dir = config.db_dir
        if not os.path.exists(self.db_dir):
            os.makedirs(self.db_dir)

        # Define path for all database
        self.databases = {
            "CAZy.dmnd": CAZY_DB_URL,
            "dbCAN.hmm": HMMER_DB_URL,
            "dbCAN-sub.hmm": DBCAN_SUB_DB_URL,
            "fam-substrate-mapping.tsv": DBCAN_SUB_MAP_URL,
            "TCDB.dmnd": TCDB_DB_URL,
            "TF.hmm": TF_DB_URL,
            "STP.hmm": STP_DB_URL,
            "PUL.dmnd": PUL_DB_URL,
            "dbCAN-PUL.xlsx": PUL_MAP_URL,
            "dbCAN-PUL.tar.gz": PUL_ALL_URL,
            "peptidase_db.dmnd": PEPTIDASE_DB_URL,
            "sulfatlas_db.dmnd": SULFATLAS_DB_URL
        }

    def download_file(self):
        """Download all databases"""
        for filename, url in self.databases.items():
            output_path = os.path.join(self.db_dir, filename)

            # If file exists and no_overwrite is set, skip download
            if os.path.exists(output_path) and hasattr(self.config, 'no_overwrite') and self.config.no_overwrite:
                logging.info(f"File {filename} already exists, skipping download.")
                continue

            logging.info(f"Downloading {filename} from {url}")
            try:
                self._download_single_file(url, output_path)
                logging.info(f"{filename} successfully downloaded")

                # Extract PUL_ALL if it's the tar.gz file
                if filename == "dbCAN-PUL.tar.gz":
                    self._extract_tar_file(output_path, self.db_dir)

            except Exception as e:
                logging.error(f"{filename} download error: {e}")

    def _download_single_file(self, url, output_path):
        """Download a single file with progress bar

        Args:
            url: The URL to download from
            output_path: The path to save the file to
        """
        with requests.get(url, stream=True) as response:
            response.raise_for_status()
            total_size = int(response.headers.get('content-length', 0))

            with open(output_path, 'wb') as f, tqdm(
                desc=f"Downloading {os.path.basename(output_path)}",
                total=total_size,
                unit='iB',
                unit_scale=True,
                leave=True
            ) as bar:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        bar.update(len(chunk))

    def _extract_tar_file(self, tar_path, destination):
        """Extract a tar archive

        Args:
            tar_path: Path to the tar archive
            destination: Destination directory
        """
        if not tarfile.is_tarfile(tar_path):
            logging.error(f"File is not a valid tar archive: {tar_path}")
            return

        with tarfile.open(tar_path) as tar:
            members = tar.getmembers()

            for member in tqdm(members, desc="Extracting files"):
                tar.extract(member, path=destination)

            logging.info(f"Successfully extracted {len(members)} files to {destination}")
            os.remove(tar_path)
            logging.info(f"Removed tar file: {tar_path}")


# if __name__ == "__main__":
#     # Example usage
#     config = DBDownloaderConfig(db_dir="db")
#     downloader = DBDownloader(config)
#     downloader.download_file()
