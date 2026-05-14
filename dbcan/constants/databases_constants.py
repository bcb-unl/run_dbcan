#############constants for downloading dbCAN databases##############################

HTTP_DATABASE_SNAPSHOT = "db_current"
S3_DATABASE_RELEASE = "db_v5-2-9_5-5-2026"

BASE_URL = f"https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/{HTTP_DATABASE_SNAPSHOT}"
AWS_S3_URL = f"https://dbcan.s3.us-west-2.amazonaws.com/{S3_DATABASE_RELEASE}"


# main dbCAN database
CAZY_DB_URL = f"{BASE_URL}/CAZy.dmnd"
HMMER_DB_URL = f"{BASE_URL}/dbCAN.hmm"
DBCAN_SUB_DB_URL = f"{BASE_URL}/dbCAN_sub.hmm"
DBCAN_SUB_MAP_URL = f"{BASE_URL}/fam-substrate-mapping.tsv"

## AWS S3 URLs
CAZY_DB_S3_URL = f"{AWS_S3_URL}/CAZy.dmnd"
HMMER_DB_S3_URL = f"{AWS_S3_URL}/dbCAN.hmm"
DBCAN_SUB_DB_S3_URL = f"{AWS_S3_URL}/dbCAN_sub.hmm"
DBCAN_SUB_MAP_S3_URL = f"{AWS_S3_URL}/fam-substrate-mapping.tsv"

# TF/ STP/ TCDB/ Peptidase/ SulfAtlas for functional annotation
TCDB_DB_URL = f"{BASE_URL}/TCDB.dmnd"
TF_DB_URL = f"{BASE_URL}/TF.hmm"
STP_DB_URL = f"{BASE_URL}/STP.hmm"
TF_DIAMOND_DB_URL = f"{BASE_URL}/TF.dmnd"
PEPTIDASE_DB_URL = f"{BASE_URL}/peptidase_db.dmnd"
SULFATLAS_DB_URL = f"{BASE_URL}/sulfatlas_db.dmnd"

## AWS S3 URLs
TCDB_DB_S3_URL = f"{AWS_S3_URL}/TCDB.dmnd"
TF_DB_S3_URL = f"{AWS_S3_URL}/TF.hmm"
STP_DB_S3_URL = f"{AWS_S3_URL}/STP.hmm"
TF_DIAMOND_DB_S3_URL = f"{AWS_S3_URL}/TF.dmnd"
PEPTIDASE_DB_S3_URL = f"{AWS_S3_URL}/peptidase_db.dmnd"
SULFATLAS_DB_S3_URL = f"{AWS_S3_URL}/sulfatlas_db.dmnd"

# PUL
PUL_DB_URL = f"{BASE_URL}/PUL.dmnd"
PUL_MAP_URL = f"{BASE_URL}/dbCAN-PUL.xlsx"
PUL_ALL_URL = f"{BASE_URL}/dbCAN-PUL.tar.gz"

##AWS S3 URLs
PUL_DB_S3_URL = f"{AWS_S3_URL}/PUL.dmnd"
PUL_MAP_S3_URL = f"{AWS_S3_URL}/dbCAN-PUL.xlsx"
PUL_ALL_S3_URL = f"{AWS_S3_URL}/dbCAN-PUL.tar.gz"

# Pfam
#PFAM_DB_URL = f"{BASE_URL}/Pfam-A.hmm"
#PFAM_DB_S3_URL = f"{AWS_S3_URL}/Pfam-A.hmm"

#############dictionary of databases##############################
# all database names and urls
DATABASES_CAZYME = {
            "CAZy.dmnd": CAZY_DB_URL,
            "dbCAN.hmm": HMMER_DB_URL,
            "dbCAN-sub.hmm": DBCAN_SUB_DB_URL,
            "fam-substrate-mapping.tsv": DBCAN_SUB_MAP_URL,
        }

DATABASES_CAZYME_S3 = {
            "CAZy.dmnd": CAZY_DB_S3_URL,
            "dbCAN.hmm": HMMER_DB_S3_URL,
            "dbCAN-sub.hmm": DBCAN_SUB_DB_S3_URL,
            "fam-substrate-mapping.tsv": DBCAN_SUB_MAP_S3_URL,
        }

DATABASES_CGC = {
            "TCDB.dmnd": TCDB_DB_URL,
            "TF.hmm": TF_DB_URL,
            "TF.dmnd": TF_DIAMOND_DB_URL,
            "STP.hmm": STP_DB_URL,
            "PUL.dmnd": PUL_DB_URL,
            "dbCAN-PUL.xlsx": PUL_MAP_URL,
            "dbCAN-PUL.tar.gz": PUL_ALL_URL,
            "peptidase_db.dmnd": PEPTIDASE_DB_URL,
            "sulfatlas_db.dmnd": SULFATLAS_DB_URL,
        }

DATABASES_CGC_S3 = {
            "TCDB.dmnd": TCDB_DB_S3_URL,
            "TF.hmm": TF_DB_S3_URL,
            "TF.dmnd": TF_DIAMOND_DB_S3_URL,
            "STP.hmm": STP_DB_S3_URL,
            "PUL.dmnd": PUL_DB_S3_URL,
            "dbCAN-PUL.xlsx": PUL_MAP_S3_URL,
            "dbCAN-PUL.tar.gz": PUL_ALL_S3_URL,
            "peptidase_db.dmnd": PEPTIDASE_DB_S3_URL,
            "sulfatlas_db.dmnd": SULFATLAS_DB_S3_URL,
        }

COMPRESSED_DBCAN_PUL = "dbCAN-PUL.tar.gz"



