import pandas as pd
from dbcan.constants import HMMER_COLUMN_NAMES, OVERLAP_RATIO_THRESHOLD, CGC_RESULT_FILE
import os
import logging
from Bio import SeqIO


def process_results(results, output_file):
    if results:
        df = pd.DataFrame(results, columns=HMMER_COLUMN_NAMES)
        df.sort_values(by=['Target Name', 'Target From', 'Target To'], inplace=True)
        df_filtered = filter_overlaps(df)
        df_filtered.to_csv(output_file, index=False, sep='\t')
    else:
        df = pd.DataFrame(columns=HMMER_COLUMN_NAMES)
        df.to_csv(output_file, index=False, sep='\t')

def filter_overlaps(df):
    filtered = []
    grouped = df.groupby('Target Name')

    for name, group in grouped:
        group = group.reset_index(drop=True)
        keep = []

        for i in range(len(group)):
            if not keep:
                keep.append(group.iloc[i])
                continue

            last = keep[-1]
            current = group.iloc[i]
            overlap = min(last['Target To'], current['Target To']) - max(last['Target From'], current['Target From'])
            if overlap > 0:
                overlap_ratio_last = overlap / (last['Target To'] - last['Target From'])
                overlap_ratio_current = overlap / (current['Target To'] - current['Target From'])

                if overlap_ratio_last > OVERLAP_RATIO_THRESHOLD or overlap_ratio_current > OVERLAP_RATIO_THRESHOLD:
                    if last['i-Evalue'] > current['i-Evalue']:
                        keep[-1] = current
                else:
                    keep.append(current)
            else:
                keep.append(current)

        filtered.extend(keep)

    return pd.DataFrame(filtered)

def process_cgc_sig_results(tc_config, tfdiamond_config, tf_config, stp_config, sulfatase_config, peptidase_config):
    """combine TCDB, TF and STP results into one file"""
    try:
        columns = ['Annotate Name', 'Annotate Length', 'Target Name', 'Target Length',
                'i-Evalue', 'Annotate From', 'Annotate To', 'Target From', 'Target To',
                'Coverage', 'Annotate File Name']

        # Get output directory from any of the config objects
        output_dir = getattr(tc_config, 'output_dir', None)
        if not output_dir:
            output_dir = getattr(tfdiamond_config, 'output_dir', None)
        if not output_dir:
            output_dir = getattr(tf_config, 'output_dir', None)
        if not output_dir:
            output_dir = getattr(stp_config, 'output_dir', '.')

        # Define standard file paths based on output_dir
        output_files = {
            'TC': os.path.join(output_dir, 'diamond.out.tc'),
            'TF': os.path.join(output_dir, 'diamond.out.tf'),
            'TF': os.path.join(output_dir, 'TF_hmm_results.tsv'),
            'STP': os.path.join(output_dir, 'STP_hmm_results.tsv'),
            'Sulfatase': os.path.join(output_dir, 'diamond.out.sulfatlas'),
            'Peptidase': os.path.join(output_dir, 'diamond.out.peptidase'),
        }

        # check files exist and are not empty
        dataframes = []

        for name, file_path in output_files.items():
            if not file_path or not os.path.exists(file_path):
                logging.warning(f"{name} output file not found: {file_path}")
                continue

            if os.path.getsize(file_path) == 0:
                logging.warning(f"{name} output file is empty: {file_path}")
                continue

            try:
                df = pd.read_csv(file_path, names=columns, header=0, sep='\t')
                df['Type'] = name  # add a new column to identify the type
                dataframes.append(df)
                logging.info(f"Loaded {len(df)} {name} annotations from {file_path}")
            except Exception as e:
                logging.error(f"Error reading {name} output file: {e}")

        if not dataframes:
            logging.warning("No valid CGC annotation data found")
            # generate empty file
            output_file = os.path.join(output_dir, 'total_cgc_info.tsv')
            pd.DataFrame(columns=columns + ['Type']).to_csv(output_file, index=False, sep='\t')
            return

        # combine dataframes
        total_function_annotation_df = pd.concat(dataframes, ignore_index=True)

        # filter overlaps
        filtered_df = filter_overlaps(total_function_annotation_df)

        # save combined and filtered results
        output_file = os.path.join(output_dir, 'total_cgc_info.tsv')
        filtered_df.to_csv(output_file, index=False, sep='\t')
        logging.info(f"Saved {len(filtered_df)} CGC annotations to {output_file}")

    except Exception as e:
        logging.error(f"Error processing CGC signature results: {e}")
        import traceback
        traceback.print_exc()

def process_cgc_null_pfam_annotation(Pfam_config):
    """process CGC null pfam annotation results"""
    try:
        output_dir = getattr(Pfam_config, 'output_dir', '.')
        pfam_hmm_output = os.path.join(output_dir, 'pfam_hmm_results.tsv')

        if not os.path.exists(pfam_hmm_output) or os.path.getsize(pfam_hmm_output) == 0:
            logging.warning(f"PFAM HMM output file not found or empty: {pfam_hmm_output}")
            # generate empty file
            pd.DataFrame(columns=HMMER_COLUMN_NAMES).to_csv(pfam_hmm_output, index=False, sep='\t')
            return

        df = pd.read_csv(pfam_hmm_output, names=HMMER_COLUMN_NAMES, header=0, sep='\t')
        df.sort_values(by=['Target Name', 'Target From', 'Target To'], inplace=True)
        df_filtered = filter_overlaps(df)
        df_filtered.to_csv(pfam_hmm_output, index=False, sep='\t')
        logging.info(f"Processed PFAM HMM results saved to: {pfam_hmm_output}")
    except Exception as e:
        logging.error(f"Error processing CGC null PFAM annotation: {e}")
        import traceback
        traceback.print_exc()

def extract_null_protein_ids_from_cgc(cgc_standard_out_file):
    null_protein_ids = set()
    with open(cgc_standard_out_file) as f:
        for line in f:
            if line.startswith("CGC#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 4 and fields[1].lower() == "null":
                null_protein_ids.add(fields[3])
    return null_protein_ids

def extract_fasta_by_protein_ids(input_faa, output_faa, protein_ids):
    count = 0
    with open(output_faa, "w") as out_handle:
        for record in SeqIO.parse(input_faa, "fasta"):
            if record.id in protein_ids:
                SeqIO.write(record, out_handle, "fasta")
                count += 1
    return count
def extract_null_fasta_from_cgc(cgc_standard_out_file, input_faa, output_faa):
    null_protein_ids = extract_null_protein_ids_from_cgc(cgc_standard_out_file)
    count = extract_fasta_by_protein_ids(input_faa, output_faa, null_protein_ids)
    logging.info(f"Extracted {count} null protein sequences to {output_faa}")

def annotate_cgc_null_with_pfam_and_gff(cgc_standard_out_file, pfam_hmm_result_file, gff_file, output_cgc_file, output_gff_file):
    """
    Annotate CGC null entries with Pfam annotations in both cgc_standard_out.tsv and cgc.gff files.
    """
    # 1.read pfam_hmm_result_file to build a mapping
    pfam_map = {}
    pfam_df = pd.read_csv(pfam_hmm_result_file, sep='\t', header=0)
    for _, row in pfam_df.iterrows():
        protein_id = str(row['Target Name'])
        annotation = str(row['Annotate Name'])
        pfam_map[protein_id] = annotation

    # 2. parse cgc_standard_out.tsv
    with open(cgc_standard_out_file) as fin, open(output_cgc_file, 'w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            # type is null and has pfam annotation
            if len(fields) >= 8 and fields[1].lower() == "null":
                protein_id = fields[3]
                if protein_id in pfam_map:
                    fields[1] = "Pfam"
                    fields[7] = pfam_map[protein_id]  # Gene Annotation列
                fout.write('\t'.join(fields) + '\n')
                continue
            fields = line.rstrip('\n').split('\t')
            # type is null and has pfam annotation
            if len(fields) >= 8 and fields[1].lower() == "null":
                protein_id = fields[3]
                if protein_id in pfam_map:
                    fields[1] = "Pfam"
                    fields[7] = pfam_map[protein_id]  # Gene Annotation列
            fout.write('\t'.join(fields) + '\n')

    # 3. process cgc.gff
    with open(gff_file) as fin, open(output_gff_file, 'w') as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                fout.write(line)
                continue
            attr = fields[8]
            attrs = attr.split(';')
            attr_dict = {}
            for item in attrs:
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr_dict[k] = v
            protein_id = attr_dict.get('protein_id', None)
            if protein_id and attr_dict.get('CGC_annotation', '').lower() == 'null' and protein_id in pfam_map:
                attr_dict['CGC_annotation'] = f"Pfam|{pfam_map[protein_id]}"
            new_attr = ';'.join([f"{k}={v}" for k, v in attr_dict.items()])
            fields[8] = new_attr
            fout.write('\t'.join(fields) + '\n')
