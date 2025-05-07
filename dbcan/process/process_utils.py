import pandas as pd
from dbcan.constants import HMMER_COLUMN_NAMES, OVERLAP_RATIO_THRESHOLD
import os
import logging

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

def process_cgc_sig_results(tc_config, tf_config, stp_config, sulfatase_config, peptidase_config):
    """combine TCDB, TF and STP results into one file"""
    try:
        columns = ['Annotate Name', 'Annotate Length', 'Target Name', 'Target Length',
                'i-Evalue', 'Annotate From', 'Annotate To', 'Target From', 'Target To',
                'Coverage', 'Annotate File Name']

        # Get output directory from any of the config objects
        output_dir = getattr(tc_config, 'output_dir', None)
        if not output_dir:
            output_dir = getattr(tf_config, 'output_dir', None)
        if not output_dir:
            output_dir = getattr(stp_config, 'output_dir', '.')

        # Define standard file paths based on output_dir
        output_files = {
            'TC': os.path.join(output_dir, 'diamond.out.tc'),
            'TF': os.path.join(output_dir, 'TF_hmm_results.tsv'),
            'STP': os.path.join(output_dir, 'STP_hmm_results.tsv'),
            'Sulfatase': os.path.join(output_dir, 'diamond.out.sulfatlas'),
            'Peptidase': os.path.join(output_dir, 'diamond.out.peptidase')

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
