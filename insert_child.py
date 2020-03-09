#!/usr/bin/env python3
"""
This program inserts WP3-mapped CHILDdb data from a provided json file
into the postgres database used by the reference Beacon 2.0 implementation.
"""
import argparse
import datetime
import json
import sys
import dateutil.relativedelta
import psycopg2 as pg
import vcf


def ingest(connection, cursor, data, samples, vcffile):
    """
    Given the Postgres connection, a cursor, and the data, insert the
    relevant data into the postgres database.
    """
    now = datetime.datetime.now()

    cursor.execute("""INSERT INTO beacon_dataset_table(stable_id, description, access_type, reference_genome, variant_cnt, call_cnt, sample_cnt)
                      VALUES (%s, %s, %s, %s, %s, %s, %s)
                      RETURNING id""", ("childdemo", "CHILD db demo for AGM", "PUBLIC", "GRCh37", 100, 100, len(samples)))
    dataset_id = cursor.fetchone()[0]

    padded_samples = samples + [None]*len(data)
    sample_to_id = {}

    for idnum, (patient, sample) in enumerate(zip(data, padded_samples)):
        demographic = patient["demographic"]
        stable_id = f"CHILD{(idnum+1):03d}"
        sex = demographic["biologicalSex"].lower()
        dob = datetime.datetime.strptime(demographic["age"], "%d-%b-%y")
        age_of_onset = dateutil.relativedelta.relativedelta(now, dob).years

        diseases = patient["diseases"]
        disease = ""
        if "bloodRelatedDisorders" in diseases and (not diseases["bloodRelatedDisorders"][0] in ["Never", "No"]):
            disease = "bloodRelatedDisorders"
        elif "respiratorySystem" in diseases and (not diseases["respiratorySystem"][0] in ["Never", "No"]):
            disease = "respiratorySystem"
        elif "circulatorySystem" in diseases and (not diseases["circulatorySystem"][0] in ["Never", "No"]):
            disease = "circulatorySystem"

        cursor.execute("""INSERT INTO patient_table(stable_id, sex, age_of_onset, disease)
                          VALUES (%s, %s, %s, %s)
                          RETURNING id""", (stable_id, sex, age_of_onset, disease))
        patient_id = cursor.fetchone()[0]
        print(patient_id, stable_id, sex, age_of_onset, disease)

        if sample:
            cursor.execute("""INSERT INTO beacon_sample_table(stable_id, sex, tissue, patient_id)
                            VALUES (%s, %s, %s, %s)
                            RETURNING id""", (sample, sex, "blood", patient_id))
            sample_id = cursor.fetchone()[0]
            sample_to_id[sample] = sample_id

            cursor.execute("""INSERT INTO beacon_dataset_sample_table(dataset_id, sample_id)
                            VALUES (%s, %s)""", (dataset_id, sample_id))

        connection.commit()

    nsamples = len(samples)
    count, callcount = 0, 0
    if nsamples > 0:
        vcf_reader = vcf.Reader(open(vcffile, 'r'))
        for record in vcf_reader:
            has_var = [call.sample for call in record.samples if call.is_variant]
            if not has_var:
                continue
            startpos = record.POS
            chrom, ref, alt = record.CHROM, record.REF, str(record.ALT[0])
            endpos = startpos + len(ref) - 1
            cursor.execute("""INSERT INTO beacon_data_table(dataset_id, variant_id, chromosome, reference, alternate, start, "end", call_cnt, sample_cnt, matching_sample_cnt, frequency)
                              VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                              RETURNING id""",
                              (dataset_id, record.ID, chrom, ref, alt, startpos, endpos, nsamples, nsamples, len(has_var), 1.*len(has_var)/nsamples))
            variant_id = cursor.fetchone()[0]

            samples_w_variant = [(variant_id, sample_to_id[samp]) for samp in has_var]
            cursor.executemany("""INSERT INTO beacon_data_sample_table(data_id, sample_id)
                                  VALUES(%s,%s)""", samples_w_variant)
            count += 1
            callcount += len(has_var)
            if count % 100 == 0:
                print(count, callcount)

            if count % 1000 == 0:
                break

        cursor.execute("""UPDATE beacon_dataset_table
                          SET variant_cnt=%s, call_cnt=%s
                          WHERE id=%s""",
                          (count, callcount, dataset_id))
        cursor.execute("""INSERT INTO dataset_access_level_table
                          (dataset_id, parent_field, field, access_level)
                          VALUES (%s, 'accessLevelSummary', '-', 'PUBLIC')""",
                          (dataset_id,))
        cursor.execute("""INSERT INTO beacon_dataset_consent_code_table
                          (dataset_id, consent_code_id , additional_constraint, version) 
                          VALUES (%s, 1, null, 'v1.0')""",
                          (dataset_id,))
        connection.commit()


def clear_tables(connection, cursor):
    """
    Rmove initial data if any from Beacon 2.x tables
    """
    cursor.execute("""TRUNCATE TABLE patient_table, beacon_data_table, beacon_sample_table,
                      beacon_dataset_table, beacon_dataset_sample_table, beacon_data_sample_table,
                      tmp_data_sample_table, tmp_sample_table,
                      beacon_dataset_consent_code_table, dataset_access_level_table""")
    connection.commit()

    return


def vcf_samples(vcffile):
    """
    Samples from vcf file if present
    """
    try:
        vcf_reader = vcf.Reader(open(vcffile, 'r'))
        return vcf_reader.samples
    except:
        print(f"Could not read vcffile {vcffile}: continuing without vcf data")

    return []


def main():
    """
    Parse arguments, make database connection, read file, and start ingest
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--server", help="Postgres server hostname", default="0.0.0.0")
    parser.add_argument("--port", help="Postgres server port", default=5049)
    parser.add_argument("--username", help="Postgres username", default="beacon")
    parser.add_argument("--password", help="Postgres password", default="secretpassword")
    parser.add_argument("--database", help="Postgres database", default="beacon")
    parser.add_argument("--vcffile", help="VCF file", default=None)
    parser.add_argument("datafile", default="./child.json")
    args = parser.parse_args()

    connection = None
    try:
        connection = pg.connect(user=args.username,
                                password=args.password,
                                host=args.server,
                                port=args.port,
                                database=args.database)
        cursor = connection.cursor()
    except (Exception, pg.Error) as error:
        print("Error while connecting to PostgreSQL", error, file=sys.stderr)
        if connection:
            cursor.close()
            connection.close()

    if not connection:
        return

    with open(args.datafile, 'r') as infile:
        data = json.load(infile)

    if not data:
        print(f"Error reading data file {args.datafile}.", file=sys.stderr)
        return

    if args.vcffile:
        samples = vcf_samples(args.vcffile)
    else:
        samples = []

    clear_tables(connection, cursor)
    ingest(connection, cursor, data, samples, args.vcffile)

if __name__ == "__main__":
    main()
