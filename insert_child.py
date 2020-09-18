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


def ingest(connection, cursor, data, samples, vcffile, dataset_name, dataset_description, full_access=False):
    """
    Given the Postgres connection, a cursor, and the data, insert the
    relevant data into the postgres database.
    """
    now = datetime.datetime.now()

    if full_access:
        access_type = "REGISTERED"
    else:
        access_type = "CONTROLLED"

    cursor.execute("SELECT MAX(id) from dataset_table;")
    last_id = cursor.fetchone()[0]

    cursor.execute("""INSERT INTO dataset_table(id, stable_id, description, access_type,
                                                reference_genome, variant_cnt, call_cnt,
                                                sample_cnt)
                      VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                      RETURNING id""", (last_id+1, f"{dataset_name}_{access_type}",
                                        f"{dataset_description} {access_type} access",
                                        access_type, "GRCh37", 100, 100, len(samples)))
    dataset_id = cursor.fetchone()[0]

    padded_samples = samples + [None]*len(data)
    sample_to_id = {}

    for idnum, (patient, sample) in enumerate(zip(data, padded_samples)):
        demographic = patient["demographic"]
        stable_id = f"Sythetic_CHILD{(idnum+1):03d}"

        sex = demographic["biologicalSex"].lower()
        ethnicity = demographic.get("ethnicity", "")
        if isinstance(ethnicity, list):
            ethnicity = ethnicity[0]
        dob = datetime.datetime.strptime(demographic["age"], "%d-%b-%y")
        age = dateutil.relativedelta.relativedelta(now, dob).years

        geographic_origin = ""
        if "residence" in demographic:
            geographic_origin = demographic("residence")
            if isinstance(geographic_origin, list):
                geographic_origin = geographic_origin[0]

        cursor.execute("SELECT MAX(id) from individual;")
        last_id = cursor.fetchone()[0]
        cursor.execute("""INSERT INTO individual(id, stable_id, sex, ethnicity, geographic_origin)
                          VALUES (%s, %s, %s, %s, %s)
                          RETURNING id""", (last_id + 1, stable_id, sex, ethnicity, geographic_origin))
        patient_id = cursor.fetchone()[0]
        print(patient_id, stable_id, sex, ethnicity, geographic_origin)

        diseases = patient["diseases"]
        for disease, present in diseases.items():
            if isinstance(present, str):
                present = [present]
            has_disease = any([p not in ["No", "Never"] for p in present])

            if has_disease:
                disease_name = disease
                if disease_name == "oncological":
                    disease_name = "Cancer"

                cursor.execute("SELECT MAX(id) from individual_disease_table;")
                last_id = cursor.fetchone()[0]
                cursor.execute("""INSERT INTO individual_disease_table(id, individual_id, disease_id, age)
                               VALUES (%s, %s, %s, %s)""",
                               (last_id+1, patient_id, disease, dob))

        if sample:
            cursor.execute("SELECT MAX(id) from sample_table;")
            last_id = cursor.fetchone()[0]
            cursor.execute("""INSERT INTO sample_table(id, stable_id, individual_id, individual_age_at_collection, collection_date)
                            VALUES (%s, %s, %s, %s, %s)
                            RETURNING id""", (last_id + 1, sample, patient_id, age, now))
            sample_id = cursor.fetchone()[0]
            sample_to_id[sample] = sample_id
            assert last_id + 1 == sample_id
            

            cursor.execute("SELECT MAX(id) from dataset_sample_table;")
            last_id = cursor.fetchone()[0]
            cursor.execute("""INSERT INTO dataset_sample_table(id, dataset_id, sample_id)
                            VALUES (%s, %s, %s)""", (last_id + 1, dataset_id, sample_id))

        connection.commit()

    nsamples = len(samples)
    count, callcount = 0, 0
    cursor.execute("SELECT MAX(id) from variant_table;")
    last_variant_table_id = cursor.fetchone()[0]
    cur_var_id = last_variant_table_id + 1
    if nsamples > 0:
        vcf_reader = vcf.Reader(open(vcffile, 'r'))
        for record in vcf_reader:
            has_var = [call.sample for call in record.samples if call.is_variant]
            if not has_var:
                continue
            startpos = record.POS
            chrom, ref, alt = record.CHROM, record.REF, str(record.ALT[0])
            endpos = startpos + len(ref) - 1

            cursor.execute("""INSERT INTO variant_table(id, dataset_id, variant_id, refseq, reference, alternate, start, "end", call_cnt, sample_cnt, matching_sample_cnt, frequency)
                              VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                              RETURNING id""",
                              (cur_var_id, dataset_id, record.ID, chrom, ref, alt, startpos, endpos, nsamples, nsamples, len(has_var), 1.*len(has_var)/nsamples))
            variant_id = cursor.fetchone()[0]
            cur_var_id = cur_var_id + 1

            samples_w_variant = [(variant_id, sample_to_id[samp]) for samp in has_var if samp in sample_to_id]
            cursor.executemany("""INSERT INTO variant_sample_table(variant_id, sample_id)
                                  VALUES(%s,%s)""", samples_w_variant)
            count += 1
            callcount += len(has_var)
            if count % 100 == 0:
                print(count, callcount)

            if count % 1000 == 0:
                break

        cursor.execute("""INSERT INTO dataset_access_level_table
                          (dataset_id, parent_field, field, access_level)
                          VALUES (%s, 'accessLevelSummary', '-', 'PUBLIC')""",
                          (dataset_id,))
        cursor.execute("""INSERT INTO dataset_consent_code_table
                          (dataset_id, consent_code_id , additional_constraint, version) 
                          VALUES (%s, 1, null, 'v1.0')""",
                          (dataset_id,))
        connection.commit()


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
    parser.add_argument("--server", help="Postgres server hostname", default="localhost")
    parser.add_argument("--port", help="Postgres server port", default=5432)
    parser.add_argument("--username", help="Postgres username", default="beacon")
    parser.add_argument("--database", help="Postgres database", default="beacon")
    parser.add_argument("--password", help="Postgres password", default="secretpassword")
    parser.add_argument("--vcffile", help="VCF file", default=None)
    parser.add_argument("-n", "--dataset_name", default="synthetic_childdemo")
    parser.add_argument("-D", "--dataset_description", default="Synthetic CHILD demo for MTR")
    parser.add_argument("-F", "--full_access", action="store_true")
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

    ingest(connection, cursor, data, samples, args.vcffile,
           args.dataset_name, args.dataset_description,
           full_access=args.full_access)

if __name__ == "__main__":
    main()
