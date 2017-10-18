import json
import requests
import sqlite3
import glob
import gzip
import bz2
from urllib import parse
from copy import deepcopy
from pathlib import Path

from vcf2clinvar import *
from vcf2clinvar.common import *
import myvariant

def _next_line(filebuffer):
    try:
        next_line = filebuffer.readline()
    except AttributeError:
        next_line = filebuffer.next()
    try:
        next_line = next_line.decode('utf-8')
        return next_line
    except AttributeError:
        return next_line

def _gennotes_id(chrom, pos, ref, var):
    return 'b37-{}-{}-{}-{}'.format(chrom, pos, ref, var)

def _hgvs_id(chrom, pos, ref, var):
    return parse.unquote(myvariant.format_hgvs(chrom, pos, ref, var))

def create_schema(conn):
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS clinvar
        (chrom int, pos int, id text, ref text, alt text, full_data text)
    ''')
    conn.commit()
    c.execute('''
        CREATE INDEX chrom_match ON clinvar (chrom, pos)
    ''')
    conn.commit()
    c.execute('''
        CREATE UNIQUE INDEX id_match ON clinvar (id)
    ''')
    conn.commit()
    c.close()

def parse_clinvar_data(path, conn):
    c = conn.cursor()
    print('Parsing clinvar data...')
    clinvar_file = gzip.open(path)
    clinvar_cur_line = _next_line(clinvar_file)
    while clinvar_cur_line.startswith('#'):
        clinvar_cur_line = _next_line(clinvar_file)
    clinvar_data = []
    while clinvar_cur_line:
        clinvar_vcf_line = ClinVarVCFLine(vcf_line=clinvar_cur_line)
        clinvar_data.append(clinvar_vcf_line)
        clinvar_cur_line = _next_line(clinvar_file)
        c.execute('''
            INSERT INTO clinvar VALUES (?, ?, ?, ?, ?, ?)
        ''', (CHROM_INDEX[clinvar_vcf_line.chrom],
                clinvar_vcf_line.start,
                clinvar_vcf_line.dbsnp_id if clinvar_vcf_line.dbsnp_id != '.' else None,
                clinvar_vcf_line.ref_allele,
                str(clinvar_vcf_line.alt_alleles),
                json.dumps(clinvar_vcf_line.as_dict())))
        conn.commit()
    print('Done!')
    c.close()
    return clinvar_data

def parse_user_vcf_data(user):
    print('Parsing ' + user['user']['username'] + ' VCF data...')
    vcf_file = bz2.BZ2File('user_vcf_data/' + user['local_filename'] + '.vcf.bz2', 'r')
    vcf_cur_line = _next_line(vcf_file)
    while vcf_cur_line.startswith('#'):
        vcf_cur_line = _next_line(vcf_file)
    vcf_data = []
    while vcf_cur_line:
        fields = vcf_cur_line.strip().split('\t')
        vcf_data.append({
            'chrom': fields[0],
            'pos': fields[1],
            'id': fields[2],
            'ref_allele': fields[3],
            'alt_allele': fields[4],
            'qual': fields[5],
            'filter': fields[6],
            'info': fields[7],
            'format': fields[8],
            '23andme_data': fields[9]
        })
        vcf_cur_line = _next_line(vcf_file)
    print('Done!')
    return vcf_data

def download_23andme_vcf():
    vcf_results = []
    dl_url = 'https://www.openhumans.org/api/public-data/?source=twenty_three_and_me'
    meta = {
        'next': ''
    }
    print('Fetching OpenHumans 23andMe users VCF metadata...')
    while dl_url != None:
        meta = requests.get(dl_url).json()
        for result in meta['results']:
            if 'vcf' in result['metadata']['tags']:
                vcf_results.append(result)
        dl_url = meta['next']
    print('Done! ' + str(len(vcf_results)) + ' results fetched.')
    print('Downloading users VCF data...')
    num_dl = 0
    for result in vcf_results:
        data = requests.get(result['download_url'], stream=True)
        print('Downloading ' + result['user']['username'] + ' VCF data...')
        filename = result['user']['username'] + '_' + result['user']['id'] + \
                    '_' + result['created'] + '_23andme_data'
        if Path('user_vcf_data/' + filename + '.vcf.bz2').is_file():
            print('Already downloaded, skipping...')
        else:
            with open('user_vcf_data/' + filename + '.vcf.bz2', 'wb') as f:
                for chunk in data.iter_content(128):
                    f.write(chunk)
            print('Done!')
        print('Adding to user metadata...')
        result['local_filename'] = filename
        print('Done!')
        num_dl += 1
    print('Done! Downloaded ' + str(num_dl) + ' users VCF data.')
    return vcf_results

def map_23andme_clinvar(user_data, conn):
    print('Mapping user 23andMe data with ClinVar...')
    c = conn.cursor()
    mv = myvariant.MyVariantInfo()
    # Cache MyVariantInfo requests
    mv.set_caching('./myvariant_cache', verbose=False)
    # Can definitely be optimized to reduce database or HTTP requests.
    for user in user_data:
        parsed = parse_user_vcf_data(user)
        mapped = []
        print('Mapping ' + user['user']['username'] + ' data to ClinVar...')
        print(len(parsed))
        for var in parsed:
            c.execute('''
                SELECT * FROM clinvar WHERE (chrom=? AND pos=? AND alt LIKE ?) OR (id=?)
            ''', (CHROM_INDEX[str(var['chrom'])], int(var['pos']), var['alt_allele'], var['id']))
            for result in c.fetchall():
                var['clinvar_data'] = json.loads(result[-1])
                var['gennotes_id'] = None
                var['gennotes_data'] = None
                var['hgvs_id'] = None
                var['mv_data'] = None
                if (not var['alt_allele'] == '.'):
                    var['gennotes_id'] = _gennotes_id(var['chrom'], var['pos'], var['ref_allele'], var['alt_allele'])
                    results = requests.get('https://gennotes.herokuapp.com/api/variant/', params={'variant_list': json.dumps([var['gennotes_id']])})
                    var['gennotes_data'] = results.json()
                    var['hgvs_id'] = _hgvs_id(var['chrom'], var['pos'], var['ref_allele'], var['alt_allele'])
                    try:
                        mv_data = mv.getvariant(var['hgvs_id'], fields=['clinvar', 'dbsnp', 'exac'], verbose=False)
                        var['mv_data'] = mv_data
                    except Exception as e:
                        print(var['alt_allele'])
                        print(e)
                mapped.append(var)
        with open('mapped_user_vcf/' + user['local_filename'] + '.json', 'w') as f:
            json.dump(mapped, f, indent=4)
    c.close()

def save_final_user_data(user_data):
    print('Writing user metadata to file...')
    with open('openhumans_23andme_metadata.json', 'w') as f:
        json.dump(user_data, f, indent=4)
    print('Done!')

if __name__ == '__main__':
    conn = sqlite3.connect('./clinvar.db')
    if not Path('./clinvar.db').is_file():
        create_schema(conn)
    user_data = download_23andme_vcf()
    if not Path('./clinvar.db').is_file():
        parse_clinvar_data('./clinvar.vcf.gz', conn)
    map_23andme_clinvar(user_data, conn)
    save_final_user_data(user_data)
    conn.close()

