import hashlib
import os
import shutil

from subprocess import Popen, PIPE, STDOUT

SCRIPT_DIR = COMPRESSED_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts')
COMPRESSED_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'compressed')
TASK_RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'task_results')
RESULT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static', 'results')

default_maternal = '7b59a43b06a6a3d507321dd05a4b5bc4'
default_fetal = '10fe5ad6852fd0fb1f40c77717f33151'


# Example function to call R script from util using passed arguments

def impute(self, maternal_file=None, fetal_file=None, intercept='0.1352'):
    # Read in all the data
    self.update_state(state='PROGRESS', meta={'progress': 0.1, 'total': 1, 'status': 'Loading maternal file'})
    command = ['/usr/bin/Rscript', f'{SCRIPT_DIR}/maternal_fetal.R',
               os.path.join(COMPRESSED_DIR, f"{maternal_file}.pkl"), os.path.join(COMPRESSED_DIR, f"{fetal_file}.pkl"),
               str(intercept), os.path.join(TASK_RESULT_DIR, f'{self.request.id}.tsv')]
    print(command)
    process = Popen(command, stdout=PIPE, stderr=STDOUT)

    print(process.stdout.read().decode('utf-8'))
    # Wait for the process to finish

    print(process.wait())
    self.update_state(state='PROGRESS', meta={'progress': 0.6, 'total': 1, 'status': 'Imputing'})

    filename = f"{self.request.id}.tsv"

    md5_hash = hashlib.md5()
    with open(os.path.join(TASK_RESULT_DIR, filename), "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)

    if not os.path.exists(os.path.join(RESULT_DIR, f"{md5_hash.hexdigest()}.tsv")):
        shutil.move(os.path.join(TASK_RESULT_DIR, filename), os.path.join(RESULT_DIR, f"{md5_hash.hexdigest()}.tsv"))
    else:
        os.remove(os.path.join(TASK_RESULT_DIR, filename))

    self.update_state(state='SUCCESS',
                      meta={'progress': 1, 'total': 1, 'status': 'Done', 'result': md5_hash.hexdigest()})

    return {'filename': f"{filename}.tsv", 'md5_hash': md5_hash.hexdigest(), 'progress': 1, 'total': 1,
            'status': 'Done', 'result': md5_hash.hexdigest()}
