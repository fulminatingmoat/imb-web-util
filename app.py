import hashlib
import os

import pandas as pd
from flask import Flask, render_template, request, url_for, jsonify, Response, send_from_directory

from ext import make_celery
from lib import dish_impute, maternal_fetal_impute, call_r_script

# Uploaded files will be stored under static in the current directory

CHUNK_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'chunks')
COMPLETED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'completed')
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'results')
COMPRESSED_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static', 'compressed')

app = Flask(__name__)

# Task queue messaging queue and result storage is done with Celery, backed by sqlite, separate redis server not needed
app.config.update(
    broker_url='sqla+sqlite:///celerybroker.sqlite',
    result_backend='db+sqlite:///celeryresult.sqlite',
)

celery = make_celery(app)


# Prevent circular import by decorating the function here in a wrapper
@celery.task(bind=True)
def gwas2hla_task(*args, **kwargs):
    return dish_impute(*args, **kwargs)


@celery.task(bind=True)
def maternal_fetal_task(*args, **kwargs):
    return maternal_fetal_impute(*args, **kwargs)


@celery.task(bind=True)
def call_r_task(*args, **kwargs):
    return call_r_script(*args, **kwargs)


# Convert uploaded file to pickle, to reduce file size
def convert_to_pickle(file_id, file_name):
    extension = file_name.split('.')[-1]
    gz = False
    if extension == 'gz':
        extension = file_name.split('.')[-2]
        gz = True
    filename = f"{file_id}.{extension}"
    filepath = os.path.join(COMPLETED_DIR, filename + '.gz' if gz else '')
    pickle_filepath = os.path.join(COMPRESSED_DIR, f"{file_id}.pkl")
    if os.path.exists(filepath) and not os.path.exists(pickle_filepath):
        match extension:
            case 'csv':
                pd.read_csv(filepath, sep=',').to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))
            case 'tsv':
                pd.read_csv(filepath, sep='\t').to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))
            case 'txt':
                pd.read_csv(filepath, sep=r'\s|\t|,').to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))
            case 'feather':
                pd.read_feather(filepath).to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))
            case 'parquet' | 'pq':
                pd.read_parquet(filepath).to_pickle(os.path.join(COMPRESSED_DIR, f"{file_id}.pkl"))


# serves the dish utility webpage at https://example.com/
@app.route('/')
def dish_util():
    return render_template('index.html')


# serves the maternal fetal utility webpage at https://example.com/maternal_fetal
@app.route('/maternal_fetal')
def maternal_fetal():
    return render_template('maternal_fetal.html')


# Return if we have a chunk for the given file
def validate_chunk(flow_dict):
    chunk_num = int(flow_dict.get('flowChunkNumber'))
    chunk_size = int(flow_dict.get('flowChunkSize'))
    total_size = int(flow_dict.get('flowTotalSize'))
    file_name = flow_dict.get('flowFilename')
    file_id = flow_dict.get('flowIdentifier')
    current_chunk_size = int(flow_dict.get('flowCurrentChunkSize'))
    total_chunks = int(flow_dict.get('flowTotalChunks'))
    upload_token = flow_dict.get('flowUploadToken')
    if chunk_num > total_chunks:
        return True
    # f"{file_id}.{file_name.split('.')[-1]}"
    chunk = os.path.join(CHUNK_DIR, file_id, f"{file_id}.{chunk_num}")
    completed_file = os.path.join(COMPLETED_DIR, "{file_id}.{file_name.split('.')[-1]}")
    compressed_file = os.path.join(COMPRESSED_DIR, f"{file_id}.pkl")
    # If pickled file exists, we have already processed this file
    if os.path.exists(compressed_file):
        if os.path.exists(completed_file):
            os.remove(completed_file)
        return True
    try:
        # If chunk file doesn't exist and completed file isn't there, we need the chunk
        if f"{file_id}.{chunk_num}" not in os.listdir(os.path.join(CHUNK_DIR, file_id)):
            if completed_file not in os.listdir(COMPLETED_DIR):
                return False
            # If the completed file is the wrong size, need all chunks to rebuild
            if os.path.getsize(completed_file) != total_size:
                return False
    except FileNotFoundError:
        return False
    # If chunk is wrong size, need the chunk again
    if os.path.getsize(chunk) != current_chunk_size:
        return False
    return True


# Server side component of the upload flow, if get request, returns if the current chunk is valid
# post request, saves the chunk to disk if chunk isn't valid
@app.route('/upload', methods=['GET', 'POST'])
def results():
    if request.method == 'POST':
        file = request.files['file']
        # If we need the chunk, save it to disk
        if not validate_chunk(request.form):
            chunk_name = f"{request.form['flowIdentifier']}.{request.form['flowChunkNumber']}"
            # Make sure the folder for the chunks exists
            try:
                os.mkdir(os.path.join(CHUNK_DIR, request.form['flowIdentifier']))
            except FileExistsError:
                pass
            filename = os.path.join(CHUNK_DIR, request.form['flowIdentifier'], chunk_name)
            file.save(filename)
            if chunk_name in os.listdir(os.path.join(CHUNK_DIR, request.form['flowIdentifier'])):
                return Response(status=200)
            else:
                return Response(status=203)
        else:
            return Response(status=203)
    else:
        # Get requests query if we need the file
        if validate_chunk(request.args):
            return Response(status=200)
        else:
            return Response(status=204)


# Sent at the end of the upload, we can now process the file by combining the chunks
@app.route('/finished', methods=['POST'])
def combine():
    json = request.get_json()
    filename = f"{json['hash']}.{json['file'].split('.')[-1]}"
    with open(os.path.join(COMPLETED_DIR, filename), "wb") as f:
        for chunk in range(1, len(os.listdir(os.path.join(CHUNK_DIR, json['hash']))) + 1):
            with open(os.path.join(CHUNK_DIR, json['hash'], f"{json['hash']}.{chunk}"), "rb") as c:
                f.write(c.read())
    # Check that the file has the expected hash
    md5_hash = hashlib.md5()
    with open(os.path.join(COMPLETED_DIR, filename), "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    if md5_hash.hexdigest() == json['hash']:
        # Save the file as a pickled file if the hash is correct
        convert_to_pickle(md5_hash.hexdigest(), filename)
        return Response(status=200)
    else:
        return Response(status=529)


# If the file is already saved, prevent upload from client
@app.route('/check/<file_hash>', methods=['GET'])
def check_file_exists(file_hash):
    if f'{file_hash}.pkl' in os.listdir(COMPRESSED_DIR):
        return Response(status=200)
    if file_hash in list(map(lambda x: x.split('.')[0], os.listdir(COMPLETED_DIR))):
        return Response(status=200)
    else:
        return Response(status=204)


# Run the analysis, input files referenced by hash
@app.route('/analyze', methods=['POST'])
def analyze():
    json = request.get_json()
    match json['function']:
        case 'GWAS2HLA':
            result = gwas2hla_task.delay(reference_file=json['reference_hash'], map_file=json['reference_map_hash'],
                                         summary_file=json['summary_hash'],
                                         genome_build=json['genome_build'], custom=json['custom'], maf=json['maf'])
            return jsonify({'result': url_for('task_status', task_id=result.id)}), 202, {
                'Location': url_for('task_status', task_id=result.id)}
        case 'MATERNALFETAL':
            result = maternal_fetal_task.delay(maternal_file=json['maternal_hash'],
                                               fetal_file=json['fetal_hash'], intercept=json['LD_score_intercept'])
            return jsonify({'result': url_for('task_status', task_id=result.id)}), 202, {
                'Location': url_for('task_status', task_id=result.id)}


# Check the status of the analysis
@app.route('/status/<task_id>')
def task_status(task_id):
    task = gwas2hla_task.AsyncResult(task_id)
    if task.state == 'PENDING':
        # job did not start yet
        response = {
            'state': task.state,
            'current': 0,
            'total': 1,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'progress': task.info.get('progress', 0),
            'total': task.info.get('total', 1),
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            response['result'] = task.info['result']
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'current': 1,
            'total': 1,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)


# Retrieve result file from hash
@app.route('/result/<result_hash>')
def serve_result(result_hash):
    return send_from_directory(os.path.join(RESULT_DIR), f"{result_hash}.tsv", as_attachment=True)


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=False)
