import os
import shutil

path = 'docs/notebooks/'
test_path = 'tests/notebooks/'


def run_notebook(prefix):
    # if exists, save and overwrite
    nb_filename = prefix + '.ipynb'  # ''
    config_filename = prefix + '.config.json'  # By default, benchmark config, but overridden for tests
    config_filename_tmp = prefix + '.config.tmp.json'

    # Save content in memory in tmp file
    existing_config_to_restore = os.path.exists(path + config_filename)
    if existing_config_to_restore:
        shutil.copy(path + config_filename, path + config_filename_tmp)

    # Overwritting
    shutil.copy(test_path + config_filename, path + config_filename)
    os.chdir('./docs/notebooks')
    try:
        os.system('jupyter nbconvert --execute {:s} --to html --ExecutePreprocessor.kernel_name=python'.format(
            nb_filename)
        )
    except BaseException:
        pass
    os.chdir('../..')

    if existing_config_to_restore:
        shutil.move(path + config_filename_tmp, path + config_filename)


def test_notebooks():
    for prefix in ['annotation', 'scRNA_and_smFISH', 'data_loading', 'basic_tutorial', 'scVI_reproducibility']:
        print(prefix)
        run_notebook(prefix)
