"""
Predicts solubilities based on input smiles, returning a csv with the predicted values
"""
#import packages
from multiprocessing import Pool
import os
import shutil

from typing import List, Tuple


from tqdm import tqdm
import chemprop
from chemprop.data import get_smiles
from chemprop.features import get_features_generator, load_features, save_features
from chemprop.utils import makedirs




#define functions: featurizes data in same manner as done for the model
def load_temp(temp_dir: str) -> Tuple[List[List[float]], int]:
    """
    Loads all features saved as .npz files in load_dir.
    Assumes temporary files are named in order 0.npz, 1.npz, ...
    :param temp_dir: Directory in which temporary .npz files containing features are stored.
    :return: A tuple with a list of molecule features, each molecule's features a list of floats,
    and the number of temporary files.
    """
    features = []
    temp_num = 0
    temp_path = os.path.join(temp_dir, f'{temp_num}.npz')

    while os.path.exists(temp_path):
        features.extend(load_features(temp_path))
        temp_num += 1
        temp_path = os.path.join(temp_dir, f'{temp_num}.npz')

    return features, temp_num


def generate_and_save_features(data_path: str, save_path: str, smiles_column: str=None, features_generator: str = 'rdkit_2d_normalized', save_frequency: int = 10000, restart: bool = True, sequential: bool = False ):
    """
    Computes and saves features for a dataset of molecules as a 2D array in a .npz file.
    :param args: Arguments.
    """
    # Create directory for save_path
    makedirs(save_path, isfile=True)

    # Get data and features function
    smiles = get_smiles(path=data_path, smiles_columns=smiles_column, flatten=True)
    features_generator = get_features_generator(features_generator)
    temp_save_dir = save_path + '_temp'

    # Load partially complete data
    if restart:
        if os.path.exists(save_path):
            os.remove(save_path)
        if os.path.exists(temp_save_dir):
            shutil.rmtree(temp_save_dir)
    else:
        if os.path.exists(save_path):
            raise ValueError(f'"{save_path}" already exists and args.restart is False.')

        if os.path.exists(temp_save_dir):
            features, temp_num = load_temp(temp_save_dir)

    if not os.path.exists(temp_save_dir):
        makedirs(temp_save_dir)
        features, temp_num = [], 0

    # Build features map function
    smiles = smiles[len(features):] # restrict to data for which features have not been computed yet

    if sequential:
        features_map = map(features_generator, smiles)
    else:
        features_map = Pool().imap(features_generator, smiles)

    # Get features
    temp_features = []
    for i, feats in tqdm(enumerate(features_map), total=len(smiles)):
        temp_features.append(feats)

        # Save temporary features every save_frequency
        if (i > 0 and (i + 1) % save_frequency == 0) or i == len(smiles) - 1:
            save_features(os.path.join(temp_save_dir, f'{temp_num}.npz'), temp_features)
            features.extend(temp_features)
            temp_features = []
            temp_num += 1

    try:
        # Save all features
        save_features(save_path, features)

        # Remove temporary features
        shutil.rmtree(temp_save_dir)
    except OverflowError:
        print('Features array is too large to save as a single file. Instead keeping features as a directory of files.')

# Featurize input data
#Import and save featurised data
def app_exe(config):
    """
    Predicts values for input dataframe
    """
    generate_and_save_features(data_path=config["input_files"], save_path='pred.npz', features_generator = 'rdkit_2d_normalized')

    arguments = [
        '--test_path', config["input_files"],
        '--preds_path', 'predicted_values.csv',
        '--features_path', 'pred.npz',
        '--no_features_scaling',
        '--checkpoint_dir', 'solubility_model'
    ]

    args = chemprop.args.PredictArgs().parse_args(arguments)
    chemprop.train.make_predictions(args=args)
    os.remove('pred.npz')
    output_csv = {'title': 'Predicted Values', 'file': 'predicted_values.csv'}
    return {'download':[output_csv]}

def main():
    config={"input_files": "input_smiles.csv"}
    results = app_exe(config)

if __name__ ==  '__main__':
    #set correct working directory
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    main()
