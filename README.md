# DrugWithdrawn
This repository contains the code for Predicting the withdrawl of drugs from the market as descrived in the paper "Pretrained Transformer Models for Predicting the
Withdrawal of Drugs from the Market".

## Table of Contents

- [Requirements](#requirements)
- [Data](#data)
- [Training](#training)
- [Feature Extraction](#extracting-features-from-pre-trained-models)

## Requirements
* [Huggingface Transformers](https://huggingface.co/docs/transformers/index)- used for ChemBERTA v1 and v2
* [T5Chem](https://github.com/HelloJocelynLu/t5chem)
* [Chemprop](https://github.com/chemprop/chemprop)
* [XGBoost](https://xgboost.readthedocs.io/en/stable/)
* [Scikit-learn](scikit-learn.org/)
* [Mordred](https://github.com/mordred-descriptor/mordred)- used for XGBoost and SVM models
* [DeepChem](https://deepchem.io/)- used for XGBoost and SVM models
* [Mol2Vec](https://github.com/samoturk/mol2vec)- used for XGBoost and SVM models

## Data
The data is provided under the ``split`` folder. There 4 for evaluation schemes in general and 2 are used in the paper
1. Agree No Dups - All data sources agree on the labels, a molecule must appear either in the train or test set
2. Agree Dups - All data sources agree on the labels, a molecule can appear in both the train and test set
3. No Agree No Dups - Data sources might have different labels for the same molecule, a molecule must appear either in the train or test set
4. No Agree Dups - Data sources might have different labels for the same molecule, a molecule can appear in both the train and test set

In each evaluation scheme we employ a leave one out validation, We remove entirely one of the data sources and use it as our test set,
other data sources are used for training.

Each folder inside the evaluation indicates which data source is used for testing. i.e. if the folder named DrugBank, then DrugBank is used as the test set and ChEMBL, NCATS and WITHDRAWN are used as training set.

## Training

each notebook is used to train a different model from the paper.
* ChemBERTa- is used to train and evaluate [ChemBERTa v1](https://huggingface.co/seyonec/PubChem10M_SMILES_BPE_450k)
* ChemBERTa-2- is used to train and evaluate [ChemBERTa v2 MTR](https://huggingface.co/DeepChem/ChemBERTa-77M-MTR) [ChemBERTa v2 MLM](https://huggingface.co/DeepChem/ChemBERTa-77M-MLM)
* PaperModels- is used to train and evaluate XGBoost and SVM models using different handcrafted features
* T5chem_model- is used to train and evaluate [T5chem](https://github.com/HelloJocelynLu/t5chem) (follow the instruction in the repo to how to download the pre-trained model that is used for fine-tune)

To run training for each of notebooks and follow the code written inside

## Extracting features from pre-trained models

If one desire to use pre-trained models such as ChemBERTa and Chemprop to extract fingerprint vectors for each SMILES and then use a classical machine learning algorithm for prediciton (such in paper models mixed features section or SMILES only section).

all the pre-trained features already exists in their respective data directories.

But if one wish to extract their own features it can be done as follows:

1. Extract ChemBERTa features.ipynb notebook will extract ChemBERTa features, the get_embeddings method receives a DataFrame which needs to contain a column called ``smiles``. then we load a pre-trained ChemBERTa model and iterate over the dataset and extract the hidden representation of the model for each smile then save it back to the same location as the original data.

2. To extract Chemprop features, please follow their documentation on how to train a model here [training](https://github.com/chemprop/chemprop/#training) and how to extract Chemprop fingerprints here [fingerprints](https://github.com/chemprop/chemprop/#encode-fingerprint-latent-representation).
A pre-trained Chemprop model is already provised in the Chemprop directory in additional to features extract from a pre-trained models on Tox21 dataset, this model was trained on Tox21 dataset which tries to predict toxological properties of molecules (which could indicate a reason for withdrawing a drug from the market, which is a very close task). and the Tox21 dataset is also provided in the same directory.
If one wish to use a different dataset, they can download it from the moleculeNet repository here: [datasets](https://moleculenet.org/datasets)

3. Mol2Vec features can be extracted through the PaperModels.ipynb notebook and then manually save them to the disk, additionaly same as previous extract features, they are already provided in the data directories for use and no need to extract them again.
