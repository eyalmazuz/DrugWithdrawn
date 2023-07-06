# DrugWithdrawn
This repository contains the code for Predicting the withdrawl of drugs from the market as descrived in the paper "Pretrained Transformer Models for Predicting the
Withdrawal of Drugs from the Market".

## Table of Contents

- [Requirements](#requirements)
- [Data](#data)
- [Training](#training)

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
