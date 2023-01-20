import numpy as np
from sklearn.metrics import roc_auc_score, f1_score, classification_report, confusion_matrix, average_precision_score, precision_recall_curve, auc, cohen_kappa_score, roc_curve


def senstivity(y_test, test_preds_class):
    tn, fp, fn, tp = confusion_matrix(y_test, test_preds_class).ravel()
    sensitivity = tp / (tp + fn)
    return sensitivity
def specficity(y_test, test_preds_class):
    tn, fp, fn, tp = confusion_matrix(y_test, test_preds_class).ravel()
    specificity = tn / (tn + fp)
    return specificity

def find_best_threshold(y_true, y_prob):
    # calculate roc curves
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    # get the best threshold
    J = tpr - fpr
    ix = np.argmax(J)
    best_thresh = thresholds[ix]
    return best_thresh

class Scorer():
    def __init__(self):
   
        self.results = {'Name':[],'Model':[],'AUPR':[], 'AUC':[], 'Kappa':[],'Threshold':[],'Modalities':[],'CV':[],'Precision':[], 'Recall':[],'F1-score':[],'Precision-macro':[], 'Recall-macro':[],'F1-score-macro':[], 'Precision-weighted':[], 'Recall-weighted':[],'F1-score-weighted':[], 'Specificity':[],'Sensitivity':[] }


    def add_score(self,modalities,clfName,y_test,test_preds,cv_i,best_threshold):
        auc_score = roc_auc_score(y_test, test_preds)

        self.results['Modalities'].append(modalities)
        self.results['Name'].append(clfName)
        self.results['Model'].append(clfName)
        self.results['CV'].append(cv_i)
        self.results['Threshold'].append(best_threshold)
        self.results['AUC'].append(auc_score)
        test_preds_class = [0 if x < best_threshold else 1 for x in test_preds]
        
        precision, recall, _ = precision_recall_curve(y_test, test_preds)
        aupr = auc(recall, precision)
        self.results['AUPR'].append(aupr)

        report = classification_report(y_test, test_preds_class, output_dict=True)

        self.results['Precision'].append(report['1']['precision'])
        self.results['Recall'].append(report['1']['recall'])  # ==sensitivity
        self.results['F1-score'].append(report['1']['f1-score'])
        
        self.results['Precision-macro'].append(report['macro avg']['precision'])
        self.results['Recall-macro'].append(report['macro avg']['recall'])  # ==sensitivity
        self.results['F1-score-macro'].append(report['macro avg']['f1-score'])
        
        self.results['Precision-weighted'].append(report['weighted avg']['precision'])
        self.results['Recall-weighted'].append(report['weighted avg']['recall'])  # ==sensitivity
        self.results['F1-score-weighted'].append(report['weighted avg']['f1-score'])

        specificity = specficity(y_test, test_preds_class)
        sensitivity = senstivity(y_test, test_preds_class)

        kappa = cohen_kappa_score(y_test, test_preds_class)
        self.results['Kappa'].append(kappa)
        self.results['Specificity'].append(specificity)
        self.results['Sensitivity'].append(sensitivity)