#!/usr/bin/env python

import argparse, time

from rdkit import Chem

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from matplotlib import pyplot


colors = ['red', 'blue', 'green', 'orange', 'cyan', 'purple', 'yellow', 'olive']

def is_active(mol, field_name, field_value, true_value, false_value):
    if mol.HasProp(field_name):
        if field_value:
            val = mol.GetProp(field_name)
            if val == field_value:
                return true_value
            else:
                return false_value
        else:
            return true_value
    else:
        return false_value

def parseValuesTxt(actives, file_name,
                    separator=None, name_column_idx=0, score_column_idx=1, no_header=False, descending=False):

    y = []
    scores = []
    errors = 0
    count = 0
    
    if name_column_idx == None:
        name_column_idx = 0
    if score_column_idx == None:
        score_column_idx = 1
    
    with open(file_name, 'r') as data:
        if not no_header:
            data.readline()
        while True:
            line = data.readline()
            if not line:
                break
                
            tokens = line.split(separator)
            name = tokens[int(name_column_idx)]
            score = tokens[int(score_column_idx)]
            if name in actives:
                y.append(1)
            else:
                y.append(0)
            
            if descending:    
                scores.append(float(score))
            else:
                scores.append(float(score) * -1.0)
            count += 1
        
    return y, scores
         

def parseValuesSdf(actives, poses_file_name,
                    active_field_name=None, active_field_value=None,
                    inactive_field_name=None, inactive_field_value=None, 
                    name_field_name=None, score_field_name=None,
                    descending=False):

    if not actives and not active_field_name and not inactive_field_name:
        raise ValueError('Must specify one of actives-file-name, active-field-name or inactive-field-name')

    if active_field_name and inactive_field_name:
        raise ValueError('Must specify one of active_field_name or inactive_field_name, not both')
    
    if actives and not name_field_name:
        raise ValueError('When using an actives file and data is a SD file you must specify the name-field-name argument')

    supplr = Chem.SDMolSupplier(poses_file_name)
    y = []
    scores = []
    errors = 0
    count = 0
    for mol in supplr:
        if mol.HasProp(score_field_name):
            score = mol.GetDoubleProp(score_field_name)
            if descending:    
                scores.append(score)
            else:
                scores.append(score * -1.0)
        else:
            print('No score field for record', count)
            continue;
            
        if actives:
            if mol.HasProp(name_field_name):
                name = mol.GetProp(name_field_name)
                if name in actives:
                    y.append(1)
                else:
                    y.append(0)
            else:
                print('No name field for record', count)
        elif active_field_name:
            is_active_val = is_active(mol, active_field_name, active_field_value, 1, 0)
            y.append(is_active_val)
        elif inactive_field_name:
            is_active_val = is_active(mol, inactive_field_name, inactive_field_value, 0, 1)
            y.append(is_active_val)
             
        count += 1
        
    return y, scores


def add_curve(index, actives, data_file_name, label, 
                active_field_name=None, active_field_value=None,
                inactive_field_name=None, inactive_field_value=None,
                name_field=None, score_field=None, descending=False,
                color=None, no_header=False, separator=None):
    
    if color:
        col = color
    else:
        col = colors[index]
    
    if data_file_name.endswith('.sdf') or data_file_name.endswith('.sd'):
        # SDF
        y, scores = parseValuesSdf(actives, data_file_name, 
            active_field_name=active_field_name, active_field_value=active_field_value, 
            inactive_field_name=inactive_field_name, inactive_field_value=inactive_field_value,
            name_field_name=name_field, score_field_name=score_field, descending=descending)
    else:
        # CSV etc.
        y, scores = parseValuesTxt(actives, data_file_name,
            name_column_idx=name_field, score_column_idx=score_field,
            no_header=no_header, separator=separator, descending=descending)
        
    fpr, tpr, thresholds = roc_curve(y, scores, pos_label=1)
    auc = roc_auc_score(y, scores)
    
    pyplot.plot(fpr, tpr, col, label=label + ' (AUC=' + str(round(auc, 2)) + ')')
    
    return y, scores
    

def read_actives(filename):
    with open(filename, 'r') as f:
        actives = f.read().splitlines()
    return actives
        

def main():

    # Example:
    #   python3 calculate-roc-curves.py --actives-file-name actives.txt\
    #     -p1 results_rdock/results_1poseperlig.sdf  --name-field-name1 _Name --score-field-name1 SCORE.norm -l1 rDock

    ### command line args definitions #########################################

    parser = argparse.ArgumentParser(description='Prepare rDock docking')
    parser.add_argument('--no-diagonal', action='store_true', help="Don't show the diagonal line showing random performance")
    parser.add_argument('-a', '--actives-file-name', help='File with actives')
    parser.add_argument('-o', '--output-file-name', help='File name for output')
    parser.add_argument('-s', '--figure-size', type=float, default=5, help='Figure size in inches')
    
    max_curves = 9
    for n in range(1, max_curves):
        s = str(n)
        parser.add_argument('-p'+s, '--poses-file-name'+s, help='Poses SDF')
        parser.add_argument('--active-field-name'+s, help='Field name used to determine if record is an active')
        parser.add_argument('--active-field-value'+s, help='Optional field value used to determine if record is an active')
        parser.add_argument('--inactive-field-name'+s, help='Field name used to determine if record is an inactive')
        parser.add_argument('--inactive-field-value'+s, help='Optional field value used to determine if record is an inactive')
        parser.add_argument('--name-field'+s, 
            help='Optional field that contains the compound name (SDF) or index (TXT) that is present when using the --actives-file-name argument')
        parser.add_argument('--score-field'+s, help='Field name (SDF) or index (TXT) used for the score')
        parser.add_argument('--separator'+s, help='Separator for TXT files. Default is whitespace')
        parser.add_argument('--no-header'+s, action='store_true', help='No header line is present for TXT files')
        parser.add_argument('--descending'+s, action='store_true', help='Rank the scores in descending order')
        parser.add_argument('-l'+s, '--label'+s, help='Label for curve')
        parser.add_argument('-c'+s, '--color'+s, help='Color for the cuve')
   
    args = parser.parse_args()
    #print("calculate-roc-curves: ", args)
       
    if args.actives_file_name:
        actives = read_actives(args.actives_file_name)
    else:
        actives = None

    t0 = time.time()
    for n in range(1, max_curves):
        s = str(n)
        if getattr(args, 'poses_file_name' + s, None):
            print('Processing', s, getattr(args, 'poses_file_name' +s))
            y, scores = add_curve(
                n - 1, actives,
                getattr(args, 'poses_file_name' +s), 
                active_field_name=getattr(args, 'active_field_name' +s, None),
                active_field_value=getattr(args, 'active_field_value' +s, None),
                inactive_field_name=getattr(args, 'inactive_field_name' +s, None),
                inactive_field_value=getattr(args, 'inactive_field_value' +s, None),
                name_field=getattr(args, 'name_field' +s, None),
                score_field=getattr(args, 'score_field' +s, None),
                descending=getattr(args, 'descending' +s),
                label=getattr(args, 'label' +s, None),
                color=getattr(args, 'color' +s, None),
                no_header=getattr(args, 'no_header' +s),
                separator=getattr(args, 'separator' +s, None)
                )
    
    if not args.no_diagonal:
        scores0 = [0 for _ in range(len(y))]
        fpr0, tpr0, thresholds = roc_curve(y, scores0, pos_label=1)
        pyplot.plot(fpr0, tpr0, linestyle='--', color='grey', label='Random')
    
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    # show the legend
    pyplot.legend(loc='lower right')
    # show or save the plot
    if args.output_file_name:
        fig = pyplot.gcf()
        fig.set_size_inches(args.figure_size, args.figure_size)
        fig.savefig(args.output_file_name)
    else:
        pyplot.show()
    
    
    t1 = time.time()
    
    
    
    
if __name__ == "__main__":
    main()
    
