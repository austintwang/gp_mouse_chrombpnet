import pandas as pd

def parse_scores(scores_path):
    scores_data = pd.read_csv(scores_path, sep='\t', header=0)
    nlq = scores_data["-log10q"]
    f10 = nlq >= 1
    f1 = nlq >= 2
    n = len(nlq)
    
    return scores_data, f10, f1, n