# %%
import os
import pandas as pd
import re

experimentInfoDir = r"C:\Users\dfoss\Documents\Projects\RaderLab\RaderLabCode\ExperimentInfo"

all_experiments = {}
all_variables = {}

#_Replicate is necessary for the numbering replicates other then that the column data can be what ever you want.
data = pd.read_csv(r"C:\Users\dfoss\Documents\Projects\RaderLab\RNASeqData\AllCombinedRNASeqData.csv", index_col=0)

# %%

def py_initializeExpInfo():
    for root, dirs, files in os.walk(experimentInfoDir):
        for file in files:
            all_experiments[file[:-4]] = {}

            with open(os.path.join(root,file), "r") as csv_file:
                for line in csv_file:
                    key, val = line.split(",")
                    all_experiments[file[:-4]][key] = val.rstrip("\n")
                    if key in all_variables:
                        all_variables[key].add(val.rstrip("\n"))
                    else:
                        all_variables[key] = set([val.rstrip("\n")])


# %%
def py_getExp(filter=None):

    if not filter:
        # global all_experiments
        return sorted(list(all_experiments.keys()))
    
    #Right now filter is just OR of all given elements
    keep = set()
    for val in filter:
        for key in all_experiments.keys():
            if val in key:
                keep.add(val)
    
    return sorted(list(keep))



# %%
def py_getExpInfo(key):
    return(all_experiments[key])




# %%
def py_getExpVariables():
    return sorted(list(all_variables))
# %%
def py_getExpVariableOptions(variable):
    return sorted(list(all_variables[variable]))

# %%

def py_getColumnsOfExp(experiments):
    re_exp = "|".join(experiments)
    col_keep = []
    for col in data.columns:
        if re.match(re_exp, col):
            col_keep.append(col)
    
    col_keep.sort()
    return col_keep
# %%
def py_getExperimentData(experiments):

    col_keep = py_getColumnsOfExp(experiments)

    return data[col_keep]

# %%
def py_getAnnotation(experiments):
    
    col_keep = py_getColumnsOfExp(experiments)

    annontation = pd.DataFrame()
    for exp in experiments:
        for col in col_keep:
            if exp in col:
                hold_s = pd.Series(all_experiments[exp])
                hold_s.name = col
                annontation = annontation.append(hold_s)
    
    annontation = annontation.dropna(axis=1)
    return annontation


# %%
def py_combineVariables(annotation, variables):
    variables = sorted(variables)
    annotation["_".join(annotation)] = annotation[variables].agg("_".join, axis=1) 
    return annotation
# %%
py_initializeExpInfo()
# %%
print(py_getExp())
# %%
print(py_getExpInfo('LP_49h'))

#%%
print(py_getExpVariables())

#%%

print(py_getExpVariableOptions("time"))
# %%
py_getExperimentData(["LP_49h", "RM_49h"])
# %%
df = py_getAnnotation(["LP_49h", "RM_49h"])
df
# %%
py_combineVariables(df, ["time","condition"])
# %%
