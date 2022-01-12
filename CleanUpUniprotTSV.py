#%%
import pandas as pd
import pickle as pl
#%%

dataUniprot = pd.read_csv("Annotations\Cmerolae_Organism_UniProt_Pull_2021-08-02_master.txt", delimiter="\t")
# dataUniprot = dataUniprot[dataUniprot["Gene names  (ORF )"] != None]
# %%
dataUniprot.shape
#%%
# dataUniprot = dataUniprot[dataUniprot["Gene names  (ORF )"].notnull()]
# dataUniprot.shape
# %%
dataUniprot.columns
# %%
columnsOfInterest = ["Unnamed: 0", "Gene names  (ORF )", 'Protein names', 'Function [CC]']
#%%
dataUniprot = dataUniprot[columnsOfInterest]
# %%

def stringAfterFind(s, q):
    return s[s.find(q)+len(q):]

genesSeen = set()
newData = []
for index, row in dataUniprot.iterrows():
    genes = set()
    
    genes.add(row["Unnamed: 0"])

    ORFcolumn = str(row["Gene names  (ORF )"])

    if "CYME_" in ORFcolumn:
        geneList = ORFcolumn.split(" ")
        for gene in geneList:
            genes.add(stringAfterFind(gene, "CYME_"))

    for gene in genes:
        if not gene in genesSeen:

            if not (pd.isna(row["Function [CC]"])):
                function = stringAfterFind(row["Function [CC]"], "FUNCTION: ")
            else:
                function = row["Function [CC]"]
            newData.append([gene, row["Protein names"], function, genes - {gene}])
        
        genesSeen.add(gene)


    
            # print(stringAfterFind(gene, "CYME_"))
        # print("ORF")
        # print(row["Unnamed: 0"])
        # print(row["Gene names  (ORF )"].split(" "))


    # if len(row["Unnamed: 0"].split(" ")) > 1:
    #     print("Unamed: 0")
    #     print(row["Unnamed: 0"])
    
# %%
df = pd.DataFrame(newData, columns=["GeneID", "Protein", "Function", "Eqivalents"])
# %%
df.to_csv("uniprot_allgenes_2022-01-12.csv", index=False)
# %%
