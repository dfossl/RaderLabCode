import csv
import pickle
import pandas as pd
import numpy as np
import random
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import matplotlib
from venn import venn
import os
import global_python_meta
random.seed(123)
np.random.seed(123)



#
# asdf = [1,2,3]
#
# for i in asdf:
#     print(i < 2)

# date_rng = pandas.date_range('2015-01-01', periods=200, freq='D')
#
# df1 = pandas.DataFrame(numpy.random.randn(100, 3), columns='A B C'.split(), index=date_rng[:100])
#
# print(df1.index[0])
#
# for index, row in df1.iterrows():
#     print(index)

annotationDir = global_python_meta.annotationDir

# links each go term to how many genes mapp to that pathway
pickle_in = open(os.path.join(annotationDir, "goTerm2NumberofGenes_Dict.pickle"), "rb")
dict_GOterm2NumberGenes = pickle.load(pickle_in)

# Links each gene to a list of pathway IDs
pickle_in = open(os.path.join(annotationDir, "gene2GOTermList_Dict.pickle"), "rb")
dict_GeneToListPathway = pickle.load(pickle_in)

# Gets description of any goterm
pickle_in = open(os.path.join(annotationDir, "Term2Description_UNIVERSAL.pickle"), "rb")
dict_GoTerm2description = pickle.load(pickle_in)

pickle_in = open(os.path.join(annotationDir, "GeneID_2_Kegg.pl"), "rb")
dict_GeneID_2_Kegg = pickle.load(pickle_in)

pickle_in = open(os.path.join(annotationDir, "dict_KEGG_pathway_2_genelist.pl"), "rb")
dict_KEGG_pathway_2_genelist = pickle.load(pickle_in)

pickle_in = open(os.path.join(annotationDir, "list_allKeggGenesInPathways.pl"), "rb")
list_allKeggGenesInPathways = pickle.load(pickle_in)

genesToProteindf = pd.read_csv(os.path.join(annotationDir, "uniprot_allgenes_2022-01-12.csv"), index_col=0)


def sortCombinedFileUpandDowngenesToList(df_LFC, score_thres, inputDir):
    groups = set()
    for col in df_LFC.columns:
        holderList = col.rsplit("_")
        # print(holderList)
        holderList.pop(0)
        group = "_".join(holderList)
        groups.add(group)
    # return groups
    dictGroupToUpregGeneList = {}
    dictGroupToDownregGeneList = {}
    for group in groups:
        dictGroupToUpregGeneList[group] = []
        dictGroupToDownregGeneList[group] = []
    # print(dictGroupToUpregGeneList)
    # print(dictGroupToDownregGeneList)

    for index, row in df_LFC.iterrows():
        for group in groups:
            if row["LFC_" + group] >= 0 and row["Score_" + group] >= score_thres:
                dictGroupToUpregGeneList[group].append(index)
            elif row["LFC_" + group] < 0 and row["Score_" + group] >= score_thres:
                dictGroupToDownregGeneList[group].append(index)

    fileNames = []
    for group in groups:
        filename = os.path.join(inputDir, group + "_" + str(score_thres) + "_listUpreg.pickle")
        fileNames.append(filename)
        pickle_out = open(filename, "wb")
        pickle.dump(dictGroupToUpregGeneList[group], pickle_out)
        pickle_out.close()

    for group in groups:
        filename = os.path.join(inputDir, group + "_" + str(score_thres) + "_listDownreg.pickle")
        fileNames.append(filename)
        pickle_out = open(filename, "wb")
        pickle.dump(dictGroupToDownregGeneList[group], pickle_out)
        pickle_out.close()

    return fileNames


def sortDeseqResultUpandDowngenesToCSV(df_LFC, thres, inputDir):
    df_Upreg = pd.DataFrame(columns = df_LFC.columns)
    df_Downreg = pd.DataFrame(columns = df_LFC.columns)

    for index, row in df_LFC.iterrows():
        if row["log2FoldChange"] >= thres:
            df_Upreg.append(other=row)
        elif row["log2FoldChange"] <= -1*thres:
            df_Downreg.append(other=row)

    df_Upreg.to_csv(os.path.join(inputDir, "OnlyUpreg.csv"))
    df_Downreg.to_csv(os.path.join(inputDir, "OnlyDownreg.csv"))


def sortDesqResultsUpandDownToList(df_LFC, thres, filename, inputDir):
    UpregList = []
    DownregList = []

    for index, row in df_LFC.iterrows():
        if row["log2FoldChange"] >= thres:
            UpregList.append(index)
        elif row["log2FoldChange"] <= -1*thres:
            DownregList.append(index)

    upregFileName = os.path.join(inputDir, filename +"_listUpreg.pickle")
    pickle_out = open(upregFileName, "wb")
    pickle.dump(UpregList, pickle_out)
    pickle_out.close()

    downregFileName = os.path.join(inputDir, filename + "_listDownreg.pickle")
    pickle_out = open(downregFileName, "wb")
    pickle.dump(DownregList, pickle_out)
    pickle_out.close()

    return [upregFileName, downregFileName]

def pythonEnrichmentOnPickleGeneList(filename):
    pickle_in = open(filename, "rb")
    geneList = pickle.load(pickle_in)
    holderList = filename.rsplit(".")
    holderList.pop(-1)
    enrichFilename = ".".join(holderList)
    outfilename = enrichFilename + '_EnrichmentResults.csv'
    with open(outfilename, 'w', newline="") as f:
        w = csv.DictWriter(f, ["GoTerm", "Decription", "Pvalue", "Genes", "k,T,m,q"])
        w.writeheader()

        T = 0
        m = 0
        q = 0
        k = 0
        for GOterm in dict_GOterm2NumberGenes:
            # print(GOterm)
            if GOterm in dict_GoTerm2description:
                pathwayDescription = dict_GoTerm2description[GOterm]
            else:
                pathwayDescription = "fix"
            # print(pathwayDescription)
            DEGgenesInPathway = 0
            listOFGenesInGOterm = []
            for gene in geneList:
                if gene in dict_GeneToListPathway:
                    if GOterm in dict_GeneToListPathway[gene]:
                        DEGgenesInPathway = DEGgenesInPathway + 1
                        listOFGenesInGOterm.append(gene)

            T = 5009
            m = len(geneList)
            q = int(dict_GOterm2NumberGenes[GOterm])
            k = DEGgenesInPathway - 1
            result = hypergeom.sf(k, T, m, q)
            values = str(k) + " " + str(T) + " " + str(m) + " " + str(q)
            w.writerow(
                {'GoTerm': GOterm, "Decription": pathwayDescription, 'Pvalue': result, 'Genes': listOFGenesInGOterm,
                 "k,T,m,q": values})
    return outfilename

def pythonKeggEnrichmentOnRVector(geneList, filename):
    """
    Should work
    """
    outfilename = filename + '_KEGGEnrichmentResults.csv'
    with open(outfilename, 'w', newline="") as f:
        w = csv.DictWriter(f, ["Pathway", "Decription", "Pvalue", "Genes", "k,T,m,q"])
        w.writeheader()
        T = len(list_allKeggGenesInPathways)
        m = len(geneList)        
        q = 0
        k = 0

        for kegg_path in dict_KEGG_pathway_2_genelist:
            pathwayDescription = kegg_path.split(' ',1)[1]


            # DEGgenesInPathway = 0
            listOFGenesInPath = []

            for gene in geneList:
                if gene in dict_KEGG_pathway_2_genelist[kegg_path]:
                    # DEGgenesInPathway = DEGgenesInPathway + 1
                    listOFGenesInPath.append(gene)



            q = len(dict_KEGG_pathway_2_genelist[kegg_path])
            k = len(listOFGenesInPath) - 1
            result = hypergeom.sf(k, T, m, q)
            values = str(k) + " " + str(T) + " " + str(m) + " " + str(q)
            w.writerow(
                {'Pathway': kegg_path.split(' ',1)[0], "Decription": pathwayDescription, 'Pvalue': result, 'Genes': listOFGenesInPath,
                 "k,T,m,q": values})
    return outfilename

def py_getKEGGDataframe(df):
    """
    Change index to KEGG?

    Args:
        df ([type]): [description]

    Returns:
        [type]: [description]
    """
    # return 1
    new_indeices = set(list_allKeggGenesInPathways).intersection(set(df.index))
    new_df = df.loc[new_indeices,:]

    old_indx = list(new_df.index)
    for i in range(len(old_indx)):
        old_indx[i] = dict_GeneID_2_Kegg[old_indx[i]]
    
    new_df.index = old_indx

    return new_df

def py_GetAllPathwaysOfSig(filenames, p_cutoff=0.05):
    
    pathways = set()
    for name in filenames:
        df = pd.read_csv(name)

        for index, row in df.iterrows():
            if row["Pvalue"] < p_cutoff:
                pathways.add(row["Pathway"])
    
    return list(pathways)


def pythonEnrichmentOnRVector(geneList, filename):

    outfilename = filename + '_EnrichmentResults.csv'
    with open(outfilename, 'w', newline="") as f:
        w = csv.DictWriter(f, ["GoTerm", "Decription", "Pvalue", "Genes", "k,T,m,q"])
        w.writeheader()
        T = 0
        m = 0
        q = 0
        k = 0
        for GOterm in dict_GOterm2NumberGenes:
            # print(GOterm)
            if GOterm in dict_GoTerm2description:
                pathwayDescription = dict_GoTerm2description[GOterm]
            else:
                pathwayDescription = "fix"
            # print(pathwayDescription)
            DEGgenesInPathway = 0
            listOFGenesInGOterm = []
            for gene in geneList:
                if gene in dict_GeneToListPathway:
                    if GOterm in dict_GeneToListPathway[gene]:
                        DEGgenesInPathway = DEGgenesInPathway + 1
                        listOFGenesInGOterm.append(gene)

            T = 5009
            m = len(geneList)
            q = int(dict_GOterm2NumberGenes[GOterm])
            k = DEGgenesInPathway - 1
            result = hypergeom.sf(k, T, m, q)
            values = str(k) + " " + str(T) + " " + str(m) + " " + str(q)
            w.writerow(
                {'GoTerm': GOterm, "Decription": pathwayDescription, 'Pvalue': result, 'Genes': listOFGenesInGOterm,
                 "k,T,m,q": values})
    return outfilename


"""
df is a dataframe that is returned from DESeq2::results

variable is the statistic you want to color the dots on

comparison, a string for the final file name

threshold the value for the vaiable you want to use as the threshold for coloring.

if a percentile is provided the the threshold is created from that percentile
"""
def py_generateMAPlot(df, outputfiledir, variable, comparison, threshold=None, percentile=None):

    percentile = percentile*100

    df = df.dropna()
    n = len(df.index)
    values = list(df[variable].values)
    # values.sort()

    if (threshold is None) and (percentile is None):
        return("Error")
    elif threshold:
        pass
    elif percentile:
        threshold = np.nanpercentile(values, percentile)

    colors = []

    for index, row in df.iterrows():
        value = row[variable]

        if value > threshold:
            LFC = row["log2FoldChange"]
            if LFC > 0:
                colors.append("UP")
            else:
                colors.append("DOWN")
        else:
            colors.append("Remaining")

    df["color"] = colors
    x = (df.loc[:,"baseMean"].values)
    x = np.log(x)
    y = (df.loc[:,"log2FoldChange"].values)

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(figsize=(5*3, 5))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    # the scatter plot:
    targets = ["Remaining", "UP", "DOWN"]
    colors = ["grey", "r", "skyblue"]
    for target, color in zip(targets, colors):
        indicesToKeep = df['color'] == target
        ax_scatter.scatter(np.log(df.loc[indicesToKeep, 'baseMean'])
                   , df.loc[indicesToKeep, 'log2FoldChange']
                   , c=color
                   , s=8)
        # print(t)
    children = ax_scatter.get_children()
    ax_scatter.legend(handles=[children[1], children[2]], labels=["Up","Down"])
    # ax_scatter.legend(["UP", "DOWN"])


    # now determine nice limits by hand:
    binwidth = 0.25
    # print(np.abs([x, y]))
    # print(np.abs([x, y]).max())
    x_lim = np.ceil(np.abs(x).max() / binwidth) * binwidth
    # x_lim = 20000

    y_lim = np.ceil(np.abs(y).max() / binwidth) * binwidth

    ax_scatter.set_xlim((-0, x_lim))
    ax_scatter.set_ylim((-y_lim, y_lim))

    ax_scatter.set_xlabel("Natural log Counts")
    ax_scatter.set_ylabel("LFC")



    x_bins = np.arange(-0, x_lim + binwidth, binwidth)
    # print(x_bins)
    y_bins = np.arange(-y_lim, y_lim + binwidth, binwidth)

    ax_histx.hist(x, bins=x_bins, histtype='stepfilled', fc="grey", ec="grey")
    ax_histy.hist(y, bins=y_bins, histtype='stepfilled', orientation='horizontal', fc="grey", ec="grey")

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())

    plt.savefig(os.path.join(outputfiledir, comparison+"MA_plot.png"), dpi=300)
    plt.close('all')



"""
df is datafraem from DESeq2::results
"""
def py_pvalueVSbasemean(df, outputfiledir, comparison, lfc_transform = "max"):

    df = df.dropna()
    # print(df)
    holder = -np.log(df["padj"])

    largest = holder[holder != np.Inf].max()
    # holder = holder.replace(np.inf, largest)
    holder = np.where(holder==np.inf, largest, holder)

    lfc_values = df.loc[:,"log2FoldChange"].abs().values

    if lfc_transform == "max":
        percentile_99 = (np.nanpercentile(lfc_values, 99))
        label = f"Absolute LFC* ({np.sum(lfc_values > percentile_99)})"

        df.loc[df['log2FoldChange'].abs() > percentile_99, "log2FoldChange"] = percentile_99
    elif lfc_transform == "arctanh":
        df['log2FoldChange'] = np.arctanh(df['log2FoldChange'])
        label = "Absolute arctanh(LFC)"
    elif lfc_transform == "log":
        df['log2FoldChange'] = np.log(df['log2FoldChange'].abs())
        label = "log(Absolute LFC)"

    # print(df)

    fig, ax = plt.subplots()
    width = 12
    fig.set_size_inches(width,.33*width)
    fig.set_tight_layout(True)
    sctr = ax.scatter(np.log(df["baseMean"]), holder, c=abs(df["log2FoldChange"]), cmap='Reds')
    ax.set_xlabel("log(baseMean)")
    ax.set_ylabel("-log(padj)")
    plt.colorbar(sctr, ax=ax, label=label)
    # plt.show()
    plt.savefig(os.path.join(outputfiledir, comparison + "_LFC_padj_basemean.png"),dpi=250)
    plt.close(fig)


"""
makes histrogram for a given variable of the given dataframe.

Subset values based on percentile.
"""
def py_makeHistPlot(df, outputfiledir, variable, percentile, comparison):
    percentile = percentile*100
    holder = (df[variable])

    threshold = np.nanpercentile(holder, percentile)

    holder = holder[holder > threshold]

    largest = holder[holder != np.Inf].max()
    holder = holder.replace(np.inf, largest)

    sns.distplot((holder), hist=True, kde=False,
             color = 'darkblue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})

    plt.title(comparison + " " + str(percentile) + " percentile")

    plt.savefig(os.path.join(outputfiledir, "histogram_"+comparison+"_"+str(percentile)+"_percentile_"+variable+".png"))
    plt.close("all")


def py_makeHistPlotFromList(df, outputfiledir, variable, index , label):

    holder = df.loc[index, variable]

    largest = holder[holder != np.Inf].max()
    holder = holder.replace(np.inf, largest)

    sns.distplot((holder), hist=True,
             kde=False,
             color = 'darkblue',
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})

    plt.title(label+" "+variable)

    plt.savefig(os.path.join(outputfiledir, "histogram_"+label+"_"+"_percentile_"+variable+".png"))
    plt.close("all")



"""
Takes the combined results dataframe from all comparisons and makes a density plot for given columns.

columns is a list of columns to include in the density Plot
"""
def py_densityplotCompareColumns(df, columns, percentile, comparison, xlabel, outputfiledir, legend=False, expLabel=""):

    percentile = percentile*100

    return_data = {}
    i = 0
    #This is needed since single vectors from R get passed as strings
    #Since we know that no column name should be less then 2 we can check
    #if this is the case we know its a single column string.
    #not best way to handle but works for now
    if len(columns[0]) < 2:
        column = columns
	# return(percentile)
        holder = (df[column])

        # return holder

        threshold = np.nanpercentile(holder, percentile)

        # return threshold

        holder = holder[holder > threshold]

        largest = holder[holder != np.Inf].max()
        holder = holder.replace(np.inf, largest)

        ax =sns.distplot(holder, hist = False, kde = True,
                     kde_kws = {'linewidth': 3},
                     label = column)


        return_data[column+"_x"] = ax.get_lines()[i].get_data()[0]
        return_data[column+"_y"] = ax.get_lines()[i].get_data()[1]
    else:
        for column in columns:

            # return(percentile)
            holder = (df[column])
    
            # return holder
    
            threshold = np.nanpercentile(holder, percentile)
    
            # return threshold
    
            holder = holder[holder > threshold]
     
            largest = holder[holder != np.Inf].max()
            holder = holder.replace(np.inf, largest)
    
            ax =sns.distplot(holder, hist = False, kde = True,
                         kde_kws = {'linewidth': 3},
                         label = column)


            return_data[column+"_x"] = ax.get_lines()[i].get_data()[0]
            return_data[column+"_y"] = ax.get_lines()[i].get_data()[1]
            
            i += 1
        
    # Plot formatting
    if legend:
        
        ax.legend(prop={'size': 14}, title = 'Group')

        # Shrink current axis by 20%
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        leg_plt = plt.legend(prop={'size': 16}, title = 'Group')
        plt.legend().set_visible(False)
    
    plt.title('Density Plot ' + comparison + " " + str(percentile) + " percentile")
    plt.xlabel(xlabel)
    plt.ylabel('Density')
    plt.tight_layout()
    if expLabel:
        plt.savefig(os.path.join(outputfiledir, expLabel+"_denisty_compare_"+comparison+"_"+str(percentile)+"_percentile_"+xlabel+".png"))
    else:
        plt.savefig(os.path.join(outputfiledir, "denisty_compare_"+comparison+"_"+str(percentile)+"_percentile_"+xlabel+".png"))

    
    plt.close("all")

    return return_data


"""
Gets the indexes from a df that satisfy the precentile threshold for a particular variable.
"""
def py_sortcolumns(df, percentile, variable):

    percentile = percentile*100

    index_keep = set()

    for column in df.columns:
        if variable.lower() in column.lower():
            keep = df[df[column] > np.nanpercentile(df[column], percentile)]
            index_keep = index_keep.union(set(keep.index))

    return(list(index_keep))


"""
returns df with only the genes specified.
"""
def py_subsetDataOnGenes(df, genesToKeep):

    df = df.loc[genesToKeep,:]
    return df


"""
returns a dataframe thats normalized on a per row basis, depending on a refernece column.
"""
def py_normalizeOnVariable(df, variable, reference_column):

    columns = []
    for column in df.columns:
        if variable.lower() in column.lower():
            columns.append(column)

    df = df[columns]

    for index, row in df.iterrows():

        ref_LFC = row[reference_column]
        # df.set_value(index, reference_column, 0)
        df.at[index, reference_column] = 0

        for column in columns:
            if column != reference_column:
                # df.set_value(index, column, (row[column])-(ref_LFC))
                df.at[index, column] = (row[column])-(ref_LFC)

    return df


"""
Genrates a elbow plot for the provided data.
"""
def py_elbowPlot(df, variable, outputdir, label, max_k=20):

    max_k = int(max_k)
    # df = df.loc[genesToKeep,:]

    columns = []
    for column in df.columns:
        if variable.lower() in column.lower():
            columns.append(column)

    df = df[columns]

    # return(list(df.index))

    mms = MinMaxScaler()
    mms.fit(df)
    data_transformed = mms.transform(df)
    #
    Sum_of_squared_distances = []
    K = range(1,max_k)
    for k in K:
        km = KMeans(n_clusters=k, random_state=123)
        km = km.fit(data_transformed)
        Sum_of_squared_distances.append(km.inertia_)
    plt.plot(K, Sum_of_squared_distances, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Sum_of_squared_distances')
    plt.title('Elbow Method For Optimal k')
    plt.savefig(os.path.join(outputdir, variable+"_"+label+"_elbowPlot.png"))
    plt.close("all")


def py_optimalK(data, variable, nrefs=3, maxClusters=15):
    """
    Calculates KMeans optimal K using Gap Statistic from Tibshirani, Walther, Hastie
    Params:
        data: ndarry of shape (n_samples, n_features)
        nrefs: number of sample reference datasets to create
        maxClusters: Maximum number of clusters to test for
    Returns: (gaps, optimalK)
    """

    nrefs=int(nrefs)
    maxClusters=int(maxClusters)

    # data = data.loc[genesToKeep,:]

    columns = []
    for column in data.columns:
        if variable.lower() in column.lower():
            columns.append(column)

    data = data[columns]


    gaps = np.zeros((len(range(1, maxClusters)),))
    resultsdf = pd.DataFrame({'clusterCount': [], 'gap': []})
    for gap_index, k in enumerate(range(1, maxClusters)):

        # Holder for reference dispersion results
        refDisps = np.zeros(nrefs)

        # For n references, generate random sample and perform kmeans getting resulting dispersion of each loop
        for i in range(nrefs):
            # Create new random reference set
            randomReference = np.random.random_sample(size=data.shape)

            # Fit to it
            km = KMeans(k, random_state=123)
            km.fit(randomReference)

            refDisp = km.inertia_
            refDisps[i] = refDisp

        # Fit cluster to original data and create dispersion
        km = KMeans(k, random_state=123)
        km.fit(data)

        origDisp = km.inertia_

        # Calculate gap statistic
        gap = np.log(np.mean(refDisps)) - np.log(origDisp)

        # Assign this loop's gap statistic to gaps
        gaps[gap_index] = gap

        resultsdf = resultsdf.append({'clusterCount': k, 'gap': gap}, ignore_index=True)

    return([resultsdf, gaps.argmax() + 1])  # Plus 1 because index of 0 means 1 cluster is optimal, index 2 = 3 clusters are optimal


"""
generates kmeans on provided data and columns
"""
def py_performKmeans(df, variable, k):

    k = int(k)

    # df = df.loc[genesToKeep,:]

    columns = []
    for column in df.columns:
        if variable.lower() in column.lower():
            columns.append(column)

    df = df[columns]


    x = df.loc[:,:].values

    scaler = StandardScaler()
    scaler.fit(x)
    x = scaler.transform(x)
    # print(x)
    # n_clusters = 17
    model = KMeans(n_clusters=k, random_state=123)
    model.fit(x)
    clust_labels = model.predict(x)

    # print(clust_labels)
    df["cluster"] = pd.Series(clust_labels, index=df.index)

    return(df)
    # df.to_csv(str(n_clusters) + "_Union_top10percent_kmeans.csv")


def py_visualizeKmeans(df, variable, outputdir, label=""):

    n_clusters = df["cluster"].max() + 1

    columns = []
    for column in df.columns:
        if variable.lower() in column.lower():
            columns.append(column)

    for i in range(0,n_clusters):
        df2 = df[df.loc[:,"cluster"] == i]
        df2 = df2[columns]



        data= {"timepoint":[], "LFC":[], "gene":[]}

        # counter = 0
        for index, row in df2.iterrows():
            # counter += 1
            for column in df2.columns:
                data["timepoint"].append(column)
                data["LFC"].append(row[column])
                data["gene"].append(index)
            # if counter > 50:
            #     break
            # print(data)


        ax = sns.pointplot(x="timepoint",y="LFC",hue="gene",data=data)
        ax.get_legend().remove()
        if label != "":
            plt.savefig(os.path.join(outputdir, str(n_clusters) + "_subset_cluster_" + str(i) +"_"+label+"_.png"))
        else:
            plt.savefig(os.path.join(outputdir, str(n_clusters) + "_subset_cluster_" + str(i) + ".png"))

        plt.close("all")


"""
performs enrichment on kmeans clusters
"""
def py_enrichmentOnKmeans(df):
    n_clusters = df["cluster"].max() + 1

    for i in range(0, n_clusters):
        genes = df[df["cluster"] == i].index
        label = str(n_clusters) + "_cluster_" + str(i)
        pythonEnrichmentOnRVector(genes, label)
        # sort.generateFileForKegg(df, i)


def py_countPlotKmeans(df, outputdir, label=""):

    n_clusters = int(df["cluster"].max() + 1)

    # return n_clusters


    # dict_group_2_columns = dict()
    # for index, row in coldata.iterrows():
    #     if row[reference_column] not in dict_group_2_columns.keys():
    #         dict_group_2_columns[row[reference_column]] = [index]
    #     else:
    #         dict_group_2_columns[row[reference_column]].append(index)
    #
    # """
    # join df_kmeans and df_counts
    # """
    # df_kmeans = df_kmeans["cluster"]
    #
    # df = df_counts.join(df_kmeans, how="inner")
    #
    #
    # reference_columns = dict_group_2_columns[reference_entry]

    # return dict_group_2_columns
    #
    for i in range(0,n_clusters):
        df2 = df[df.loc[:,"cluster"] == i]
        # return df2
        # df2 = df2[columns]
    #
    #     data= {"group":[], "count":[], "timepoint":[]}
    #
    #     for index, row in df2.iterrows():
    #         # counter += 1
    #
    #         ref_count_avg = np.mean(row[reference_columns])
    #
    #         data["group"].append(coldata.loc[reference_columns[0], color_column])
    #         data["count"].append(0)
    #         data["timepoint"].append(coldata.loc[reference_columns[0], x_axis_group_column])
    #
    #
    #         for group in dict_group_2_columns.keys():
    #             # print(group)
    #
    #             if reference_entry == group:
    #                 continue
    #             else:
    #                 columns_toUse = dict_group_2_columns[group]
    #                 # print(columns_toUse)
    #                 counts = np.mean(row[columns_toUse])
    #             #
    #                 data["group"].append(coldata.loc[columns_toUse[0], color_column])
    #                 data["count"].append(int(counts - ref_count_avg))
    #                 data["timepoint"].append(coldata.loc[columns_toUse[0], x_axis_group_column])
    # #
    #
        # return data
        ax = sns.pointplot(x="timepoint",y="count",hue="group",data=df2)
        # ax.get_legend().remove()
        if label != "":
            plt.savefig(os.path.join(outputdir, str(n_clusters) + "_subset_countsNormalized_cluster_" + str(i) +"_" +label+".png"))
        else:
            plt.savefig(os.path.join(outputdir, str(n_clusters) + "_subset_countsNormalized_cluster_" + str(i) + ".png"))
        plt.close("all")


def py_visualizeKmeans_Normalized(df, variable, outputdir, label=""):

    n_clusters = df["cluster"].max() + 1

    columns = []
    for column in df.columns:
        if variable.lower() in column.lower():
            columns.append(column)

    for i in range(0,n_clusters):
        df2 = df[df.loc[:,"cluster"] == i]
        df2 = df2[columns]


        data= {"timepoint":[], "LFC":[]}

        # counter = 0
        for index, row in df2.iterrows():
            # counter += 1
            for column in df2.columns:
                data["timepoint"].append(column)
                data["LFC"].append(row[column])
            # if counter > 50:
            #     break
            # print(data)


        ax = sns.pointplot(x="timepoint",y="LFC",data=data)
        # ax.get_legend().remove()
        if label != "":
            plt.savefig(os.path.join(outputdir, str(n_clusters) + "_subset_cluster_" + str(i) +"_"+label+"_.png"))
        else:
            plt.savefig(os.path.join(outputdir, str(n_clusters) + "_subset_LFC_Normalized_cluster_" + str(i) + "_"+label+".png"))

        plt.close("all")
    # n_clusters = df["cluster"].max() + 1
    #
    # columns = []
    # for column in df.columns:
    #     if variable.lower() in column.lower():
    #         columns.append(column)
    #
    # for i in range(0,n_clusters):
    #     df2 = df[df.loc[:,"cluster"] == i]
    #     df2 = df2[columns]
    #
    #
    #
    #     data= {"timepoint":[], "LFC":[]}
    #
    #     # counter = 0
    #     for index, row in df2.iterrows():
    #         # counter += 1
    #         data["timepoint"].append(reference_column)
    #         ref_LFC = row[reference_column]
    #         data["LFC"].append(0)
    #         for column in df2.columns:
    #             if column != reference_column:
    #                 data["timepoint"].append(column)
    #                 data["LFC"].append(int(row[column] - ref_LFC))
    #         # if counter > 50:
    #         #     break
    #         # print(data)
    #
    #
    #     ax = sns.pointplot(x="timepoint",y="LFC",data=data)
    #     # ax.get_legend().remove()
    #     plt.savefig(outputdir+"/"+str(n_clusters) + "_subset_LFC_Normalized_cluster_" + str(i) + ".png")
    #     plt.close("all")


def innerJoinOnGenesAndCluster(df, df_kmeans):

        df_kmeans = df_kmeans["cluster"]
        df = df.join(df_kmeans, how="inner")

        return df


"""
takes LFC and normalizes each row on a reference column.
"""
def py_generateNormalized_LFC(df, variable, reference_column):

    df2 = py_normalizeOnVariable(df, variable, reference_column)
    df2 = innerJoinOnGenesAndCluster(df2, df)

    return df2







"""
returns data frame normalized on count reference column.

This dataframe has a weird format but it makes it easier to plot needed figures.
"""
def py_generateNormalized_Count(df_kmeans, df_counts, coldata, reference_column, reference_entry, color_column, x_axis_group_column, outputdir):

    n_clusters = df_kmeans["cluster"].max() + 1


    dict_group_2_columns = dict()
    for index, row in coldata.iterrows():
        if row[reference_column] not in dict_group_2_columns.keys():
            dict_group_2_columns[row[reference_column]] = [index]
        else:
            dict_group_2_columns[row[reference_column]].append(index)

    """
    join df_kmeans and df_counts
    """
    df_kmeans = df_kmeans["cluster"]

    df = df_counts.join(df_kmeans, how="inner")


    reference_columns = dict_group_2_columns[reference_entry]

    # return dict_group_2_columns
    #

        # df2 = df2[columns]

    data= {"group":[], "count":[], "timepoint":[], "gene":[], "cluster":[]}

    for index, row in df.iterrows():
        # counter += 1

        ref_count_avg = np.mean(row[reference_columns])

        data["group"].append(coldata.loc[reference_columns[0], color_column])
        # data["count"].append(0)
        data["count"].append(1)
        data["timepoint"].append(coldata.loc[reference_columns[0], x_axis_group_column])
        data["gene"].append(index)
        data["cluster"].append(row["cluster"])


        for group in dict_group_2_columns.keys():
            # print(group)

            if reference_entry == group:
                continue
            else:
                columns_toUse = dict_group_2_columns[group]
                # print(columns_toUse)
                counts = np.mean(row[columns_toUse])
            #
                data["group"].append(coldata.loc[columns_toUse[0], color_column])
                # data["count"].append(int((counts)-(ref_count_avg)))
                data["count"].append(float((counts+.01)/(ref_count_avg+.01)))
                data["timepoint"].append(coldata.loc[columns_toUse[0], x_axis_group_column])
                data["gene"].append(index)
                data["cluster"].append(row["cluster"])

    #
    #
        # return data

    df3 = pd.DataFrame.from_dict(data)

    return df3

        # df3.to_csv(outputdir+"/"+str(n_clusters) + "_subset_countsNormalized_Data.csv")
        # ax = sns.pointplot(x="timepoint",y="count",hue="group",data=df3)
        # # ax.get_legend().remove()
        # plt.savefig(outputdir+"/"+str(n_clusters) + "_subset_countsNormalized_cluster.png")
        # plt.close("all")

def py_getMaxRange(df_kmeans):

    n_clusters = df_kmeans["cluster"].max() + 1

    max_range = 0
    for i in range(0, n_clusters):
        df_temp = df_kmeans[df_kmeans.loc[:,"cluster"] == i]

        max = 0
        min = 0
        for column in df_temp.columns:

            if column != "cluster":

                temp_max = df_temp[column].max()
                if temp_max > max:
                    max = temp_max

                temp_min = df_temp[column].min()
                if temp_min < min:
                    min = temp_min

        temp_max_range = max-min

        if temp_max_range > max_range:
            max_range = temp_max_range

    return np.ceil(max_range)

def py_getMaxRangeOfCluster(df_kmeans, cluster):


    df_temp = df_kmeans[df_kmeans.loc[:,"cluster"] == cluster]

    max_range = 0
    max = 0
    min = 0
    for column in df_temp.columns:

        if column != "cluster":

            temp_max = df_temp[column].max()
            if temp_max > max:
                max = temp_max

            temp_min = df_temp[column].min()
            if temp_min < min:
                min = temp_min

    temp_max_range = max-min

    if temp_max_range > max_range:
        max_range = temp_max_range

    return [np.ceil(max_range), max, min]

"""
generates a combined plot to vizualize the kmeans results.
"""
def py_kmeansGraphsCombined(df_kmeans_lfc, df_kmeans_norm, df_counts, outputdir, estimator="median", label=""):


    max_range = py_getMaxRange(df_kmeans_lfc)


    n_clusters = df_kmeans_lfc["cluster"].max() + 1

    columns_lfc = []
    for column in df_kmeans_lfc.columns:
        if "lfc" in column.lower():
            columns_lfc.append(column)

    columns_norm = []
    for column in df_kmeans_norm.columns:
        if "lfc" in column.lower():
            columns_norm.append(column)

    for i in range(0,n_clusters):

        """
        making kmeans lfc data
        """
        df2 = df_kmeans_lfc[df_kmeans_lfc.loc[:,"cluster"] == i]
        df2 = df2[columns_lfc]



        data= {"timepoint":[], "LFC":[], "gene":[]}

        # counter = 0
        for index, row in df2.iterrows():
            # counter += 1
            for column in df2.columns:
                data["timepoint"].append(column)
                data["LFC"].append(row[column])
                data["gene"].append(index)
            # if counter > 50:
            #     break
            # print(data)

        '''
        making kmeans normalized lfc data.
        '''
        df_norm = df_kmeans_norm[df_kmeans_norm.loc[:,"cluster"] == i]
        df_norm = df_norm[columns_norm]


        data_norm = {"timepoint":[], "LFC":[]}

        # counter = 0
        for index, row in df_norm.iterrows():
            # counter += 1
            for column in df_norm.columns:
                data_norm["timepoint"].append(column)
                data_norm["LFC"].append(row[column])


        data_count = df_counts[df_counts.loc[:,"cluster"] == i]






        x = 5
        fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(x*4.67, x))

        ax1, ax2, ax3, ax4 = ax.flatten()

        ax = sns.pointplot(x="timepoint",y="LFC",hue="gene",data=data,ax=ax1)
        ax.get_legend().remove()
        ax1.set_title("LFC")

        # ax = sns.pointplot(x="timepoint",y="LFC",data=data_norm,ax=ax2,estimator=np.median)

        ax2.set_title("LFC Median")
        ax3.set_title("LFC Median Scaled")

        cluster_MaxRange = py_getMaxRangeOfCluster(df_kmeans_lfc, i)
        scale = (max_range - cluster_MaxRange[0])/2
        ax3.set_ylim([cluster_MaxRange[2]-scale, cluster_MaxRange[1]+scale])

        # ax = sns.pointplot(x="timepoint",y="count",hue="group",data=data_count, ax=ax3,estimator=np.median)
        ax4.set_title("Counts Normalized")

        if estimator == "median":
            ax = sns.pointplot(x="timepoint",y="LFC",data=data_norm,ax=ax2,estimator=np.median)
            ax = sns.pointplot(x="timepoint",y="LFC",data=data_norm,ax=ax3,estimator=np.median)
            ax = sns.pointplot(x="timepoint",y="count",hue="group",data=data_count, ax=ax4,estimator=np.median)
        elif estimator == "mean":
            ax = sns.pointplot(x="timepoint",y="LFC",data=data_norm,ax=ax2,estimator=np.mean)
            ax = sns.pointplot(x="timepoint",y="LFC",data=data_norm,ax=ax3,estimator=np.mean)
            ax = sns.pointplot(x="timepoint",y="count",hue="group",data=data_count, ax=ax4,estimator=np.mean)


        fig.suptitle("Cluster " + str(i))


        # plt.tight_layout()
        if label != "":
            plt.savefig(os.path.join(outputdir, "kmeans_all_results_"+str(i)+"_"+label+".png"), dpi=250)
        else:
            plt.savefig(os.path.join(outputdir, "kmeans_all_results_"+str(i)+".png"), dpi=250)

        plt.close("all")

def py_enrichmentVisualization(df, outputdir, filename, layout=None, bottom=.5, left=.5):

    labels = []
    for column in df.columns:
        name_list = column.split("_")
        comparison = name_list[1]
        reg = name_list[2]
        labels.append(f'{comparison}_{reg}')

    cmap = matplotlib.colors.ListedColormap(['white', '#f51000'])

    if 4 >= len(labels):
        x_axis_size = 4*2
    else:
        x_axis_size = len(labels)*2


    fig, ax = plt.subplots(figsize=(x_axis_size,(.2*len(df.index))))
    plt.pcolor(df, cmap= cmap, edgecolor="black", linewidth=.5)
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), labels)
    plt.xticks(rotation=90)
    plt.margins(.9, tight=None)
    # 
    if layout == "tight":
        plt.tight_layout()
    else:
        plt.subplots_adjust(left=left, right=.99, top=.99, bottom=bottom)
    plt.savefig(os.path.join(outputdir, filename+"test.png"), dpi=250)
    # plt.savefig(filename+".png", dpi=250)

    plt.close("all")



def py_subsetEnrichmentResults(df, percentile=None, regulation=None, columns=None):

    if not columns:
        columns = []
        is_per = True
        id_reg = True

        for column in df.columns:
            if not percentile or (percentile == float(column.split("_")[0])):
                is_per = True
            else:
                is_per = False

            if not regulation or (regulation.lower() in column.lower()):
                is_reg = True
            else:
                is_reg = False

            if is_reg and is_per:
                columns.append(column)

    df = df[columns]


    index_to_keep = []
    for index, row in df.iterrows():
        keep = False
        for column in df.columns:
            keep = keep or row[column]

        if keep:
            index_to_keep.append(index)
        else:
            pass

    df = df.loc[index_to_keep,:]

    return df


def py_intersectOfColumns_GoTerms(df, percentile=None, regulation=None, columns=None):

    subset = py_subsetEnrichmentResults(df, percentile, regulation, columns)

    firstEntry=True

    for column in subset.columns:
        if firstEntry:
            goterms = set(subset.loc[subset[column]].index)
            firstEntry = False
        else:
            goterms = goterms.intersection(set(subset.loc[subset[column]].index))

        # print(goterms)

    return goterms



def py_venDiagramGoTerms(df_True_False, outputdir, filename):

    dict_Holder = {}
    for column in df_True_False.columns:
        dict_Holder[column] = set(df_True_False.loc[df_True_False[column]].index)

    venn(dict_Holder)
    plt.savefig(os.path.join(outputdir, filename+".png"), dpi=250)
    plt.close("all")

def py_venn(list_sets, list_labels, outputdir, filename):

    if len(list_sets) != len(list_labels):
        return("error, labels must equal lists")

    dict_Holder = {}
    for i in range(0,len(list_labels)):
        dict_Holder[list_labels[i]] = set(list_sets[i])

    venn(dict_Holder)
    plt.savefig(os.path.join(outputdir, filename+".png"), dpi=250)
    plt.close("all")


def binaryGroupGenerator(n_sets):
    for i in range(1, 2**n_sets):
        yield bin(i)[2:].zfill(n_sets)

def get_all_set_combinations(list_sets, list_labels):

    result_dict = {}
    n_sets = len(list_sets)

    if len(list_sets) != len(list_labels):
        return("error, labels must equal lists")

    dict_Holder = {}
    for i in range(0,len(list_labels)):
        dict_Holder[list_labels[i]] = set(list_sets[i])
    
    set_list = list(dict_Holder.values())

    set_labels = list(dict_Holder.keys())
    
    for logic in binaryGroupGenerator(n_sets):

        included_sets = []
        included_labels = []
        excluded_sets = []
        excluded_labels = []

        for i in range(n_sets):
            if logic[i] == "1":
                included_sets.append(set_list[i])

                included_labels.append(set_labels[i])
            elif logic[i] == "0":
                excluded_sets.append(set_list[i])

                excluded_labels.append(set_labels[i])
                
        label = " ".join(included_labels)
        
        result_dict[label] = set.intersection(*included_sets) - set.union(set(), *excluded_sets)

    
    return result_dict
    





def py_exclusiveSet(df_True_False, column_include):

    if isinstance(column_include, str):
        goterms = set(df_True_False.loc[df_True_False[column_include]].index)

        for column in df_True_False.columns:
            if column != column_include:
                terms = set(df_True_False.loc[df_True_False[column]].index)
                goterms = goterms - terms
    else:

        firstEntry = True
        for column in column_include:
            if firstEntry:
                goterms = set(df_True_False.loc[df_True_False[column]].index)
                firstEntry = False
            else:
                goterms = goterms.intersection(set(df_True_False.loc[df_True_False[column]].index))


        for column in df_True_False.columns:
            if column not in column_include:
                terms = set(df_True_False.loc[df_True_False[column]].index)
                goterms = goterms - terms

    return goterms


# this was just being used to see why the results are the way they are on the normilizaiton.
def py_exploringTheNormalization(df_kmeans, df_counts, coldata, reference_column, reference_entry):


    dict_group_2_columns = dict()
    for index, row in coldata.iterrows():
        if row[reference_column] not in dict_group_2_columns.keys():
            dict_group_2_columns[row[reference_column]] = [index]
        else:
            dict_group_2_columns[row[reference_column]].append(index)


    df_kmeans = df_kmeans["cluster"]

    df = df_counts.join(df_kmeans, how="inner")


    reference_columns = dict_group_2_columns[reference_entry]


    data = dict()

    data["cluster"] = []
    for key in dict_group_2_columns.keys():
        data[key] = []
        data[key+"_normMinus"] = []
        data[key+"_normDivide"] = []


    for index, row in df.iterrows():
        # counter += 1

        ref_count_avg = np.mean(row[reference_columns])

        data[reference_entry].append(ref_count_avg)
        data[reference_entry+"_normMinus"].append(0)
        data[reference_entry+"_normDivide"].append(1)
        data["cluster"].append(row["cluster"])


        for group in dict_group_2_columns.keys():
            # print(group)

            if reference_entry == group:
                continue
            else:
                columns_toUse = dict_group_2_columns[group]
                # print(columns_toUse)
                counts = np.mean(row[columns_toUse])
            #
                data[group].append(counts)
                data[group+"_normMinus"].append(int((counts)-(ref_count_avg)))
                data[group+"_normDivide"].append(float((counts)/(ref_count_avg)))

    #
    #
        # return data

    df3 = pd.DataFrame.from_dict(data)

    return df3


def py_lookingAtNormalizaiton(df, cluster):

    df = df[df.loc[:,"cluster"] == cluster]

    return df.describe()




def py_histogramForColumnsAndCluster(df, cluster):

    df = df[df.loc[:,"cluster"] == cluster]


    for column in df.columns:
        holder = df[column]

        largest = holder[holder != np.Inf].max()
        holder = holder.replace(np.inf, largest)

        sns.distplot((holder), hist=True, kde=False,
                 color = 'darkblue',
                 hist_kws={'edgecolor':'black'},
                 kde_kws={'linewidth': 4})

        plt.title(column+" "+str(cluster))

        plt.savefig("histogram_"+column+"_"+str(cluster)+".png")
        plt.close("all")



def py_subsetColumnOnThreshold(df, column=None, percentile=None, threshold=None, direction="above"):

    if percentile:
        percentile = percentile*100
        threshold = np.nanpercentile(df[column], percentile)
    elif not threshold:
        return "error"

    if direction == "above":
        keep = df[df[column] > threshold]
    elif direction == "below":
        keep = df[df[column] < threshold]

    return keep.index

def py_Difference(set1, set2):

    set1 = set(set1)
    set2 = set(set2)

    return set1 - set2

def py_Intersect(set1, set2):
    set1 = set(set1)
    set2 = set(set2)

    return set1.intersection(set2)

def py_Union(set1, set2):
    set1 = set(set1)
    set2 = set(set2)

    return set1.union(set2)

def py_plotCounts(cts, coldata, referenceColumn, gene, outputdir, groups=None, label=""):

    data = {"group":[], "count":[], "column":[]}

    for index, row in coldata.iterrows():

        group = row[referenceColumn]

        if groups == None:
            data["count"].append(cts.loc[gene,index])
            data["group"].append(group)
            data["column"].append(index)

        else:
            if group in groups:
                data["count"].append(cts.loc[gene,index])
                data["group"].append(group)
                data["column"].append(index)




    x = 5
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(x, x))


    ax = sns.pointplot(x="group", y="count",hue="column", data=data,ax=ax,ci="sd")
    ax.get_legend().remove()
    ax.set_title("count plot " + gene)

    plt.savefig(os.path.join(outputdir, "countplot_"+gene+"_"+label+".png"))
    plt.close("all")


def py_subsetDifferentialExpressionResults(df, percentile, regulation, columns = None, variable = None):


    if regulation == "up":
        xor_bool = False
    else:
        xor_bool = True



    percentile = percentile*100

    index_keep = set()
    genesToKeep = []


    if columns == None and variable == None:
        return ("Error, need column or variable")
    elif columns != None and variable != None:
        return ("Error either column or variable not both")
    elif variable != None:
        columns = []
        for column in df.columns:
            if variable.lower() in column.lower():
                columns.append(column)


    if type(columns) == type(""):
        # pass
        # print(columns)
        # print(type(columns))
        lfc_column = "LFC_" + columns.split("_")[1]
        keep = df[df[columns] > np.nanpercentile(df[columns], percentile)]
        #
        for index, row in keep.iterrows():
            if (row[lfc_column] > 0) != xor_bool:
                # print(row[columns] > 0)
                genesToKeep.append(index)
        #
        index_keep = index_keep.union(set(genesToKeep))
    else:
        for column in columns:
            keep = df[df[column] > np.nanpercentile(df[column], percentile)]

            lfc_column = "LFC_" + column.split("_")[1]

            for index, row in keep.iterrows():
                if (row[lfc_column] > 0) != xor_bool:
                    genesToKeep.append(index)

            index_keep = index_keep.union(set(genesToKeep))

    return(df.loc[index_keep,columns])


def py_intersectOfColumns_SigGenes(df, percentile, regulation, columns = None, variable = None):

    # subset = py_subsetDifferentialExpressionResults(df=df, percentile=percentile, regulation=regulation, columns=columns, variable=variable)
    #
    # firstEntry = True
    #
    # percentile = percentile*100
    # for column in subset.columns:
    #
    #     keep = subset[subset[column] > np.nanpercentile(subset[column], percentile)]
    #     if firstEntry:
    #         genes = set(keep.index)
    #         firstEntry = False
    #     else:
    #         genes = genes.intersection(set(keep.index))


    if columns == None and variable == None:
        return ("Error, need column or variable")
    elif columns != None and variable != None:
        return ("Error either column or variable not both")
    elif variable != None:
        columns = []
        for column in df.columns:
            if variable.lower() in column.lower():
                columns.append(column)

    firstEntry = True
    for column in columns:

        genes = py_subsetDifferentialExpressionResults(df=df, columns=column, percentile=percentile, regulation=regulation)

        if firstEntry:
            gene_intersect = set(genes.index)
            firstEntry = False
        else:
            gene_intersect = gene_intersect.intersection(set(genes.index))

        # print(goterms)

    return list(gene_intersect)


def py_geneVennDiagram(df, columns, percentile, regulation, outputdir, filename):

    dict_Holder = {}
    for column in columns:

        genes = py_subsetDifferentialExpressionResults(df=df, columns=column, percentile=percentile, regulation=regulation)
        print(genes)
        dict_Holder[column] = set(genes.index)

    venn(dict_Holder)
    plt.savefig(os.path.join(outputdir, filename+".png"), dpi=250)
    plt.close("all")

# def py_test(df, variable, threshold=None, percentile=None):
#     df = df.dropna()
#     n = len(df.index)
#     values = list(df[variable].values)
#     values.sort()
#
#     if (threshold is None) and (percentile is None):
#         return("Error")
#     elif threshold:
#         pass
#     elif percentile:
#         per_Index = int(np.ceil(n*percentile))
#         threshold = values[per_Index]
#
#     print(threshold)
#
#     return([(np.nanpercentile(values, 90)), (np.percentile(values, 90))])

def py_getlistofGenes_Padj_range(res, max_=np.inf, min_=-np.inf):
    holder = res[(-np.log(res["padj"]) > min_) & (-np.log(res["padj"]) < max_)]


    return list(holder.index)


"""
Takes a dataframe with c merolae genes and returns a dataframe with the columns attached.
"""
def py_addProteinColumnToDataframe(df):
    proteins = genesToProteindf[["Protein", "Function"]].astype(str)
    return proteins.join(df, how="right")
