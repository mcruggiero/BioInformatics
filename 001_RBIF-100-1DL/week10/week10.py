
# Data Analytics
import numpy as np
import pandas as pd

# Ploting Libraries
import seaborn as sns
import matplotlib.pyplot as plt

# Clustering Libraries
from sklearn.cluster import KMeans, DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

# Confusion matrix
from sklearn.metrics import confusion_matrix


plt.style.use('fivethirtyeight')

class Week10:
    """
    Part 1: load & stats
    Part 2: analysis
    Part 3: cluster
    """
    def __init__(self, clinical_file):
        
        df = self.load(clinical_file)
        cluster_animals = self.analysis(clinical_file, df)
        self.cluster(clinical_file, cluster_animals)
        
    def stats(self, clinical_file, codename):
        """
        Returns mean and standard deviation of codename text file, used in load
        """
        data = clinical_file + "/diversityScores/{}.diversity.txt".format(codename)
        
        # Convert txt files to float list
        text_list = open(data).readlines()
        data_list = [float(x) for x in text_list]
        
        # Calculate statistics
        mean = np.mean(data_list)
        std = np.std(data_list)
    
        return [mean, std] 
    
    def load(self, clinical_file):
        """
        Creates .txt with statistical information and return dataframe
        """
        clinical_data = clinical_file + "clinical_data.txt"
        df = pd.read_csv(clinical_data, sep='\t')
        
        # Creates new columns for statistics
        df["mean"] = df["code_name"].apply(lambda a: self.stats(clinical_file, a)[0])
        df["std"] = df["code_name"].apply(lambda a: self.stats(clinical_file, a)[1])
        
        # Export data frame to csv
        export_name = "clinical_data.stats.txt"
        df.to_csv(export_name, sep='\t')
        
        return df
        
    def analysis(self, clinical_file, df):
        """
        Finds and plots first, second, and last place animal information
        """
        # Sort list and grab values
        sorted_mean_list =  sorted(df["mean"].tolist())
        first_place  = sorted_mean_list[-1]
        second_place = sorted_mean_list[-2]
        last_place   = sorted_mean_list[0]

        animal_dataframes = []
        i = 0
        for place in [first_place, second_place, last_place]:
            # Grab animal name from dataframe
            animal = df[df["mean"] == place]["code_name"].values[0]

            # Convert txt files to float list
            animal_text = clinical_file + "/distanceFiles/{}.distance.txt".format(animal)
            animal_data = pd.read_csv(animal_text, header=None).rename(columns={0: 'x1', 1: 'x2'})
            
            # Plot Values
            animal_data.plot(kind='scatter', x='x1', y='x2', title = animal.title())
            plt.savefig( '{}.png'.format(clinical_file, animal), dpi=300)
            
            # Make dataframe for merger with target labels
            animal_data["target"] = animal
            animal_data["num_target"] = i
            i += 1
            animal_dataframes.append(animal_data)
        
        # Concatenate dataframes and reset index
        cluster_animals = pd.concat(animal_dataframes).reset_index(drop=True)
        return cluster_animals
        
    def cluster(self, clinical_file, cluster_animals):
        """
        Apply K-means clustering to graphs
        """
        # Will not apply train-test split
        X = cluster_animals.loc[:, ['x1', 'x2']]
        y = cluster_animals["num_target"]
        
        # Scale X to effectively cluster
        sc = StandardScaler()
        X_sc = sc.fit_transform(X)

        # k-means clustering test: see notes in jupyter notebook
        scores = []
        for k in range(2, 7):
            cl = KMeans(n_clusters=k)
            cl.fit(X_sc)
            inertia = cl.inertia_
            sil = silhouette_score(X_sc, cl.labels_)
            scores.append([k, inertia, sil])

        score_df = pd.DataFrame(scores)
        score_df.columns = ['k', 'inertia', 'silhouette']
        
        # Plot analysis for best cluster
        fig, axes = plt.subplots(1, 2, figsize=(14, 7));
        fig.suptitle('Test for best clustering of animal data', fontsize=16)
        axes[0].plot(score_df.k, score_df.inertia);
        axes[0].set_title('Inertia over k');
        axes[1].plot(score_df.k, score_df.silhouette);
        axes[1].set_title('Silhouette Score over k');
        plt.savefig('clustering_choice.png', dpi=300);
        
        # Grab best cluster choice from silhoette
        best_silhoette = score_df["silhouette"].max()
        clusters = score_df[score_df["silhouette"] == best_silhoette]["k"].values[0]
        
        km = KMeans(n_clusters=clusters)
        km.fit(X_sc)
        cluster_animals['prediction'] = km.labels_
        
        # Create a dataframe for cluster_centers (centroids)
        centroids = pd.DataFrame(
            sc.inverse_transform(km.cluster_centers_),
            columns=["x1", "x2"])
       
        # Plot scatter by cluster / color, and centroids
        colors = ["red", "green", "blue"]
        cluster_animals['color'] = cluster_animals['prediction'].map(lambda p: colors[p])

        ax = cluster_animals.plot(    
            kind="scatter", 
            x="x1", y="x2",
            figsize=(10,8),
            c = cluster_animals['color'])

        centroids.plot(
            kind="scatter", 
            x="x1", y="x2", 
            marker="*", c=["r", "g", "b"], s=550,
            ax=ax);
        
        plt.savefig('clustering.png', dpi=300)
        
        # Extra-extra credit, confusion matrix
        y_actual = y
        y_pred   = cluster_animals['prediction']
        
        df_confusion = pd.crosstab(y_actual, y_pred, rownames=[''], colnames=['Predicted'], margins=True)
        
        df_confusion.to_csv("confusion_matrix.txt", sep="\t")
        
if __name__== "__main__":
    #### 
    # Important: Match your file path here. Also ensure that your .txts are tab-delimited
    ####
    clinical_file = "./inputfiles/"
    
    Week10(clinical_file)
