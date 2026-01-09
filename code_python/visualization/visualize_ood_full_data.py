import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from visualization_setting import set_plot_style, save_img_path
from typing import Callable, Optional, Union, Dict, Tuple
# from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score
from scipy.stats import wasserstein_distance_nd,spearmanr
from sklearn.preprocessing import StandardScaler
import json
from matplotlib.lines import Line2D
from visualize_ood_learning_curve import (
                                        get_score, 
                                        ensure_long_path,
                                        process_learning_curve_scores,
                                        plot_bar_ood_iid,
                                        get_comparison_of_features
                                        )

set_plot_style(tick_size=16)
MARKERS = ['o', '^', 's', 'D', 'v', '<', '>', 'p', '*', 'X', 'h']

HERE: Path = Path(__file__).resolve().parent
DATASETS: Path = HERE.parent.parent / 'datasets'
training_df_dir: Path = DATASETS/ 'training_dataset'
RG_DATA = pd.read_pickle(training_df_dir/'Rg data with clusters aging imputed.pkl')
RG_DATA['Polymers cluster'] = RG_DATA['canonical_name'] 
target = 'OOD_target_log Rg (nm)'
results_path = HERE.parent.parent / 'results'/ target


#---------Covariance Features---------
covarience_features = ['Xn', 'Mw (g/mol)', 'PDI', 'Concentration (mg/ml)', 
                        'Temperature SANS/SLS/DLS/SEC (K)',
                        'polymer dP', 'polymer dD' , 'polymer dH',
                        'solvent dP', 'solvent dD', 'solvent dH',
                        "Dark/light", "Aging time (hour)", "To Aging Temperature (K)",
                        "Sonication/Stirring/heating Temperature (K)", "Merged Stirring /sonication/heating time(min)"
                        ]

sd_caler = StandardScaler()
cov_continuous_vector = RG_DATA[covarience_features].to_numpy(dtype=float)
cov_continuous_vector_scaled = sd_caler.fit_transform(cov_continuous_vector)
cov_mordred_vector = pd.DataFrame(RG_DATA['Trimer_Mordred'].tolist()).to_numpy(dtype=float)
cov_mordred_vector_scaled = sd_caler.fit_transform(cov_mordred_vector)
cov_ECFP_vector = np.array(RG_DATA['Trimer_ECFP6_count_512bits'].tolist(), dtype=int)
cv_MACCS_vector = np.array(RG_DATA['Trimer_MACCS'].tolist())
cov_mordred_and_continuous_vector = np.concatenate([cov_mordred_vector, cov_continuous_vector], axis=1)
cov_mordred_and_continuous_vector_scaled = sd_caler.fit_transform(cov_mordred_and_continuous_vector)

naming: dict = {
        'mordred vector': cov_mordred_vector_scaled,
        'ECFP vector': cov_ECFP_vector,
        'all features vector': cov_continuous_vector_scaled,
        'combined mordred-numerical vector': cov_mordred_and_continuous_vector_scaled,
        'MACCS vector': cv_MACCS_vector,
}


#----------- Helper Functions --------------
def get_ood_iid(scores,ml_score_metric: str,algorithm: str):
    data = []
    for cluster, _ in scores.items():
        if cluster.startswith("CO_") or cluster.startswith("ID_") or cluster.startswith("IID_"):
            # print(cluster)
            # cluster_id = cluster.split("_")[1]
            mean_ood, std_ood = get_score(scores, f"{cluster}", ml_score_metric)
            data.append({
            "Cluster": cluster,
            "Model": algorithm,
            f"Score": mean_ood,
            f"Std": std_ood,
        })
    return data


def weighted_jaccard(u, v,eps: float = 1e-6):
    min_sum = np.sum(np.minimum(u, v))
    max_sum = np.sum(np.maximum(u, v))
    return 1 - ((min_sum+eps) / (max_sum+eps))


def calculate_wasserstein_distance(feature_space:np.ndarray, 
                                   cluster_col:str,
                                   data:pd.DataFrame
                                   ):
    all_cat =  data[cluster_col].to_numpy()
    unique_all_clusters, _ = np.unique(all_cat, return_counts=True)
    wasserstein_scores = {}
    if cluster_col =='substructure cluster':
        side_chain_cat = data['Side Chain Cluster'].to_numpy()
        unique_side_clusters, label_counts = np.unique(side_chain_cat, return_counts=True)
        for test_label in unique_all_clusters:
            test_mask = all_cat == test_label
            train_mask = ~test_mask
            distance_scores= wasserstein_distance_nd(feature_space[test_mask],feature_space[train_mask])
            wasserstein_scores[test_label] = distance_scores
        for test_label in unique_side_clusters:
            test_mask = side_chain_cat == test_label
            train_mask = ~test_mask
            distance_scores= wasserstein_distance_nd(feature_space[test_mask],feature_space[train_mask])
            wasserstein_scores['Polar'] = distance_scores
            break

    else:
        for test_label in unique_all_clusters:
            test_mask = all_cat == test_label
            train_mask = ~test_mask
            distance_scores= wasserstein_distance_nd(feature_space[test_mask],feature_space[train_mask])
            wasserstein_scores[test_label] = distance_scores

    return wasserstein_scores


def make_accumulating_scores(scores, ml_score_metric: str, 
                             co_vector, cluster_types: str,
                            algorithm: str, is_equal_size:bool=False,
                            is_ood_iid_distance:bool=True) -> list:


    data = []
    wasserstein_scores = calculate_wasserstein_distance(feature_space=naming[co_vector],
                                                        cluster_col=cluster_types,
                                                        data=RG_DATA)
    
    all_clusters = RG_DATA[cluster_types].unique().tolist()

    all_clusters = [
        cluster_id for cluster_id in all_clusters 
        if (RG_DATA[cluster_types] == cluster_id).sum() >= 2
    ]
    if cluster_types == 'substructure cluster':
        all_clusters.append('Polar')

    for cluster_id in all_clusters:

        if is_equal_size==True:
            _, scores_at_equal_training_size = process_learning_curve_scores(scores,ml_score_metric)
            cluster_scores = scores_at_equal_training_size[scores_at_equal_training_size['Cluster'] == f"CO_{cluster_id}"]
            cluster_test_scores = cluster_scores[cluster_scores["Score Type"] == "Test"]
            mean = cluster_test_scores["Score"].iloc[0]
            std = cluster_test_scores["Std"].iloc[0]
        else:
            if is_ood_iid_distance==True:
                mean_ood, std_ood = get_score(scores, f"CO_{cluster_id}", ml_score_metric)
                mean_ID, std_ID = get_score(scores, f"ID_{cluster_id}", ml_score_metric)
                mean = abs(mean_ood - mean_ID)/abs(mean_ID)*100
                std = 0
            else:
                mean_ood, std_ood = get_score(scores, f"CO_{cluster_id}", ml_score_metric)
                mean = mean_ood
                std = std_ood
        data.append({
            "Cluster": cluster_id,
            "Model": algorithm,
            "wasserstein distance": wasserstein_scores.get(cluster_id, np.nan),
            f"{ml_score_metric}_mean": mean,
            f"{ml_score_metric}_std": std
        })
    return data


#----------- Plotting Function --------------

def plot_OOD_Score_vs_distance(df, ml_score_metric: str,
                                score_metric: str = 'rmse',
                                saving_path: Optional[Path] = None,
                                file_name: str = None,
                                is_equal_size: bool = False,
                                is_ood_iid_distance: bool = True,
                                figsize=(6, 5),
                                fontsize: int = 14,  
                                ) -> None:

    models = list({row['Model'] for row in df})
    clusters = list({row['Cluster'] for row in df})
    print(df)

    palette = sns.color_palette("tab20", len(clusters)) if len(clusters) > 8 else sns.color_palette("Set2", len(clusters))
    color_map = {cluster: palette[i] for i, cluster in enumerate(sorted(clusters))}
    marker_map = {model: MARKERS[i % len(MARKERS)] for i, model in enumerate(sorted(models))}

    legend_elements = [
        Line2D([0], [0], marker='_', color=color_map[cluster], linestyle='None',
               label=f"Cluster {cluster}", markersize=20, markeredgewidth=5)
        for cluster in sorted(clusters)
    ] + [
        Line2D([0], [0],
               marker=marker_map[model],
               linestyle='None',
               label=model,
               markersize=10,
               markerfacecolor='none' if marker_map[model] == 'o' else 'k',
               markeredgecolor='k')
        for model in sorted(models)
    ]

    plt.figure(figsize=figsize)

    for row in df:
        marker = marker_map[row["Model"]]
        is_circle = marker == 'o'

        plt.errorbar(
            row['wasserstein distance'],
            row[f"{ml_score_metric}_mean"],
            yerr=row[f"{ml_score_metric}_std"],
            fmt=marker,
            markerfacecolor='none' if is_circle else color_map[row["Cluster"]],
            markeredgecolor=color_map[row["Cluster"]],
            markersize=10,
            linestyle='none'
        )

    if is_equal_size or not is_ood_iid_distance:
        y_label = f"{ml_score_metric.upper()} Score"
    else:
        y_label = f"Performance Degradation % ({score_metric.upper()})"

    y_max = max(row[f"{ml_score_metric}_mean"] for row in df)
    y_max_tick = np.ceil(y_max * 5) / 5
    yticks = np.arange(0,  500, 100)
    plt.yticks(yticks)
    plt.ylabel(y_label, fontsize=fontsize-1.5, fontweight='bold')
    plt.xlabel("Wasserstein Test-Train Distance", fontsize=fontsize-.5, fontweight='bold')
    plt.legend(handles=legend_elements, loc="best", fontsize=fontsize - 3)  # keep legend slightly smaller
    plt.tick_params(axis='both', labelsize=fontsize-1)

    plt.tight_layout(rect=[0, 0, 1, 1.05])
    plt.tight_layout()
    if saving_path and file_name:
        save_img_path(saving_path, f"{file_name}_distance_metric-Wasserstein.png")
    plt.show()
    plt.close()




# comparison_of_features_lc = {
#     '(Mordred-Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_RF_Standard':'Mordred+continuous',
#     '(Mordred-Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_RF_Standard':'Mordred+continuous+aging',
#     '(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_RF_Standard':'continuous',
#     '(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_RF_Standard':'continuous+aging',
# }




if __name__ == "__main__":
    cluster_list = [
                    # 'KM4 ECFP6_Count_512bit cluster',	
                    # 'KM3 Mordred cluster',
                    # 'HBD3 MACCS cluster',
                    # 'KM5 polymer_solvent HSP and polysize cluster',
                    # 'KM4 polymer_solvent HSP and polysize cluster',
                    # 'substructure cluster',
                    # 'KM4 polymer_solvent HSP cluster',
                    # 'KM4 Mordred_Polysize cluster',
                    'Polymers cluster',
                    ]
    for cluster in cluster_list:
        score_metrics = ["rmse"]
        for accuracy_metric in score_metrics:
            ## Plot OOD-IID vs distance for mix of fingerprints and numerical

            # for fp in ['MACCS', 'Mordred', 'ECFP3.count.512']:
            #     co_vectors = ['numerical vector']
            #     if fp == 'Mordred':
            #         co_vectors.extend(['mordred vector', 'combined mordred-numerical vector'])
            #     if fp == 'MACCS':
            #         co_vectors.append('MACCS vector')

            #     if fp == 'ECFP3.count.512':
            #         co_vectors.append('ECFP vector')
                
            #     for co_vector in co_vectors:
            #         combined_data = []
            #         for model in ['XGBR', 'NGB', 'RF']:
            #                 scores_folder_path = results_path / cluster / 'Trimer_scaler'
            #                 score_file = scores_folder_path / f'({fp}-Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_{model}_Standard_scores.json'
            #                 score_file = ensure_long_path(score_file)
            #                 if not os.path.exists(score_file):
            #                     print(f"File not found: {score_file}")
            #                     continue 

            #                 with open(score_file, "r") as f:
            #                     scores = json.load(f)

            #                 model_data = make_accumulating_scores(scores, accuracy_metric, co_vector, cluster, model,is_ood_iid_distance=False)
            #                 combined_data.extend(model_data)
                                    
            #         saving_folder = scores_folder_path/ f'scores vs distance (full data)'/ f"{co_vector}"
            #         plot_OOD_Score_vs_distance(combined_data, accuracy_metric, co_vector=co_vector,
            #                         saving_path=saving_folder, file_name=f"fingerprint-{fp}_metric-{accuracy_metric}",is_ood_iid_distance=False)
            #         print('scores vs distance (full data)')
                            
                    #  Plot OOD-IID vs distance for equal training size
                    #         score_file = scores_folder_path / f'({fp}-Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH)_{model}_hypOFF_Standard_lc_scores.json'
                    #         score_file = ensure_long_path(score_file)
                    #         if not os.path.exists(score_file):
                    #             print(f"File not found: {score_file}")
                    #             continue 

                    #         with open(score_file, "r") as f:
                    #             scores = json.load(f)

                    #         model_data = make_accumulating_scores(scores, accuracy_metric, co_vector, cluster, model,is_equal_size=True)
                    #         combined_data.extend(model_data)


                    # saving_folder = scores_folder_path/ f'scores vs distance at equal training size'/ f"{co_vector}"
                    # plot_OOD_Score_vs_distance(combined_data, accuracy_metric, co_vector=co_vector,
                    #                 saving_path=saving_folder, file_name=f"fingerprint-{fp}_metric-{accuracy_metric}", is_equal_size=True)
                    # print('Plot equal size scores')




            # Plot OOD-IID vs distance for only numerical
            co_vector = 'all features vector'
            
            combined_data = []
            ood_iid_bar_combined_models = []
            OOD_IID_distance = True
            ncol = 0
            for model in ['RF', 'XGBR']:
                scores_folder_path = results_path / cluster / 'scaler'
                score_file = scores_folder_path / f'(Xn-Mw-PDI-concentration-temperature-polymer dP-polymer dD-polymer dH-solvent dP-solvent dD-solvent dH-light exposure-aging time-aging temperature-prep temperature-prep time)_{model}_Standard_scores.json'
                score_file = ensure_long_path(score_file)
                with open(score_file, "r") as f:
                    scores = json.load(f)
                    
                model_data = make_accumulating_scores(scores, accuracy_metric, co_vector, cluster, model,is_ood_iid_distance=OOD_IID_distance)
                combined_data.extend(model_data)

                ood_iid_bar = get_ood_iid(scores, accuracy_metric, model)
                ood_iid_bar_combined_models.extend(ood_iid_bar)
                ncol+=1
   

            saving_folder = scores_folder_path/  f'scores vs distance (full data)'/ f"{co_vector}"
            f_name = f"{co_vector}_{accuracy_metric}" 
            f_name =  f"{f_name}_OOD-IID" if OOD_IID_distance else f"{f_name}_OOD"

            plot_OOD_Score_vs_distance(combined_data, accuracy_metric,
                                        saving_path=saving_folder, file_name=f_name,
                                        is_ood_iid_distance=OOD_IID_distance,figsize=(6,5),fontsize=18)
            print('Plot scores vs distance (full data)')


            saving_folder = scores_folder_path/  f'OOD-IID bar plot (full data)'
            plot_bar_ood_iid(pd.DataFrame(ood_iid_bar_combined_models), accuracy_metric,
                            saving_folder,f'metric-{accuracy_metric}', 
                            figsize=(16, 7), text_size=22,ncol=ncol)
 
            # Quantify correlation between distance and performance degradation
            df = pd.DataFrame(combined_data)
            results = {}
            for model in df['Model'].unique():
                sub = df[df['Model'] == model]
                rho, pval = spearmanr(sub['wasserstein distance'], sub['rmse_mean'])
                results[model] = {'rho': rho, 'pval': pval}

            # Show results
            for model, res in results.items():
                print(f"{model}: Spearman rho = {res['rho']:.3f}, p = {res['pval']:.3g}")

            



            ## plot bar plot for OOD-IID for comparative features
            # model = 'RF'
            # comparison_of_features_full= get_comparison_of_features(model, '_Standard')

            # combined_data = []
            # for file, file_discription in comparison_of_features_full.items():

            #     if any(keyword in file for keyword in ['Mordred', 'MACCS', 'ECFP3.count.512']):
            #         scores_folder_path = results_path / cluster / 'Trimer_scaler'
            #     else:
            #         scores_folder_path = results_path / cluster / 'scaler'
                
            #     score_file_lc = ensure_long_path(scores_folder_path / f'{file}_scores.json')
            #     predictions_file_lc = ensure_long_path(scores_folder_path / f'{file}_predictions.json')

            #     if not os.path.exists(score_file_lc) or not os.path.exists(predictions_file_lc):
            #         print(f"File not found: {file}")
            #         continue  

            #     with open(score_file_lc, "r") as f:
            #         scores = json.load(f)
            #     with open(predictions_file_lc, "r") as s:
            #         predictions = json.load(s)

            #     model_data = get_ood_iid(scores, accuracy_metric, file_discription)
            #     combined_data.extend(model_data)
            # combined_data = pd.DataFrame(combined_data)
            # saving_folder = results_path / cluster / f'OOD-IID bar plot at full data (comparison of features)'
            # # print(combined_data)
            # plot_bar_ood_iid(combined_data, 'rmse', 
            #                 saving_folder, file_name=f'{accuracy_metric}_{model}_comparison of feature',
            #                     text_size=16, figsize=(13, 8), ncol=3)

            # print("save OOD vs IID bar plot at equal training size for comparison of features")



