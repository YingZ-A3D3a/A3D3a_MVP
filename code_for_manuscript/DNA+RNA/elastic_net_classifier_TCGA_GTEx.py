####################################################################################
# copyright @A3D3a, UT MD Anderson, codes written by Ying Zhu

# These codes identify seed genes for RNAseq data for each cancer type by elastic net classifier
####################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from sklearn.model_selection import train_test_split
from scipy.stats import ranksums
from statsmodels.stats.multitest import fdrcorrection
from imblearn.over_sampling import SMOTE
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
import pickle
import seaborn as sns
from matplotlib.patches import Patch

####################################################################################
def elastic_net_classifier(save_dir,t_types):
    
    ### find differentially expressed genes of training set
    def deg_train(X_all_train,y_train):
        # concatenate training data and labels
        X_all_train_comb = pd.concat([X_all_train,y_train],axis = 1)

        p_values = []
        fold_change = []
        abs_lfc = [] # abolute log2 fold change
        mean_A = []
        mean_B = []

        grpA = X_all_train_comb.loc[X_all_train_comb.group=='primary_tumor',] # primary
        grpA = grpA.iloc[:,:-1]
        grpB = X_all_train_comb.loc[X_all_train_comb.group=='normal',] # normal
        grpB = grpB.iloc[:,:-1]

        for i in range(X_all_train.shape[1]):
            A = 2**(grpA.iloc[:,i].values)
            B = 2**(grpB.iloc[:,i].values)
            p_values.append(ranksums(A,B).pvalue)
            fold_change.append(np.mean(A)/np.mean(B))
            abs_lfc.append(abs(np.log2(np.mean(A)/np.mean(B))))
            mean_A.append(np.mean(A))
            mean_B.append(np.mean(B))

        p_df = pd.DataFrame({'p_value': p_values, 'mean_tumor':mean_A,'mean_normal':mean_B, 'Fold_Change': fold_change,'abs_log2_fc': abs_lfc}, index = X_all_train.columns)
        p_df.sort_values(by = ['p_value','abs_log2_fc'],inplace = True, ascending = [True,False])
        p_df['gene'] = p_df.index

        ### remove inf cases
        p_df = p_df.loc[p_df.abs_log2_fc!=float('inf'),]

        ### compute FDR
        fdrs = fdrcorrection(p_df.p_value.values)[1]
        p_df['FDR_BH'] = fdrs
        p_df = p_df[['gene','mean_tumor','mean_normal','Fold_Change','abs_log2_fc','p_value','FDR_BH']]

        ### save DEGs
        save_fnm= os.path.join(model_dir,'DEG_training_combatseq_TCGA_primary_vs_GTEx+NAT_022123.csv')
        p_df.to_csv(save_fnm,index = False)
        return p_df
    
    def fpr_fnr(y_true, y_pred):    
        fp = np.sum((y_pred == 1) & (y_true == 0))
        tp = np.sum((y_pred == 1) & (y_true == 1))
        fn = np.sum((y_pred == 0) & (y_true == 1))
        tn = np.sum((y_pred == 0) & (y_true == 0))
        fpr = (fp / (fp + tn))
        fnr = (fn / (fn + tp))
        return fpr,fnr
    
    def compute_test_acc(y_test,y_pred,n_scores):
        print(f"predicted y: {Counter(y_pred)}")
        print(f"test accuracy: {(y_pred==y_test).sum()/len(y_test)}")

        ytest = (y_test=='primary_tumor').astype(int)
        ypred = (y_pred=='primary_tumor').astype(int)

        fpr,fnr = fpr_fnr(ytest,ypred)
        val_acc_mean,val_acc_std = np.mean(n_scores),np.std(n_scores)
        test_acc = (y_pred==y_test).sum()/len(y_test)
        tmp = Counter(ytest)
        test_tumor_samples = tmp[1]
        test_normal_samples = tmp[0]

        perf_df = pd.DataFrame({'val_accuracy_mean':[val_acc_mean],
                                'val_accuracy_std':[val_acc_std],
                                'test_tumor_samples':[test_tumor_samples],
                                'test_normal_samples':[test_normal_samples],
                                'test_accuracy':[test_acc],
                                'FPR_test':[fpr],
                                'FNR_test':[fnr]})
        fnm = os.path.join(model_dir,'validation_and_test_accuracy.csv')
        perf_df.to_csv(fnm,index = None)

    ### make heatmap of elastic net selected genes
    def make_heatmap(df_1,g_kp):
        df_select = df_1.loc[:,g_kp]
        df_select['group'] = df.group
        df_select['sample_ID'] = df_select.index
        df_select = df_select.sort_values(by = ['group','sample_ID'])
        data0 = df_select.iloc[:,:-2]
        data0T = data0.T
        data0T.head()
        plt.rcParams["figure.figsize"] = (15,10)

        Group = df_select.group.astype(object)
        lut1 = dict(zip(Group.unique(), ['#60FD00','#ED2323','blue']))
        col_colors1 = Group.map(lut1)

        g = sns.clustermap(data0T, col_cluster = False, method = 'ward', z_score = 0, 
                    cmap = 'vlag',col_colors = col_colors1, vmin=-2, vmax=2.5,
                    yticklabels=True, xticklabels=False, cbar_pos=(.01, 0.75, .03, .08),
                    figsize=(10,(0.25*(len(data0T.index)))))
        # g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 7)
        handles = [Patch(facecolor=lut1[name]) for name in lut1]
        plt.legend(handles, lut1, title='Group',
                bbox_to_anchor=(0.06, 0.9), bbox_transform=plt.gcf().transFigure, loc='center left')
        # plt.tight_layout()
        fnm = os.path.join(model_dir,'heatmap_elastic_net_classifier_genes.png')
        plt.savefig(fnm,dpi = 512)
        plt.show()
        fnm_mat = os.path.join(model_dir,'heatmap_elastic_net_classifier_genes_data_matrix.csv')
        df_select.to_csv(fnm_mat)
    
    for t_type in t_types:
        print(t_type)
        fnm = os.path.join(save_dir,t_type,'primary_normal_TCGA_GTEx_'+t_type+'_combatseq_protein_coding_log2.csv')
        sav_dir = os.path.join(save_dir,t_type)
        df_a = pd.read_csv(fnm,index_col = 0)
        group = df_a.group

            # save the model to disk
        model_dir = os.path.join(sav_dir,'elastic net model')
        os.makedirs(model_dir, exist_ok=True)  # succeeds even if directory exists.
        
        ### select groups for comparions, combine GTEx_normal and TCGA_normal into the same normal group
        grp_cp = ['TCGA_normal','GTEx_normal']
        tmp = []
        for i in df_a.group.values:
            if i in grp_cp:
                tmp.append('normal')
            else:
                tmp.append('primary_tumor')
        df = df_a.copy()
        df.group = tmp
        
        X_all,y = df.iloc[:,:-1],df.iloc[:,-1]

        np.random.seed(12345)
        X_all_train, X_all_test, y_train, y_test = train_test_split(X_all, y,
                                                            stratify=y,
                                                            test_size=0.3)
        
        p_df = deg_train(X_all_train,y_train)
        ### set a stricter threshold
        p_df_subset = p_df.loc[(p_df.FDR_BH<5*10**(-2))&(p_df.abs_log2_fc>1),:]
        p_df_subset = p_df_subset.sort_values(by = ['FDR_BH','abs_log2_fc'],ascending = [True,False])
        p_df_subset.reset_index(inplace = True, drop = True)
        topN = 2500
        p_df_subset = p_df_subset.iloc[:topN,]
        p_df_subset.index = p_df_subset.gene

        degs = p_df_subset.gene.values.tolist()
        X_train = X_all_train.loc[:,degs]
        X_test = X_all_test.loc[:,degs]
        X_train_resampled, y_train_resampled = SMOTE(random_state = 12345).fit_resample(X_train, y_train)
        scaler = MinMaxScaler()
        scaler.fit(X_train_resampled)
        X_train_mat = scaler.transform(X_train_resampled)
        X_test_mat = scaler.transform(X_test)
        X_train_o = pd.DataFrame(X_train_mat, columns = degs)
        X_test_o = pd.DataFrame(X_test_mat, columns = degs)
        y_train_o = y_train_resampled
                
        model_method = 'hinge'
        acc_list1 = []
        x1 = [i for i in np.arange(0.01,1,0.05)]
        for alpha0 in x1:
            EN = SGDClassifier(loss='hinge', penalty='elasticnet', alpha=alpha0, l1_ratio=0.5, random_state = 100)
            cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=1)
            # evaluate the model and collect the scores
            n_scores = cross_val_score(EN, X_train_o, y_train_o, scoring='accuracy', cv=cv, n_jobs=-1)
            acc_list1.append(np.mean(n_scores))

        plt.rcParams["figure.figsize"] = (6,5)
        plt.plot(x1, acc_list1)
        for i in range(1,(len(acc_list1)-1)):
            if (acc_list1[i]-acc_list1[i-1]<0.005) & (acc_list1[i]-acc_list1[i+1]>0.02):
                thre_alpha = x1[i]
                break
        plt.axvline(x = thre_alpha, color = 'green', linestyle='--',label = 'axvline - full height')
        plt.xlabel('alpha')
        plt.ylabel('validation accuracy')
        fnm = os.path.join(model_dir,t_type+'_elastic_net_validation_accuracy.png')
        plt.title('alpha='+str(thre_alpha))
        plt.savefig(fnm,dpi = 512)
        plt.show()

        EN = SGDClassifier(loss=model_method, penalty='elasticnet', alpha=thre_alpha, l1_ratio=0.5, random_state = 100)

        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=1)
        # evaluate the model and collect the scores
        n_scores = cross_val_score(EN, X_train_o, y_train_o, scoring='accuracy', cv=cv, n_jobs=-1)
        # report the model performance
        print('CV Mean Accuracy: %.3f (%.3f)' % (np.mean(n_scores), np.std(n_scores)))

        # predict the class label
        # define a single row of input data
        EN = SGDClassifier(loss=model_method, penalty='elasticnet', alpha=thre_alpha, l1_ratio=0.5, random_state = 100)
        EN.fit(X_train_o, y_train_o)
        y_pred = EN.predict(X_test_o)
        # compute test accuracy, fpr,pnr
        compute_test_acc(y_test,y_pred,n_scores)

        ### save the model
        filename = os.path.join(model_dir,'finalized_model_elastic_net_GTEx_NAT_022123.sav')
        pickle.dump(EN, open(filename, 'wb'))

        # selected features
        genes = X_train.columns
        g_kp = []
        g_kp.extend(genes[EN.coef_[0]!=0])
        g_kp = list(set(g_kp))
        print(f'selected feature numbers: {len(g_kp)}')

        ### save a dataframe of the deg genes with the info of if they are within the classifier genes
        temp = p_df_subset.index
        elastic_net_selected = []
        for i in temp:
            if i in g_kp:
                elastic_net_selected.append(1)
            else:
                elastic_net_selected.append(0)
        p_df_subset['elastic_net_select'] = elastic_net_selected

        f_name = os.path.join(model_dir,'DEGs_with_elastic_info_GTEx+NAT_022123.csv')
        p_df_subset.to_csv(f_name,index = None)
        
        make_heatmap(df,g_kp)





