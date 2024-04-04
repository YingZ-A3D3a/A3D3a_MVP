import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr,spearmanr
# input genes of interest (gene_oes), markov model output (final_prob_markov), and genes in the seed (g_seed)

def make_plots(gene_oes,final_prob_markov,g_seed):
    plt.rcParams['figure.figsize'] = [15, 10]
    fig, ax = plt.subplots(2,3)
    colors = ["b","g","r","c","m","y","k","tomato","gray","plum","pink","purple"] # 12 colors

    ##############################################################
    # plot 1

    x = final_prob_markov.degree_in_disease.values.tolist()
    y = final_prob_markov.final_probability.values.tolist()
    colors0 = []
    for i in final_prob_markov.genes.values.tolist():
        if i not in g_seed:
            colors0.append('lightgray')
        else:
            colors0.append('lightgreen')

    ax[0,0].scatter(x,y,color = colors0)

    ax[0,0].set_xlabel('degree in the disease network')
    ax[0,0].set_ylabel('final score')

    ##############################################################
    # plot 2

    x = final_prob_markov.degree_in_disease.values.tolist()
    y = final_prob_markov.final_probability.values.tolist()
    ax[0,1].scatter(x,y,color = 'lightsteelblue')
    for i in range(len(gene_oes)):
        indx_oe = final_prob_markov.loc[final_prob_markov.genes==gene_oes[i],].index.values[0]-1 # because the index of the final_prob_markov model starts from 1
        ax[0,1].scatter(x[indx_oe],y[indx_oe],color = colors[i])
        ax[0,1].text(x[indx_oe]-.0005, y[indx_oe]+.0005, gene_oes[i], fontsize=9)
    ax[0,1].set_xlabel('degree in the disease network')
    ax[0,1].set_ylabel('final score')

    ##############################################################
    # plot 3

    x = final_prob_markov.initial_score.values.tolist()
    y = final_prob_markov.final_probability.values.tolist()
    ax[0,2].scatter(x,y,color = 'lightsteelblue')
    for i in range(len(gene_oes)):
        indx_oe = final_prob_markov.loc[final_prob_markov.genes==gene_oes[i],].index.values[0]-1
        ax[0,2].scatter(x[indx_oe],y[indx_oe],color = colors[i])
        ax[0,2].text(x[indx_oe]-.0005, y[indx_oe]+.0005, gene_oes[i], fontsize=9)
    ax[0,2].set_xlabel('initial score')
    ax[0,2].set_ylabel('final score')

    ##############################################################
    # plot 4

    x = np.log2(final_prob_markov.publications.values.tolist())
    y = final_prob_markov.final_probability.values.tolist()
    # y = final_prob_markov.max_rank.values.tolist()
    ax[1,0].scatter(x,y,color = 'lightsteelblue')
    for i in range(len(gene_oes)):
        indx_oe = final_prob_markov.loc[final_prob_markov.genes==gene_oes[i],].index.values[0]-1
        ax[1,0].scatter(x[indx_oe],y[indx_oe],color = colors[i])
        ax[1,0].text(x[indx_oe]-.0005, y[indx_oe]+.0005, gene_oes[i], fontsize=9)

    # x1_0 = x.reshape(-1,1)
    # y = np.array(y)
    # y_0 = y.reshape(-1,1)
    # reg = LinearRegression().fit(x1_0,y_0)
    # reg_pred = reg.predict(x1_0)
    # # plt.plot(x1_0,reg_pred,color="black", linewidth=1)

    # ### mse,mae, and r2 (coefficient of determination) of regression
    # print(f"mae:{mean_absolute_error(y_0,reg_pred)}")
    # print(f"rmse:{np.sqrt(mean_squared_error(y_0,reg_pred))}")
    # print(f"r-squared:{r2_score(y_0,reg_pred)}")
    # print(f"coef:{reg.coef_}, intercept:{reg.intercept_}")

    # print(pearsonr(x,y))
    # print(spearmanr(x,y))

    ax[1,0].set_xlabel('log(publication)')
    ax[1,0].set_ylabel('final score')

    ##############################################################
    # plot 5
    x = final_prob_markov.degree_in_background.values.tolist()
    y = final_prob_markov.final_probability.values.tolist()
    # y = final_prob_markov.max_rank.values.tolist()
    ax[1,1].scatter(x,y,color = 'lightsteelblue')
    for i in range(len(gene_oes)):
        indx_oe = final_prob_markov.loc[final_prob_markov.genes==gene_oes[i],].index.values[0]-1
        ax[1,1].scatter(x[indx_oe],y[indx_oe],color = colors[i])
        ax[1,1].text(x[indx_oe]-.0005, y[indx_oe]+.0005, gene_oes[i], fontsize=9)

    ax[1,1].set_xlabel('degree in background interactome')
    ax[1,1].set_ylabel('final score')
    
    ##############################################################
    # plot 6

    x1 = final_prob_markov.initial_score.values.tolist()
    x = final_prob_markov.degree_in_disease.values.tolist()
    constant1 = 0.5
    constant2 = 0
    x1 = np.array(x)*(np.power(np.array(x1)+constant2,1)+constant1)

    ax[1,2].scatter(x1,y)

    x1_0 = x1.reshape(-1,1)
    y = np.array(y)
    y_0 = y.reshape(-1,1)
    # reg = LinearRegression().fit(x1_0,y_0)
    # reg_pred = reg.predict(x1_0)
    ax[1,2].set_xlabel(f'degree*[initial_score+{constant1}]in the disease network')
    ax[1,2].set_ylabel('final score')
    

    fig.tight_layout()

#     ##############################################################
#     # plot 7

#     x = np.log2(final_prob_markov.publications.values.tolist())
#     y = -np.log2(final_prob_markov.final_rank.values.tolist())
#     # y = final_prob_markov.max_rank.values.tolist()
#     ax[2,0].scatter(x,y,color = 'lightsteelblue')
#     for i in range(len(gene_oes)):
#         indx_oe = final_prob_markov.loc[final_prob_markov.genes==gene_oes[i],].index.values[0]-1
#         ax[2,0].scatter(x[indx_oe],y[indx_oe],color = colors[i])
#         ax[2,0].text(x[indx_oe]-.0005, y[indx_oe]+.0005, gene_oes[i], fontsize=9)

# #     x1_0 = x.reshape(-1,1)
# #     y = np.array(y)
# #     y_0 = y.reshape(-1,1)
# #     reg = LinearRegression().fit(x1_0,y_0)
# #     reg_pred = reg.predict(x1_0)
# #     # plt.plot(x1_0,reg_pred,color="black", linewidth=1)

# #     ### mse,mae, and r2 (coefficient of determination) of regression
# #     print(f"mae:{mean_absolute_error(y_0,reg_pred)}")
# #     print(f"rmse:{np.sqrt(mean_squared_error(y_0,reg_pred))}")
# #     print(f"r-squared:{r2_score(y_0,reg_pred)}")
# #     print(f"coef:{reg.coef_}, intercept:{reg.intercept_}")

# #     print(pearsonr(x,y))
# #     print(spearmanr(x,y))

#     ax[2,0].set_xlabel('log(publication)')
#     ax[2,0].set_ylabel('-log(final rank)')

#     ##############################################################
#     # plot 8
#     x = final_prob_markov.degree_in_background.values.tolist()
#     y = -np.log2(final_prob_markov.final_rank.values.tolist())
#     # y = final_prob_markov.max_rank.values.tolist()
#     ax[2,1].scatter(x,y,color = 'lightsteelblue')
#     for i in range(len(gene_oes)):
#         indx_oe = final_prob_markov.loc[final_prob_markov.genes==gene_oes[i],].index.values[0]-1
#         ax[2,1].scatter(x[indx_oe],y[indx_oe],color = colors[i])
#         ax[2,1].text(x[indx_oe]-.0005, y[indx_oe]+.0005, gene_oes[i], fontsize=9)

#     ax[2,1].set_xlabel('degree in background interactome')
#     ax[2,1].set_ylabel('-log(final rank)')

#     ##############################################################
#     # plot 9

#     x = final_prob_markov.initial_score.values.tolist()
#     y = -np.log2(final_prob_markov.final_rank.values.tolist())
#     ax[2,2].scatter(x,y,color = 'lightsteelblue')
#     for i in range(len(gene_oes)):
#         indx_oe = final_prob_markov.loc[final_prob_markov.genes==gene_oes[i],].index.values[0]-1
#         ax[2,2].scatter(x[indx_oe],y[indx_oe],color = colors[i])
#         ax[2,2].text(x[indx_oe]-.0005, y[indx_oe]+.0005, gene_oes[i], fontsize=9)
#     ax[2,2].set_xlabel('initial score')
#     ax[2,2].set_ylabel('-log(final rank)')