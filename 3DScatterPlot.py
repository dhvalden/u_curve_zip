import pandas as pd
import plotly.express as px
from sklearn.mixture import BayesianGaussianMixture
from sklearn.cluster import OPTICS
import matplotlib.pyplot as plt

dfw = pd.read_csv('./data/wilac_clustering_data.csv')
dfw.mean()
dfw


fig = px.scatter_3d(dfw, x='Gini2016', y='GDPpc2016', z='InsAcceptance',
                    color='quadraticEffect', hover_name='sample',
                    color_continuous_scale='RdBu',
                    color_continuous_midpoint=0, template="plotly_dark",
                    opacity=0.9,
                    title='Scatterplot Quadratic Effect on Willingness \
                    to Participate in Activism')
fig.write_html("./plots/wilacQua.html")

dfc = pd.read_csv('./data/colac_clustering_data.csv')
dfc.mean()

fig = px.scatter_3d(dfc, x='Gini2016', y='GDPpc2016', z='InsAcceptance',
                    color='quadraticEffect', hover_name='sample',
                    color_continuous_scale='RdBu',
                    color_continuous_midpoint=0, template="plotly_dark",
                    opacity=0.9,
                    title='Scatterplot Quadratic Effect on Participation \
                    in Activism')
fig.write_html("./plots/colacQua.html")

# Fit a Dirichlet process Gaussian mixture using five components
X = dfw[['Gini2016','GDPpc2016','InsAcceptance']].to_numpy()
bgmm = BayesianGaussianMixture(n_components=3,
                               covariance_type='full').fit(X)
labels = bgmm.predict(X)
plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')
plt.savefig('./plots/bgmm.png')


clustering = OPTICS(min_samples=4, xi=.05, min_cluster_size=.05).fit(X)
labels = clustering.labels_[clustering.ordering_]
plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis')
labels
plt.savefig('./plots/optics.png')
