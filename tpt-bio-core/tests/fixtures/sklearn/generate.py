"""Generates the scikit-learn reference values checked in tests/sklearn_reference.rs.

Requires numpy and scikit-learn (verified with scikit-learn 1.9.0). The datasets
must match the Rust test: X[i][j] = ((i*37 + j*11) % 17) + 0.5*j for the spread
data, and 10*(i//4) + ((i*7 + j*3) % 5)*0.1 for the blob data.
"""

import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

np.set_printoptions(precision=17, linewidth=200)

n, d = 12, 3
X = np.array([[((i * 37 + j * 11) % 17) + 0.5 * j for j in range(d)] for i in range(n)], dtype=float)
B = np.array([[10.0 * (i // 4) + ((i * 7 + j * 3) % 5) * 0.1 for j in range(d)] for i in range(n)])

z = StandardScaler().fit_transform(X)
print("zscore row 0:", z[0].tolist(), "row 11:", z[-1].tolist())

med = np.median(X, axis=0)
mad = np.median(np.abs(X - med), axis=0)
rob = (X - med) / (1.4826 * mad)
print("robust row 0:", rob[0].tolist(), "row 11:", rob[-1].tolist())

pca = PCA(n_components=2).fit(X)
print("pca mean:", pca.mean_.tolist())
print("pca explained variance:", pca.explained_variance_.tolist())
print("pca explained variance ratio:", pca.explained_variance_ratio_.tolist())
print("pca components:", pca.components_.tolist())
print("pca scores row 0:", pca.transform(X)[0].tolist())

km = KMeans(n_clusters=3, n_init=10, random_state=0).fit(B)
print("kmeans inertia:", km.inertia_)
print("kmeans labels:", km.labels_.tolist())
