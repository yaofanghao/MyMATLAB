'''
Author: yao fanghao
Date: 2023-12-12 20:24:44
LastEditTime: 2023-12-12 20:32:31
LastEditors: yao fanghao
'''
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd

df = pd.read_excel('test_label.xlsx', header=None)
train_data = df.iloc[:, :8].values  # 前八列
train_labels = df.iloc[:, 8].values  # 第九列

# 创建PCA对象，并指定降维后的维度
pca = PCA(n_components=2)

# 使用PCA对数据进行降维
X_reduced = pca.fit_transform(train_data)

# 打印降维后的数据
print(X_reduced)

# 将降维后的数据转换为DataFrame
df = pd.DataFrame(data=X_reduced, columns=['PC1', 'PC2'])

# 导出到Excel文件
df.to_excel('PCA-test-result.xlsx', index=False)